"""In-memory Optuna search with the study's deterministic CV upper-bound pruning."""
import itertools
import warnings
import optuna


# Pruning and fixed anchor configurations --------------------------------------


def complete_trials(study):
    """Return fully evaluated trials in this Optuna study, in proposal order."""
    return study.get_trials(deepcopy=False, states=(optuna.trial.TrialState.COMPLETE,))


class BoundPruner(optuna.pruners.BasePruner):
    def prune(self, study, trial):
        """Prune trial only when its stored BA upper bound cannot match study's incumbent."""
        if trial.user_attrs.get('anchor', False) or trial.user_attrs.get('fold', 5) >= 5:
            return False
        complete = complete_trials(study)
        upper = trial.user_attrs.get('upper')
        return upper is not None and bool(complete) and upper < max(t.value for t in complete) - 1e-12


def anchors(method):
    """Return sorted C/tau/gamma anchor dictionaries for ablation method M00--M11."""
    if method not in ('M00', 'M10', 'M01', 'M11'):
        raise ValueError('Unknown ablation method')
    return [dict(Cbar=c, tau=tau, gamma=gamma) for c, tau, gamma in itertools.product(
        (.1, 1., 10., 100.), (0., .1, 1.) if method[1] == '1' else (0.,),
        (0., .1, 1.) if method[2] == '1' else (0.,))]


# Sequential search used by the R interface ------------------------------------


class Search:
    def __init__(self, method, seed, budget=100, gamma_max=10000.):
        """Create sequential TPE; method sets penalties, seed RNG, budget proposals, gamma_max domain."""
        self.method = method
        self.budget = int(budget)
        self.gamma_max = float(gamma_max)
        initial = anchors(method)
        limits = {'tau': (1e-4, 100.), 'gamma': (1e-5, self.gamma_max)}
        self.domains = {name: bounds for flag, (name, bounds) in zip(method[1:], limits.items())
                        if flag == '1'}
        if self.budget < len(initial) or self.gamma_max < 1:
            raise ValueError('Budget/domain must contain the fixed anchors')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', optuna.exceptions.ExperimentalWarning)
            sampler = optuna.samplers.TPESampler(seed=int(seed), n_startup_trials=20,
                multivariate=True, group=True, constant_liar=False)
        self.study = optuna.create_study(direction='maximize', sampler=sampler, pruner=BoundPruner())
        self.anchor_count = len(initial)
        for row in initial:
            encoded = {'Cbar': row['Cbar']}
            for name in self.domains:
                encoded[name + '_active'] = row[name] > 0
                if row[name] > 0:
                    encoded[name + '_positive'] = row[name]
            self.study.enqueue_trial(encoded)
        self.active = None

    def ask(self):
        """Request the next proposal; return None when the requested budget is exhausted."""
        if self.active is not None:
            raise RuntimeError('Finish the active proposal first')
        if sum(t.state.is_finished() for t in self.study.trials) >= self.budget:
            return None
        trial = self.study.ask()
        pars = {'Cbar': trial.suggest_float('Cbar', 1e-2, 1e3, log=True), 'tau': 0., 'gamma': 0.}
        for name, (lo, hi) in self.domains.items():
            if trial.suggest_categorical(name + '_active', [False, True]):
                pars[name] = trial.suggest_float(name + '_positive', lo, hi, log=True)
        trial.set_user_attr('parameters', pars)
        anchor = trial.number < self.anchor_count
        trial.set_user_attr('anchor', anchor)
        self.active = trial
        return {'trial': trial.number + 1, 'parameters': pars, 'anchor': anchor}

    def report(self, score, upper, fold):
        """Report partial pooled BA score, optimistic final BA upper, and completed fold (1--5)."""
        self.active.report(float(score), int(fold))
        self.active.set_user_attr('upper', float(upper))
        self.active.set_user_attr('fold', int(fold))
        return self.active.should_prune()

    def tell(self, score=None, pruned=False):
        """Finish the current proposal with full-CV score or a certified pruned status."""
        if pruned:
            if not self.active.should_prune():
                raise ValueError('Pruning is not justified by the recorded bound')
            self.study.tell(self.active, state=optuna.trial.TrialState.PRUNED)
        else:
            self.study.tell(self.active, float(score))
        self.active = None

    def best(self):
        """Return the complete trial selected by pooled BA then lexicographic parameter tie-break."""
        complete = complete_trials(self.study)
        maximum = max(t.value for t in complete)
        selected = min((t for t in complete if t.value >= maximum - 1e-12),
            key=lambda t: tuple(t.user_attrs['parameters'][p] for p in ('Cbar', 'tau', 'gamma')) + (t.number,))
        return dict(trial=selected.number + 1, score=selected.value, parameters=selected.user_attrs['parameters'])
