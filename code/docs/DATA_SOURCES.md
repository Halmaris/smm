# Data sources and input format

Data are downloaded separately under each provider's terms; no dataset or result
is distributed with this repository. The links below are the release references
used in the manuscript. Temporal tasks use the TSML/UEA and Multiverse releases
on Zenodo; image sources are MedMNIST, UCI, LIBSVM and Makerere Beans.

Supply `X` as a list of equally sized numeric matrices (or a `p × q × n` array)
and `y` as class labels. Keep the released TEST split separate. MedMNIST and Beans
TRAIN and VALIDATION form the development set. Beans is bilinearly resized to
28 × 28 without augmentation. RGB images stack their three planes into an
84 × 28 matrix with `structure = 'spatial', channels = 3`; channel graphs remain
separate. Grayscale images use one spatial graph. Temporal matrices have channels
in rows and ordered measurements in columns.

For tuning, supply fixed stratified fold identifiers 1–5. Scaling is estimated
within each training fold, never using TEST. A prepared simulation block is a
list containing `train` (`X`, `y`, `folds`) and `test` (`X`, `y`); see the generator
functions for the exact object structure. Simulation inputs are generated locally
from explicit seeds; they are not external datasets.

## Release links

| Dataset | Source |
|---|---|
| AppliancesEnergy_disc | [release](https://zenodo.org/records/18497240) |
| ArticularyWordRecognition | [release](https://zenodo.org/records/11204924) |
| AsphaltObstaclesCoordinates | [release](https://zenodo.org/records/18497270) |
| AsphaltPavementTypeCoordinates | [release](https://zenodo.org/records/18497280) |
| AsphaltRegularityCoordinates | [release](https://zenodo.org/records/18497298) |
| AtrialFibrillation | [release](https://zenodo.org/records/11206175) |
| AutomotiveRoadTrials | [release](https://zenodo.org/records/18497336) |
| BasicMotions | [release](https://zenodo.org/records/11206179) |
| Beans | [release](https://huggingface.co/datasets/AI-Lab-Makerere/beans) |
| Blink | [release](https://zenodo.org/records/18735667) |
| BreastMNIST | [release](https://zenodo.org/records/10519652) |
| ButtonPress | [release](https://zenodo.org/records/15522930) |
| CharacterTrajectories | [release](https://zenodo.org/records/18723007) |
| Cricket | [release](https://zenodo.org/records/11206185) |
| EmoPain | [release](https://zenodo.org/records/18499143) |
| Epilepsy | [release](https://zenodo.org/records/11206204) |
| ERing | [release](https://zenodo.org/records/11206210) |
| EyesOpenShut | [release](https://zenodo.org/records/18496872) |
| FingerMovements | [release](https://zenodo.org/records/11206220) |
| HandMovementDirection | [release](https://zenodo.org/records/11206224) |
| Heartbeat | [release](https://zenodo.org/records/11206229) |
| IRDS-EFL | [release](https://zenodo.org/records/18510498) |
| IRDS-EFR | [release](https://zenodo.org/records/18510501) |
| IRDS-SAL | [release](https://zenodo.org/records/18510505) |
| IRDS-SAR | [release](https://zenodo.org/records/18510507) |
| IRDS-SFE | [release](https://zenodo.org/records/18510509) |
| IRDS-SFR | [release](https://zenodo.org/records/18510518) |
| IRDS-STL | [release](https://zenodo.org/records/18510520) |
| IRDS-STR | [release](https://zenodo.org/records/18510522) |
| IRDS-SFL | [release](https://zenodo.org/records/18510515) |
| JapaneseVowels | [release](https://zenodo.org/records/18735628) |
| KERAAL-CTK | [release](https://zenodo.org/records/18510796) |
| KERAAL-ELK | [release](https://zenodo.org/records/18510806) |
| KERAAL-RTK | [release](https://zenodo.org/records/18510813) |
| KIMORE-LA-C | [release](https://zenodo.org/records/18510879) |
| KIMORE-LT-C | [release](https://zenodo.org/records/18510884) |
| KIMORE-PR-C | [release](https://zenodo.org/records/18510887) |
| KIMORE-Sq-C | [release](https://zenodo.org/records/18510889) |
| KIMORE-TR-C | [release](https://zenodo.org/records/18510893) |
| KINECAL-3WFV | [release](https://zenodo.org/records/18510895) |
| KINECAL-GGFV | [release](https://zenodo.org/records/18510901) |
| Libras | [release](https://zenodo.org/records/11206239) |
| LowCost | [release](https://zenodo.org/records/15523038) |
| NATOPS | [release](https://zenodo.org/records/11206248) |
| Optdigits | [release](https://archive.ics.uci.edu/dataset/80/optical+recognition+of+handwritten+digits) |
| PneumoniaMNIST | [release](https://zenodo.org/records/10519652) |
| RacketSports | [release](https://zenodo.org/records/11206263) |
| RetinaMNIST | [release](https://zenodo.org/records/10519652) |
| SelfRegulationSCP1 | [release](https://zenodo.org/records/11206265) |
| SelfRegulationSCP2 | [release](https://zenodo.org/records/11206269) |
| SPHERE-WUS | [release](https://zenodo.org/records/18510936) |
| UCDHE-MP | [release](https://zenodo.org/records/18510938) |
| UCDHE-MP-MC | [release](https://zenodo.org/records/18510940) |
| UCDHE-Rowing | [release](https://zenodo.org/records/18510942) |
| UCDHE-Rowing-MC | [release](https://zenodo.org/records/18510948) |
| UIPRMD-DS-C | [release](https://zenodo.org/records/18511325) |
| UIPRMD-HS-C | [release](https://zenodo.org/records/18511341) |
| UIPRMD-IL-C | [release](https://zenodo.org/records/18511343) |
| UIPRMD-SASLR-C | [release](https://zenodo.org/records/18511347) |
| UIPRMD-SL-C | [release](https://zenodo.org/records/18511349) |
| UIPRMD-SSA-C | [release](https://zenodo.org/records/18511351) |
| UIPRMD-SSE-C | [release](https://zenodo.org/records/18511355) |
| UIPRMD-SSIER-C | [release](https://zenodo.org/records/18511357) |
| UIPRMD-SSS-C | [release](https://zenodo.org/records/18511361) |
| UIPRMD-STS-C | [release](https://zenodo.org/records/18511363) |
| USPS | [release](https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/multiclass.html#usps) |
| UWaveGestureLibrary | [release](https://zenodo.org/records/11206282) |
