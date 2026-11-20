"""Verify local input files against a historical configuration's SHA-256 hashes."""
import argparse
import hashlib
import json
from pathlib import Path


def verify(config, root):
    """Return mismatches; config is a JSON path and root holds the input files."""
    hashes = json.loads(Path(config).read_text())['data_sha256']
    failures = []
    for name, expected in hashes.items():
        path = Path(root) / name
        if not path.is_file():
            failures.append(f'{name}: missing')
            continue
        with path.open('rb') as stream:
            digest = hashlib.file_digest(stream, 'sha256').hexdigest()
        if digest != expected:
            failures.append(f'{name}: checksum mismatch')
    return failures


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', type=Path)
    parser.add_argument('inputs', type=Path)
    args = parser.parse_args()
    errors = verify(args.config, args.inputs)
    print('\n'.join(errors) if errors else 'All input checksums match.')
    raise SystemExit(bool(errors))
