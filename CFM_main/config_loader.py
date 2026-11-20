'''
config_loader.py
================
Reads CFM configuration files. The format is chosen from the file extension:
``.toml`` is read with the standard-library ``tomllib`` (Python 3.11+); any
other extension (normally ``.json``) is read as JSON, so existing configs keep
working. TOML allows comments, which JSON does not.

Both formats give the same flat dictionary of strings, numbers, booleans and
lists. (TOML has no null, so a config that needs ``null`` has to stay JSON.)
'''

import json
import subprocess
import tomllib
import zlib
from pathlib import Path

import numpy as np

from _version import __version__


def load_config(path):
    '''
    Read a CFM configuration file and return it as a dictionary.

    :param path: path to a ``.toml`` or ``.json`` file. Any extension other than
        ``.toml`` is parsed as JSON.
    '''
    path = Path(path)
    text = path.read_text(encoding='utf-8')
    if path.suffix.lower() == '.toml':
        return tomllib.loads(text)
    return json.loads(text)


def provenance():
    '''
    Return the version and git revision of the CFM source as a dictionary
    (``cfm_version``, ``cfm_git``), to be stored in the configuration of a run.

    ``cfm_git`` is the output of ``git describe --tags --dirty --always`` in the
    directory of this file, so ``-dirty`` flags uncommitted changes. It is
    ``'unknown'`` when the source is not a git checkout or git is unavailable.
    '''
    try:
        rev = subprocess.run(['git', 'describe', '--tags', '--dirty', '--always'],
                             cwd=Path(__file__).resolve().parent, capture_output=True,
                             text=True, timeout=5, check=True).stdout.strip() or 'unknown'
    except Exception:
        rev = 'unknown'
    return {'cfm_version': __version__, 'cfm_git': rev}


def get_rng(c, stream=''):
    '''
    Return a numpy random Generator for one source of randomness in a run.

    The seed is the configuration key ``random_seed`` (default 0), so a run is
    repeated exactly by the configuration file alone. ``stream`` names the source
    (for example ``'srho'`` or an isotope name); each name gets its own
    independent sequence, so the draws of one source do not depend on whether
    another source is used.

    :param c: configuration dictionary
    :param stream: name of the source of randomness
    '''
    seed = int(c.get('random_seed', 0))
    return np.random.default_rng([seed, zlib.crc32(str(stream).encode())])
