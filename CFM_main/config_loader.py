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
import tomllib
import zlib
from pathlib import Path

import numpy as np


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
