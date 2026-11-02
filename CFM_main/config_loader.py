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
from pathlib import Path


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
