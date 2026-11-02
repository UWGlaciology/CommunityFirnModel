'''
Tests of config_loader.load_config: .json and .toml configs give the same
dictionary, and a model built from a .toml config is the same as one built
from the .json config.
'''

import json

import pytest

from cfm_testing import CFM_MAIN, ShortRun
from config_loader import load_config


def to_toml(c):
    '''Write a flat config dictionary as TOML text (JSON strings are valid TOML strings).'''
    def value(v):
        if isinstance(v, bool):
            return 'true' if v else 'false'
        if isinstance(v, list):
            return '[' + ', '.join(value(x) for x in v) + ']'
        return json.dumps(v)
    return ''.join(f'{k} = {value(v)}\n' for k, v in c.items())


def test_toml_with_comments_is_read(tmp_path):
    path = tmp_path / 'c.toml'
    path.write_text('''# a comment line
liquid = "bucket"      # trailing comment
RhoImp = 830           # kg m-3
ThickImp = 0.1
LeakyLenses = false
outputs = ["density", "depth"]   # list
''')
    assert load_config(path) == {'liquid': 'bucket', 'RhoImp': 830, 'ThickImp': 0.1,
                                 'LeakyLenses': False, 'outputs': ['density', 'depth']}


def test_json_is_the_default(tmp_path):
    c = {'liquid': 'bucket', 'RhoImp': 830, 'flag': True, 'outputs': ['density']}
    for name in ('c.json', 'c.JSON', 'c.cfg', 'c'):
        (tmp_path / name).write_text(json.dumps(c))
        assert load_config(tmp_path / name) == c
    assert load_config(str(tmp_path / 'c.json')) == c  # a string path works too


def test_bad_toml_raises(tmp_path):
    path = tmp_path / 'c.toml'
    path.write_text('liquid = "bucket\n')
    with pytest.raises(Exception, match='.'):
        load_config(path)


@pytest.mark.parametrize('name', ['example.json', 'example_df.json', 'example_csv.json', 'AirConfig.json'])
def test_example_configs_survive_a_toml_round_trip(tmp_path, name):
    c = json.loads((CFM_MAIN / name).read_text())
    path = tmp_path / 'c.toml'
    path.write_text(to_toml(c))
    assert load_config(path) == c


def test_model_from_toml_config_equals_model_from_json_config(tmp_path):
    run = ShortRun(tmp_path / 'run', years=1)
    toml_path = tmp_path / 'run' / 'config.toml'
    toml_path.write_text(to_toml(run.config))
    model_json = run.build()
    run.config_path = toml_path
    model_toml = run.build()
    assert model_toml.c == model_json.c


def test_example_toml_matches_example_json():
    '''example.toml is the commented copy of example.json: same keys and values, except that
    the *_options helper arrays of the json are comments in the toml.'''
    from config_validation import _is_helper_key
    j = json.loads((CFM_MAIN / 'example.json').read_text())
    j = {k: v for k, v in j.items() if not _is_helper_key(k)}
    assert load_config(CFM_MAIN / 'example.toml') == j


def test_example_toml_option_comments_list_the_json_options():
    '''Every "# options:" comment in example.toml lists the same values as the json's helper array.'''
    j = json.loads((CFM_MAIN / 'example.json').read_text())
    lines = (CFM_MAIN / 'example.toml').read_text().splitlines()
    commented = {}
    for line, nxt in zip(lines, lines[1:]):
        if line.startswith('# options:'):
            commented[nxt.split(' = ')[0]] = [x.strip() for x in line[len('# options:'):].split(',')]
    for key, options in commented.items():
        helper = j.get(key + '_options', j.get('isoOptions') if key == 'iso' else None)
        assert options == [str(x) for x in helper], key
    assert {k[:-len('_options')] for k in j if k.endswith('_options')} <= set(commented)
