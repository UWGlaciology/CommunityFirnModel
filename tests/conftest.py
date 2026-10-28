'''
pytest configuration for the Community Firn Model.

The CFM_main modules import each other by plain module name (e.g.
``from constants import *``), so CFM_main is put on sys.path here and the
tests import modules the same way the model does.

Run the tests from the repository root with:

    python -m pytest tests
'''

import sys
from pathlib import Path

CFM_MAIN = Path(__file__).resolve().parents[1] / 'CFM_main'
sys.path.insert(0, str(CFM_MAIN))
