"""Running the tool with "python -m brucellaid"."""

import importlib
import runpy
import subprocess
import sys

import pytest

from brucellaid import __version__


def test_python_m_brucellaid():
    out = subprocess.run([sys.executable, '-m', 'brucellaid', '--version'], capture_output=True, text=True)
    assert out.returncode == 0
    assert out.stdout.strip() == 'brucellaid {}'.format(__version__)


def test_the_main_module_runs_the_command(monkeypatch, capsys):
    monkeypatch.setattr(sys, 'argv', ['brucellaid', '--version'])
    with pytest.raises(SystemExit) as exit_info:
        runpy.run_module('brucellaid', run_name='__main__')
    assert exit_info.value.code == 0
    assert capsys.readouterr().out.strip() == 'brucellaid {}'.format(__version__)


def test_importing_the_main_module_runs_nothing(capsys):
    """Importing it must not run the command: only "python -m brucellaid" does."""
    module = importlib.import_module('brucellaid.__main__')
    importlib.reload(module)
    assert module.main is not None
    assert capsys.readouterr().out == ''
