"""Copyright (C) 2026 UKATC

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>."""

import os
import sys
from pathlib import Path
from typing import Final
import subprocess

import pytest
from starbug2.constants import STAR_BUG_PARAMS, PROBLEMATIC_FILTER_WARNING
from starbug2.core.star_bug_config import StarBugMainConfig
from tests.generic import TEST_PATH_STR, verify_test_data_exists, TEST_PATH

TEST_PARAM_PATH: Final[str] = os.path.join(
    TEST_PATH_STR, "../param_files/old_format.param")
starbug2_executable = Path(sys.executable).parent / "starbug2"


def test_parse_param():
    verify_test_data_exists()
    assert isinstance(StarBugMainConfig.parse_param("A=1//."), dict)

    assert StarBugMainConfig.parse_param("A = 1 //.") == {'A': 1}
    assert StarBugMainConfig.parse_param("A = B //.") == {'A': 'B'}
    assert StarBugMainConfig.parse_param("A = B //.\n") == {'A': 'B'}

    assert StarBugMainConfig.parse_param("A = //.\n") == {'A': ''}
    assert StarBugMainConfig.parse_param(" = //.") == {}

    assert StarBugMainConfig.parse_param("A=B") == {"A": "B"}
    assert StarBugMainConfig.parse_param("A=B/") == {"A": "B/"}
    assert StarBugMainConfig.parse_param("A=B/.") == {"A": "B/."}
    assert StarBugMainConfig.parse_param("A=1/.") == {"A": "1/."}

    assert StarBugMainConfig.parse_param("A      =1") == {"A": 1}
    assert StarBugMainConfig.parse_param("A=1      ") == {"A": 1}
    assert StarBugMainConfig.parse_param("A=1     a") == {"A": "1     a"}


def test_load_default_params():
    config: StarBugMainConfig = StarBugMainConfig()
    assert config.param_tag == STAR_BUG_PARAMS


def test_load_params():
    config: StarBugMainConfig = StarBugMainConfig.load_params("does_not_exist")

    subprocess.run(
        [str(starbug2_executable), "--local-param",
         f"--output={TEST_PATH}"], check=True, text=True,
        capture_output=True,)
    second_config: StarBugMainConfig = (
        StarBugMainConfig.load_params("starbug.param"))

    for value, _ in config.MAIN_PARAM_FILE_MAP.values():
        assert getattr(config, value) == getattr(second_config, value)
    assert second_config.param_tag == STAR_BUG_PARAMS
    assert config.param_tag == STAR_BUG_PARAMS

    os.remove(os.path.join(TEST_PATH_STR, "starbug.param"))


def test_update_params():
    subprocess.run(
        [str(starbug2_executable), "--local-param", f"--output={TEST_PATH}"],
        check=True)

    # change param to param1 to force a fail
    param_path = os.path.join(TEST_PATH_STR, 'starbug.param')
    content = Path(param_path).read_text()
    Path(param_path).write_text(content.replace("PARAM", "PARAM1"))

    with pytest.raises(
        TypeError,
        match="Param PARAM1 no longer works within Starbug2. Please "
              "execute starbug2 --update-param"):
        StarBugMainConfig.load_params(
            os.path.join(TEST_PATH_STR, 'starbug.param'))
    os.remove(os.path.join(TEST_PATH_STR, 'starbug.param'))


def test_update_params_no_changes():
    subprocess.run(
        [str(starbug2_executable), "--local-param", f"--output={TEST_PATH}"],
        check=True)
    subprocess.run(
        [str(starbug2_executable), "--update-param", f"--output={TEST_PATH}"],
        check=True)
    os.remove(os.path.join(TEST_PATH_STR, 'starbug.param'))


def test_update_params_old_to_new():
    subprocess.run(
        [str(starbug2_executable), f"-p {TEST_PARAM_PATH}", "--update-param"])


def test_f150w2_filter(capsys):
    config: StarBugMainConfig = StarBugMainConfig()
    config.custom_filter = "F150W2"

    # check for warning
    captured = capsys.readouterr()
    assert PROBLEMATIC_FILTER_WARNING in captured.err
