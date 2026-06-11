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
from starbug2 import param
from starbug2.constants import PARAM, STAR_BUG_PARAMS


def test_parse_param():
    assert type(param.parse_param("A=1//.")) == dict

    assert param.parse_param("A = 1 //.") == {'A': 1}
    assert param.parse_param("A = B //.") == {'A': 'B'}
    assert param.parse_param("A = B //.\n") == {'A': 'B'}

    assert param.parse_param("A = //.\n") == {'A': ''}
    assert param.parse_param(" = //.") == {}

    assert param.parse_param("A=B") == {"A": "B"}
    assert param.parse_param("A=B/") == {"A": "B/"}
    assert param.parse_param("A=B/.") == {"A": "B/."}
    assert param.parse_param("A=1/.") == {"A": "1/."}

    assert param.parse_param("A      =1")=={"A": 1}
    assert param.parse_param("A=1      ")=={"A": 1}
    assert param.parse_param("A=1     a")=={"A": "1     a"}

def test_load_default_params():
    assert param.load_default_params()!={}
    assert type(param.load_default_params()) == dict
    assert PARAM in param.load_default_params().keys()
    assert param.load_default_params().get(PARAM) == STAR_BUG_PARAMS


def test_load_params():
    assert param.load_default_params() == param.load_params(None)

    assert param.load_params("doesnotexist") == {}

    os.system("starbug2 --local-param")
    assert param.load_params("starbug.param") != {}
    assert PARAM in param.load_params("starbug.param").keys()
    assert (
        param.load_params("starbug.param").get(PARAM) == STAR_BUG_PARAMS)
    os.remove("starbug.param")

def test_update_params():
    os.system("starbug2 --local-param")
    os.system("sed -i s/PARAM/PARAM1/g starbug.param")

    assert "PARAM" not in param.load_params("starbug.param").keys()
    assert param.update_param_file("starbug.param") is None
    assert param.update_param_file("starbug.param") is None
    os.remove("starbug.param")

