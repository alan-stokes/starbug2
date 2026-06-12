import os
import pytest
from starbug2.constants import STAR_BUG_PARAMS
from starbug2.star_bug_config import StarBugMainConfig


def test_parse_param():
    assert type(StarBugMainConfig.parse_param("A=1//.")) == dict

    assert StarBugMainConfig.parse_param("A = 1 //.") == {'A': 1}
    assert StarBugMainConfig.parse_param("A = B //.") == {'A': 'B'}
    assert StarBugMainConfig.parse_param("A = B //.\n") == {'A': 'B'}

    assert StarBugMainConfig.parse_param("A = //.\n") == {'A': ''}
    assert StarBugMainConfig.parse_param(" = //.") == {}

    assert StarBugMainConfig.parse_param("A=B") == {"A": "B"}
    assert StarBugMainConfig.parse_param("A=B/") == {"A": "B/"}
    assert StarBugMainConfig.parse_param("A=B/.") == {"A": "B/."}
    assert StarBugMainConfig.parse_param("A=1/.") == {"A": "1/."}

    assert StarBugMainConfig.parse_param("A      =1")=={"A": 1}
    assert StarBugMainConfig.parse_param("A=1      ")=={"A": 1}
    assert StarBugMainConfig.parse_param("A=1     a")=={"A": "1     a"}

def test_load_default_params():
    config: StarBugMainConfig = StarBugMainConfig()
    assert config.param_tag == STAR_BUG_PARAMS


def test_load_params():
    config: StarBugMainConfig = StarBugMainConfig.load_params("does_not_exist")

    os.system("starbug2 --local-param")
    second_config: StarBugMainConfig = (
        StarBugMainConfig.load_params("starbug.param"))

    for value, _ in config.MAIN_PARAM_FILE_MAP.values():
        assert getattr(config, value) == getattr(second_config, value)
    assert second_config.param_tag == STAR_BUG_PARAMS
    assert config.param_tag == STAR_BUG_PARAMS
    os.remove("starbug.param")

def test_update_params():
    os.system("starbug2 --local-param")
    os.system("sed -i s/PARAM/PARAM1/g starbug.param")

    with pytest.raises(
            TypeError,
        match="Param PARAM1 no longer works within Starbug2. Please "
              "execute starbug2 --update-param"):
        StarBugMainConfig.load_params("starbug.param")
    os.remove("starbug.param")

