import yaml

from pan2met.conf.schema import validate


def test_default_configuration_validation():
    with open("src/pan2met/conf/default.yaml") as f:
        default_config = yaml.safe_load(f)
    assert validate(default_config)
