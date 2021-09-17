import os
import yaml
from jinja2 import Template

required_keys = ["flux-product-dir", "project-tmpdir"]

def get_config_dict():
    """return the configuration dictionary with environment variables replaced"""
    
    with open("_config-calc.yml") as fid:
        config_dict_in = yaml.safe_load(fid)    
    
    for key in required_keys:
        assert key in config_dict_in, f"config missing {key}"
        
    config_dict = {}
    for key, value in config_dict_in.items():
        if isinstance(value, str):
            t = Template(value)
            config_dict[key] = t.render(env=os.environ)
        else:
            config_dict[key] = value
    return config_dict


# get configuration dictionary
config_dict = get_config_dict()

# cache directory for big file
flux_product_dir = config_dict["flux-product-dir"]
os.makedirs(flux_product_dir, exist_ok=True)

project_tmpdir = config_dict["project-tmpdir"]
os.environ["INTAKE_LOCAL_CACHE_DIR"] = f"{project_tmpdir}/intake-cache"