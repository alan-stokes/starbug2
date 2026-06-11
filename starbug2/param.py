import os
from typing import Dict
from starbug2.star_bug_config import StarBugMainConfig
from starbug2.utils import printf,p_error,get_version


def _load_params_old(f_name: str | None) -> Dict[str, int | float| str]:
    """
    Convert a parameter file into a dictionary of options

    :param f_name: path/to/file.param
    :type f_name: str or None
    :return: dictionary of options
    :rtype: dict of string, string
    """
    if f_name is None:
        return {}

    config: Dict[str, int | float| str] = {}
    if os.path.exists(f_name):
        with open(f_name, "r") as fp:
            for line in fp.readlines():
                config.update(StarBugMainConfig.parse_param(line))
    else:
        p_error("config file \"%s\" does not exist\n" % f_name)
    return config

def update_param_file(f_name: str | None) -> None:
    """
    When the local parameter file is from an older version, add or remove the
    new or obsolete keys

    :param f_name: local file to update
    :type f_name: str
    :return: None
    """
    if f_name is None:
        return

    default_param = StarBugMainConfig()

    current_param = _load_params_old(f_name)

    if os.path.exists(f_name):
        printf("Updating \"%s\"\n" % f_name)
        fpo = open("/tmp/starbug.param",'w')

        add_keys = (
            set(default_param.MAIN_PARAM_FILE_MAP.keys()) -
            set(current_param.keys()))
        del_keys = (
            set(current_param.keys()) -
            set(default_param.MAIN_PARAM_FILE_MAP.keys()))
        if add_keys:
            printf("-> adding: %s  \n"%(', '.join(add_keys)))
        if del_keys:
            printf("-> removing: %s\n"%(', '.join(del_keys)))
        
        if not len(add_keys | del_keys):
            printf("-> No updates needed\n")

        # remove invalid params
        for key in del_keys:
            current_param.pop(key, None)

        # update config object with their settings
        default_param.update(current_param)

        # generate the new output and write to file.
        output: str = default_param.generate_default_param_file_text(
            get_version())
        fpo.write("%s\n" % output)
        fpo.close()
        os.system("mv /tmp/starbug.param %s" % f_name)
    else:
        p_error("local parameter file '%s' does not exist\n" % f_name)

