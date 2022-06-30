#!/usr/bin/env python3
import os
import sys
import subprocess
def get_virtualenv_path():
    """Used to work out path to install compiled binaries to."""
    if hasattr(sys, 'real_prefix'):
        return sys.prefix + "/bin/"

    if hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix:
        return sys.prefix + "/bin/"

    if 'conda' in sys.prefix:
        return sys.prefix + "/bin/"

    return "/usr/local/bin/"
def copy_and_move_software():
    """Used the subprocess module to compile/install the C software."""
    flash_path = os.path.dirname(os.path.abspath("./"))
    flash_path = os.path.join(flash_path , "umierrorcorrect_forensics/umierrorcorrect_forensics/FLASH-lowercase-overhang/")
    venv = get_virtualenv_path()
    # compile the software
    subprocess.check_call("make", cwd=flash_path)
    cmd = "cp -R {} {}".format(os.path.join(flash_path,"flash"),venv)
    print(flash_path)
    print(venv)
    print(cmd)
    subprocess.check_call(cmd, shell=True)

copy_and_move_software()
