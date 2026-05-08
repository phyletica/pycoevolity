#!/usr/bin/env python

import os
import subprocess
import shutil

def get_exe_dir(exe_name, dir_to_check = None):
    exe_path = None
    if not dir_to_check:
        exe_path = shutil.which(exe_name)
        if not exe_path:
            raise Exception(
                f"'{exe_name}' was not found in system's PATH"
            )
    else:
        exe_path = os.path.join(dir_to_check, exe_name)

    if not os.access(exe_path, os.X_OK):
        raise Exception(
            f"{exe_name} found at '{exe_path}', but does have execute "
            f"permissions"
        )
    return os.path.dirname(exe_path)

def run_cmd(cmd, timeout = None, cwd = None):
    result = subprocess.run(
        cmd,
        capture_output = True,
        text = True,
        check = True,
        timeout = timeout,
        cwd = cwd,
    )
    return result
