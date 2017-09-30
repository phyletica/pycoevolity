#! /usr/bin/env python

import os
import sys
import subprocess
import unittest
import logging

from pycoevolity.stringutils import random_str
from pycoevolity.tempfs import TempFileSystem
from pycoevolity.test.support import package_paths

_LOG = logging.getLogger(__name__)

class PycoevolityTestCase(unittest.TestCase):
    
    def set_up(self):
        self.temp_fs = TempFileSystem(
                parent = package_paths.test_path(),
                prefix = 'PycoevolityTestTemp-')
        self.test_id = 'pycoevolity-' + random_str()

    def tear_down(self):
        self.register_file_system()
        self.temp_fs.purge()

    def get_test_path(self, parent=None, prefix='temp'):
        return self.temp_fs.get_file_path(parent=parent, prefix=prefix)

    def get_test_subdir(self, parent=None, prefix='temp'):
        return self.temp_fs.create_subdir(parent=parent, prefix=prefix)

    def register_file(self, path):
        self.temp_fs._register_file(path)

    def register_dir(self, path):
        self.temp_fs._register_dir(path)

    def register_file_system(self):
        _LOG.debug('registering test file system...')
        for path, dirs, files, in os.walk(self.temp_fs.base_dir):
            for f in files:
                if f.startswith(self.test_id):
                    self.register_file(os.path.join(path, f))
            for d in dirs:
                if d.startswith(self.test_id):
                    self.register_dir(os.path.join(path, d))

    def exe_script(self, script_name, args, stdout = None, stderr = None,
            return_code = 0):
        script_path = package_paths.script_path(script_name)
        if isinstance(args, str):
            arg_list = args.split()
        else:
            arg_list = args
        arg_list = [str(x) for x in arg_list]
        cmd = [sys.executable, script_path] + arg_list
        _LOG.debug('Invocation:\n\t{0}'.format(' '.join(cmd)))
        p = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
        o, e  = p.communicate()
        exit_code = p.wait()
        if exit_code != return_code:
            _LOG.error("exit code {0} did not match {1}".format(exit_code,
                    return_code))
            _LOG.error("here is the stdout:\n{0}".format(o))
            _LOG.error("here is the stderr:\n{0}".format(e))
        self.assertEqual(exit_code, return_code)
        if stdout != None:
            if o != stdout:
                _LOG.error("std out did not match expected:\n{0}".format(o))
            self.assertEqual(o, stdout)
        if stderr != None:
            if e != stderr:
                _LOG.error("std error did not match expected:\n{0}".format(e))
            self.assertEqual(e, stderr)

    def assertApproxEqual(self, x, y, percent_tol=1e-6):
        eq = (((abs(x-y) / ((abs(x)+abs(y))/2.0))*100) < percent_tol)
        if not eq:
            _LOG.error('x ({0}) and y ({1}) are not equal'.format(x, y))
        self.assertTrue(eq)

