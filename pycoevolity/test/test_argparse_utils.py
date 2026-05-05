#! /usr/bin/env python

import os
import pytest
from argparse import ArgumentTypeError

from pycoevolity import argparse_utils


class TestIsPath:
    def test_valid_file(self):
        p = os.path.abspath(__file__)
        returned = argparse_utils.arg_is_path(p)
        assert returned == p

    def test_valid_dir(self):
        p = os.path.dirname(os.path.abspath(__file__))
        returned = argparse_utils.arg_is_path(p)
        assert returned == p

    def test_bad_path(self):
        p = "/not/a/likely/path/on/your/computer/4dk8093k1403"
        with pytest.raises(ArgumentTypeError):
            argparse_utils.arg_is_path(p)

class TestIsFile:
    def test_valid_file(self):
        p = os.path.abspath(__file__)
        returned = argparse_utils.arg_is_file(p)
        assert returned == p

    def test_valid_dir(self):
        p = os.path.dirname(os.path.abspath(__file__))
        with pytest.raises(ArgumentTypeError):
            argparse_utils.arg_is_file(p)

    def test_bad_path(self):
        p = "/not/a/likely/path/on/your/computer/4dk8093k1403"
        with pytest.raises(ArgumentTypeError):
            argparse_utils.arg_is_file(p)

class TestIsDir:
    def test_valid_file(self):
        p = os.path.abspath(__file__)
        with pytest.raises(ArgumentTypeError):
            argparse_utils.arg_is_dir(p)

    def test_valid_dir(self):
        p = os.path.dirname(os.path.abspath(__file__))
        returned = argparse_utils.arg_is_dir(p)
        assert returned == p

    def test_bad_path(self):
        p = "/not/a/likely/path/on/your/computer/4dk8093k1403"
        with pytest.raises(ArgumentTypeError):
            argparse_utils.arg_is_dir(p)

class TestIsDirOrNewDir:
    def test_valid_file(self):
        p = os.path.abspath(__file__)
        with pytest.raises(ArgumentTypeError):
            argparse_utils.arg_is_dir_or_new_dir(p)

    def test_valid_dir(self):
        p = os.path.dirname(os.path.abspath(__file__))
        returned = argparse_utils.arg_is_dir_or_new_dir(p)
        assert returned == p

    def test_valid_parent_dir(self):
        p = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                "a-child-unlikely-to-exist-iwouervy549873")
        returned = argparse_utils.arg_is_dir_or_new_dir(p)
        assert returned == p

    def test_bad_path(self):
        p = "/not/a/likely/path/on/your/computer/4dk8093k1403"
        with pytest.raises(ArgumentTypeError):
            argparse_utils.arg_is_dir_or_new_dir(p)

class TestIsPositiveInt:
    def test_zero(self):
        with pytest.raises(ArgumentTypeError):
            argparse_utils.arg_is_positive_int('0')

    def test_valid_ints(self):
        i = argparse_utils.arg_is_positive_int('1')
        assert i == 1
        i = argparse_utils.arg_is_positive_int('3')
        assert i == 3

    def test_neg_int(self):
        with pytest.raises(ArgumentTypeError):
            argparse_utils.arg_is_positive_int('-1')

    def test_float(self):
        with pytest.raises(ArgumentTypeError):
            argparse_utils.arg_is_positive_int('3.14')

    def test_string(self):
        with pytest.raises(ArgumentTypeError):
            argparse_utils.arg_is_positive_int('blah')
