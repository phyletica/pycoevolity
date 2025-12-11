#! /usr/bin/env python

import os

LOCAL_DIR = os.path.abspath(os.path.dirname(__file__))
TEST_DIR = os.path.dirname(LOCAL_DIR)
TEST_DATA_DIR = os.path.join(TEST_DIR, "data")
TEST_OUTPUT_DIR = os.path.join(TEST_DIR, "output")

def data_path(filename=""):
    return os.path.join(TEST_DATA_DIR, filename)

def output_path(filename=""):
    return os.path.join(TEST_OUTPUT_DIR, filename)

def test_path(filename=""):
    return os.path.join(TEST_DIR, filename)
