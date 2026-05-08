#! /usr/bin/env python

import os
import sys
import io
import shutil
import gzip
import errno
import json
import yaml


def file_names_are_unique(file_paths):
    file_names = [os.path.basename(p) for p in file_paths]
    return len(file_names) == len(set(file_names))

def load_yaml(path):
    with open(path, 'r') as stream:
        data = yaml.safe_load(stream)
    return data

def write_yaml(data, path):
    with open(path, 'w') as stream:
        yaml.safe_dump(data, stream)

def load_json(path):
    with open(path, 'r') as in_stream:
        return json.load(in_stream)

def write_json(data, path, indent = 4):
    with open(path, "w") as out_stream:
        json.dump(data, out_stream, indent = indent)

def compress_file(path, compressed_path):
    with open(path, 'rb') as in_stream:
        with gzip.open(compressed_path, 'wb') as out_stream:
            shutil.copyfileobj(in_stream, out_stream)

def decompress_file(compressed_path, path):
    with gzip.open(compressed_path, 'rb') as in_stream:
        with open(path, 'wb') as out_stream:
            shutil.copyfileobj(in_stream, out_stream)

def compress_output_path(path, output_dir = None):
    if not output_dir:
        output_dir = os.path.dirname(path)
    gz_path = os.path.join(
        output_dir,
        f"{os.path.basename(path)}.gz",
    )
    compress_file(path, gz_path)
    return gz_path

def make_directory(path):
    """
    Creates directory `path`, but suppresses error if `path` already exists.
    """
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno == errno.EEXIST:
            pass
        else:
            raise e

def is_gzipped(file_path):
    """
    Returns True if argument is a path to a gzipped file; False otherwise.

    Returns False if path does not exist:
    >>> is_gzipped('')
    False

    Returns False if regular file is empty:
    >>> import tempfile
    >>> fd, temp_path = tempfile.mkstemp()
    >>> os.close(fd)
    >>> is_gzipped(temp_path)
    False

    Returns False for regular file with content:
    >>> f = io.open(temp_path, mode = 'w', encoding='utf-8')
    >>> f.write(u'testing...') #doctest: +ELLIPSIS
    10...
    >>> f.close()
    >>> is_gzipped(temp_path)
    False
    >>> os.remove(temp_path)

    Returns False if gzipped file is empty:
    >>> import tempfile
    >>> fd, temp_path = tempfile.mkstemp()
    >>> f = gzip.open(temp_path, mode = "wt", compresslevel = 9)
    >>> f.write("")
    0
    >>> f.close()
    >>> is_gzipped(temp_path)
    False

    Returns True if file has gzipped content:
    >>> f = gzip.open(temp_path, mode = "wt", compresslevel = 9)
    >>> f.write("testing...")
    10
    >>> f.close()
    >>> is_gzipped(temp_path)
    True
    >>> os.remove(temp_path)
    """

    try:
        fs = gzip.open(file_path)
        d = next(fs)
        fs.close()
    except:
        return False
    return True


class ReadFile(object):
    """
    Obtain a text stream in read mode from a regular or gzipped file.

    Behaves like ``open`` for regular files:
    >>> test_path = os.path.join(os.path.dirname(__file__), os.path.pardir,
    ...         'test-data', 'config.yml')
    >>> if os.path.exists(test_path):
    ...     with ReadFile(test_path) as f:
    ...         l = next(f).strip()
    ...     l == "---"
    ... else:
    ...     True
    ...
    ...
    True
    
    Behaves like ``open`` for gzipped files:
    >>> test_path = os.path.join(os.path.dirname(__file__), os.path.pardir,
    ...         'test-data', 'trees', 'crocs-1.trees.gz')
    >>> if os.path.exists(test_path):
    ...     with ReadFile(test_path) as f:
    ...         l = next(f).strip()
    ...     l == "#NEXUS"
    ... else:
    ...     True
    ...
    ...
    True
    """

    open_files = set()

    def __init__(self, path):
        self.path = path
        self.gzipped = is_gzipped(self.path)
        self.encoding = 'utf-8'
        if self.gzipped:
            try:
                self.file_stream = gzip.open(filename = self.path, mode = 'rt',
                        encoding = self.encoding)
            except (TypeError, ValueError):
                self.file_stream = gzip.open(filename = self.path, mode = 'r')
        else:
            self.file_stream = io.open(self.path, mode = 'r',
                    encoding = self.encoding)
        self.__class__.open_files.add(self.path)

    def close(self):
        self.file_stream.close()
        self.__class__.open_files.remove(self.path)

    def __enter__(self):
        return self.file_stream

    def __exit__(self, type, value, traceback):
        self.close()

    def __iter__(self):
        return self.file_stream.__iter__()

    def next(self):
        return next(self.file_stream)

    def read(self):
        return self.file_stream.read()

    def readline(self):
        return self.file_stream.readline()

    def readlines(self):
        return self.file_stream.readlines()

    def _get_closed(self):
        return self.file_stream.closed

    closed = property(_get_closed)

    def seek(self, i):
        self.file_stream.seek(i)

    def flush(self):
        self.file_stream.flush()

    def fileno(self):
        return self.file_stream.fileno()
