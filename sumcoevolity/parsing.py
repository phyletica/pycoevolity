#! /usr/bin/env python

import sys
import os
import logging
import re

from sumcoevolity.fileio import ReadFile 

_LOG = logging.getLogger(__name__)

HEADER_PATTERN = re.compile(r'^\s*\D.+')
NAN_PATTERN = re.compile(r'^\s*nan\s*$', re.IGNORECASE)

def line_count(stream):
    count = 0
    for line in stream:
        count += 1
    return count

def parse_header(file_stream, sep = '\t', strict = True, seek = True):
    try:
        header_line = next(file_stream)
    except StopIteration:
        if strict:
            raise Exception('did not find header in {0}'.format(file_stream.name))
        else:
            return None
    if not HEADER_PATTERN.match(header_line):
        if strict:
            raise Exception('did not find header in {0}'.format(file_stream.name))
        else:
            return None
    header = header_line.strip().split(sep)
    if seek:
        file_stream.seek(0)
    return header

def parse_header_from_path(path, sep = '\t', strict = True, seek = True):
    with ReadFile(path) as stream:
        return parse_header(stream, sep = sep, strict = strict, seek = seek)

def spreadsheet_iter(spreadsheets, sep = '\t', header = None):
    head_line = False
    if not header:
        head_line = True
        header = parse_header_from_path(spreadsheets[0], sep = sep)
    for sheet_idx, ss in enumerate(spreadsheets):
        with ReadFile(ss) as file_stream:
            if head_line:
                h = next(file_stream).strip().split(sep)
                if header != h:
                    raise Exception('headers do not match')
            for row_idx, row in enumerate(file_stream):
                if row.strip() == '':
                    continue
                r = [el.strip() for el in row.strip().split(sep)]
                if len(r) != len(header):
                    raise Exception('row {0} of spreadsheet {1} has {2} columns, '
                            'header has {3}'.format(row_idx + 1, sheet_idx + 1,
                                    len(r), len(header)))
                yield dict(zip(header, r))

def dict_line_iter(d, sep = '\t', header = None):
    if not header:
        header = sorted(d.keys())
    if sorted(header) != sorted(d.keys()):
        raise ValueError('header does not match dict keys')
    yield '{0}\n'.format(sep.join(header))
    for i in range(len(d[header[0]])):
        yield '{0}\n'.format(sep.join([str(d[h][i]) for h in header]))

def get_dict_from_spreadsheets(spreadsheets, sep = '\t', header = None):
    ss_iter = spreadsheet_iter(spreadsheets, sep = sep, header = header)
    row_dict = next(ss_iter)
    d = dict(zip(row_dict.keys(),
            [[row_dict[k]] for k in row_dict.keys()]))
    for row_dict in ss_iter:
        for k in row_dict.keys():
            d[k].append(row_dict[k])
    return d

