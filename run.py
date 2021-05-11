#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
from PrirodaOptimizer import Conformer, PrirodaOptimizer
#from PrirodaOptimizer.conformer import Conformer
#from PrirodaOptimizer.transformer import PrirodaOptimizer


def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", default='Alanine.xyz')
    parser.add_argument("-t", "--task", default='input')
    return parser


if __name__ == "__main__":
    parser = createParser()
    namespace = parser.parse_args(sys.argv[1:])

    p = PrirodaOptimizer()
    root = './PrirodaOptimizer/data/'

    if namespace.task.lower() == "input":
        c = Conformer.from_xyz(root + namespace.file)
        print(p._prepare_input(c, tmp_dir="/tmp"))

    elif namespace.task.lower() == "output":
        with open(root + namespace.file, "r") as log_file:
            print(p._prepare_input(p._parse_output(log_file), tmp_dir="/tmp"))
