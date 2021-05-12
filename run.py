#! /usr/bin/env python
# -*- coding: utf-8 -*-
from PrirodaOptimizer import Conformer, PrirodaOptimizer


c = Conformer.from_xyz('data/methane.xyz')

p = PrirodaOptimizer(tmp_ram=-1000, n_jobs=3)
o = p.transform([c]*10)
