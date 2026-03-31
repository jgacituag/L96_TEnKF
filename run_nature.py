#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_nature.py
=============
Generate the nature run and synthetic observations required by run_example.py.

Runs the Lorenz-96 model forward and saves observations with the settings
used for the paper reference case (sigma_obs=5, density=1.0, Type-3 obs).

To generate nature runs for other parameter combinations, modify the
lists at the bottom of this script.

Usage
-----
    python run_nature.py
"""


import os
import sys
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, 'fortran', 'model'))
sys.path.insert(0, os.path.join(ROOT, 'fortran', 'da'))

import codes_experiments.nature_module as nature
import codes_experiments.default_nature_conf as conf

# ─────────────────────────────────────────────────────────
#  Output path
# ─────────────────────────────────────────────────────────
nature_dir = os.path.join(ROOT, 'data', 'Nature')
os.makedirs(nature_dir, exist_ok=True)
conf.GeneralConf['DataPath'] = nature_dir

# Fix random seed for reproducibility
conf.GeneralConf['RandomSeed'] = 10

# ─────────────────────────────────────────────────────────
#  Parameter lists to generate
#  (extend these lists to reproduce additional conditions)
# ─────────────────────────────────────────────────────────
FreqList         = [4]
SpaceDensityList = [1.0]            # 1.0 = all 40 grid points observed
ObsOperatorList  = [3]              # Type 3 = nonlinear (radar-like) obs
ObsErrorList     = [5]              # sigma_obs values; add 0.3, 1, 25 as needed

for MyFreq in FreqList:
    for MyDen in SpaceDensityList:
        for MyOO in ObsOperatorList:
            for MyOE in ObsErrorList:
                conf.ObsConf['Freq']         = MyFreq
                conf.ObsConf['SpaceDensity'] = MyDen
                conf.ObsConf['Type']         = MyOO
                conf.ObsConf['Error']        = MyOE

                exp_name = (f'Freq{MyFreq}_Den{MyDen}_Type{MyOO}_ObsErr{MyOE}')
                conf.GeneralConf['ExpName']        = exp_name
                conf.GeneralConf['NatureFileName'] = f'Paper_Nature_{exp_name}.npz'

                print(f'\nGenerating nature run: {exp_name}')
                nature.nature_run(conf)

print('\nNature runs complete.')
