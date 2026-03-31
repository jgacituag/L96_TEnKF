#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
default_nature_conf.py
======================
Default configuration for Lorenz-96 nature runs and synthetic observations.
All paths use relative references; run_nature.py overrides DataPath at runtime.
"""
import numpy as np

# =================================================================
# GENERAL SECTION
# =================================================================
GeneralConf = dict()

GeneralConf['ExpName']        = 'DefaultNature'
GeneralConf['DataPath']       = './data/Nature/'
GeneralConf['NatureFileName'] = 'Nature' + GeneralConf['ExpName'] + '.npz'
GeneralConf['RandomSeed']     = 10   # fixed for reproducibility

# =================================================================
# MODEL SECTION
# =================================================================
ModelConf = dict()

ModelConf['nx']  = 40      # Number of large-scale state variables
ModelConf['dt']  = 0.0125  # Time step (do not change)

ModelConf['Coef']  = np.array([8])
ModelConf['NCoef'] = np.size(ModelConf['Coef'])

ModelConf['FSpaceDependent'] = False
ModelConf['FSpaceAmplitude'] = np.array([1])
ModelConf['FSpaceFreq']      = np.array([1])

ModelConf['EnablePRF'] = False
ModelConf['CSigma']    = np.array([0])
ModelConf['CPhi']      = 1.0

ModelConf['EnableSRF'] = False
ModelConf['XSigma']    = 0.0
ModelConf['XPhi']      = 1.0
ModelConf['XLoc']      = np.arange(1, ModelConf['nx'] + 1)

ModelConf['TwoScaleParameters'] = np.array([10, 10, 0])  # Hint=0: single-scale
ModelConf['nxss'] = ModelConf['nx'] * 8
ModelConf['dtss'] = ModelConf['dt'] / 5

# =================================================================
# NATURE RUN SECTION
# =================================================================
NatureConf = dict()

NatureConf['NEns']    = 1      # single truth trajectory
NatureConf['RunSave'] = True
NatureConf['RunPlot'] = False
NatureConf['SPLength'] = 40    # spin-up (model time units)
NatureConf['Length']   = 4000  # run length (model time units)

# =================================================================
# OBSERVATION CONFIGURATION SECTION
# =================================================================
ObsConf = dict()

ObsConf['Freq']         = 4          # observation frequency (time steps)
ObsConf['NetworkType']  = 'regular'
ObsConf['SpaceDensity'] = 1.0        # 1.0 = all 40 grid points
ObsConf['TimeDensity']  = 1
ObsConf['Error']        = 1.0        # obs error std — overridden per experiment
ObsConf['Bias']         = 0.0
ObsConf['Type']         = 1          # 1=linear, 3=nonlinear — overridden per exp
