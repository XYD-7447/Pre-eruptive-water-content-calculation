# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 10:32:17 2022

@author: DXY
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from os import path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__) ) ) )
import config
import utils_new
import pandas as pd

# if __name__ == "__main__":
Si_primary = pd.read_csv('./data/primary_melt_composition.txt', delimiter='\s+', skiprows=1,
                      names=['SiO2', 'TiO2', 'Al2O3', 'TFe2O3', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'TOTAL'])

cpx_compos = pd.read_csv('./data/cpx_composition.txt', delimiter='\s+', skiprows=1,
                         names=['SiO2', 'TiO2', 'Al2O3', 'FeOT', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'Cr2O3'])

liq_compos = pd.read_csv('./data/liquid_composition.txt', delimiter='\s+', skiprows=1,
                         names=['SiO2', 'TiO2', 'Al2O3', 'FeOT', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'Cr2O3', 'P2O5'])

#Si_T_p=utils_new.cal_Si_T_p(1.3,Si_primary)
trend = []
for percent_H2O in np.arange(config.H2O_begin, config.H2O_end, config.H2O_interval):
    Si_T, Si_p = utils_new.cal_Si_T_p(percent_H2O, Si_primary)
    cpx_T, cpx_p = utils_new.cal_cpx_T_p(cpx_compos, liq_compos)
    left = min(Si_T) - config.E_B_T
    right = max(cpx_T) + config.E_A_T
    down = min(Si_p) - config.E_B_p
    up = max(cpx_p) + config.E_A_p
    area = utils_new.cal_overlap_area(cpx_T, cpx_p, Si_T, Si_p,
 									  config.E_A_T, config.E_A_p, config.E_B_T, config.E_B_p,
 									  left, right, down, up)
    trend.append([percent_H2O, area])
    print('%s, %s' %(round(percent_H2O, 2), round(area, 2)))

plt.scatter(np.array(trend)[:, 0], np.array(trend)[:, 1], color = 'k')
plt.xlabel('H2O wt.%')
plt.ylabel('overlap degree')
plt.show()