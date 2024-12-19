# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 10:18:31 2021

@author: DXY
"""

import numpy as np
import pandas as pd
import config

def cal_Si_T_p(percent_H2O, primary):
    # calculate mole composition
    mole_compos = pd.DataFrame(
        columns=['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O'])
    mole_compos['SiO2'] = [i / config.M["SiO2"] for i in primary['SiO2']]
    mole_compos['TiO2'] = [i / config.M["TiO2"] for i in primary['TiO2']]
    mole_compos['Al2O3'] = [i / config.M["Al2O3"] for i in primary['Al2O3']]
    mole_compos['Fe2O3'] = [i / config.M["Fe2O3"] * 2 * config.proportion_Fe3 / 2 for i in primary['TFe2O3']]
    mole_compos['FeO'] = [i / config.M["Fe2O3"] * 2 * (1 - config.proportion_Fe3) for i in primary['TFe2O3']]
    mole_compos['MnO'] = [i / config.M["MnO"] for i in primary['MnO']]
    mole_compos['MgO'] = [i / config.M["MgO"] for i in primary['MgO']]
    mole_compos['CaO'] = [i / config.M["CaO"] for i in primary['CaO']]
    mole_compos['Na2O'] = [i / config.M["Na2O"] for i in primary['Na2O']]
    mole_compos['K2O'] = [i / config.M["K2O"] for i in primary['K2O']]
    mole_compos['P2O5'] = [i / config.M["P2O5"] for i in primary['P2O5']]
    mole_compos['H2O'] = [percent_H2O * 0.01 * i / (1 - percent_H2O * 0.01) / config.M["H2O"] for i in primary['TOTAL']]

    # calculate species composition
    species = pd.DataFrame(
        columns=['Si4O8', 'Ti4O8', 'Al16/3O8', 'Fe16/3O8', 'Fe4Si2O8', 'Mg4Si2O8', 'Ca4Si2O8', 'Na2Al2Si2O8',
                 'K2Al2Si2O8', 'P16/5O8', 'H16O8'],
        index=mole_compos.index)

    for i in range(len(mole_compos['SiO2'])):
        species['Si4O8'][i] = 0.25 * (mole_compos['SiO2'][i] - 0.5 * (
                mole_compos['FeO'][i] + mole_compos['MgO'][i] + mole_compos['CaO'][i]) - mole_compos['Na2O'][i] -
                                      mole_compos['K2O'][i])

    for i in range(len(mole_compos['Al2O3'])):
        species['Al16/3O8'][i] = 0.375 * (mole_compos['Al2O3'][i] - mole_compos['Na2O'][i])

    species['Ti4O8'] = [0.25 * i for i in mole_compos['TiO2']]
    species['Fe16/3O8'] = [0.375 * i for i in mole_compos['Fe2O3']]
    species['Fe4Si2O8'] = [0.25 * i for i in mole_compos['FeO']]
    species['Mg4Si2O8'] = [0.25 * i for i in mole_compos['MgO']]
    species['Ca4Si2O8'] = [0.25 * i for i in mole_compos['CaO']]
    species['Na2Al2Si2O8'] = [i for i in mole_compos['Na2O']]
    species['K2Al2Si2O8'] = [i for i in mole_compos['K2O']]
    species['P16/5O8'] = [0.625 * i for i in mole_compos['P2O5']]
    species['H16O8'] = [0.125 * i for i in mole_compos['H2O']]
    species_percent = species.div(species.sum(axis=1), axis=0) * 100

    # calculate t and p
    temperature = []  # K
    pressure = []  # kbar

    for i in range(len(species_percent)):  # Equation 3 in Lee et al. (2009)
        temperature.append(916.45 + 13.68 * species_percent['Mg4Si2O8'][i] + 4580 /
                           species_percent['Si4O8'][i] - 0.509 * species_percent['H16O8'][i] *
                           species_percent['Mg4Si2O8'][i] + 273.15)

    for i in range(len(species_percent)):  # Equation 2 in Lee et al. (2009)
        pressure.append(10 * (
                np.log(species_percent['Si4O8'][i]) - 4.019 + 0.0165 * species_percent['Fe4Si2O8'][i] + 0.0005 *
                species_percent['Ca4Si2O8'][i] ** 2) /
                        (-770 * (temperature[i]) ** (-1) + 0.0058 * (temperature[i]) ** 0.5 - 0.003 *
                         species_percent['H16O8'][i]))

    return temperature, pressure

def cal_cpx_T_p(cpx_compos, liq_compos):
# calculate liquid mole composition
    liq = pd.DataFrame(
        columns=['SiO2', 'TiO2', 'AlO3/2', 'FeO', 'MnO', 'MgO', 'CaO', 'NaO1/2', 'KO1/2', 'CrO3/2', 'PO5/2'])
    liq['SiO2'] = [i / config.M["SiO2"] for i in liq_compos['SiO2']]
    liq['TiO2'] = [i / config.M["TiO2"] for i in liq_compos['TiO2']]
    liq['AlO3/2'] = [i * 2 / config.M["Al2O3"] for i in liq_compos['Al2O3']]
    liq['FeO'] = [i / config.M["FeO"] for i in liq_compos['FeOT']]
    liq['MnO'] = [i / config.M["MnO"] for i in liq_compos['MnO']]
    liq['MgO'] = [i / config.M["MgO"] for i in liq_compos['MgO']]
    liq['CaO'] = [i / config.M["CaO"] for i in liq_compos['CaO']]
    liq['NaO1/2'] = [i * 2 / config.M["Na2O"] for i in liq_compos['Na2O']]
    liq['KO1/2'] = [i * 2 / config.M["K2O"] for i in liq_compos['K2O']]
    liq['CrO3/2'] = [i * 2 / config.M["Cr2O3"] for i in liq_compos['Cr2O3']]
    liq['PO5/2'] = [i * 2 / config.M["P2O5"] for i in liq_compos['P2O5']]
    liq_percent = liq.div(liq.sum(axis=1), axis=0)
    liq_percent['Mg#'] = liq_percent['MgO'] / (liq_percent['MgO'] + liq_percent['FeO'])
    
    # calculate cpx composition
    cpx = pd.DataFrame(
        columns=['SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'Cr2O3'])
    cpx['SiO2'] = [i / config.M["SiO2"] for i in cpx_compos['SiO2']]
    cpx['TiO2'] = [i / config.M["TiO2"] for i in cpx_compos['TiO2']]
    cpx['Al2O3'] = [i / config.M["Al2O3"] for i in cpx_compos['Al2O3']]
    cpx['FeO'] = [i / config.M["FeO"] for i in cpx_compos['FeOT']]
    cpx['MnO'] = [i / config.M["MnO"] for i in cpx_compos['MnO']]
    cpx['MgO'] = [i / config.M["MgO"] for i in cpx_compos['MgO']]
    cpx['CaO'] = [i / config.M["CaO"] for i in cpx_compos['CaO']]
    cpx['Na2O'] = [i / config.M["Na2O"] for i in cpx_compos['Na2O']]
    cpx['K2O'] = [i / config.M["K2O"] for i in cpx_compos['K2O']]
    cpx['Cr2O3'] = [i / config.M["Cr2O3"] for i in cpx_compos['Cr2O3']]
    
    Oxy_num = pd.DataFrame()
    Oxy_num['SiO2'] = [i * 2 for i in cpx['SiO2']]
    Oxy_num['TiO2'] = [i * 2 for i in cpx['TiO2']]
    Oxy_num['Al2O3'] = [i * 3 for i in cpx['Al2O3']]
    Oxy_num['FeO'] = [i for i in cpx['FeO']]
    Oxy_num['MnO'] = [i for i in cpx['MnO']]
    Oxy_num['MgO'] = [i for i in cpx['MgO']]
    Oxy_num['CaO'] = [i for i in cpx['CaO']]
    Oxy_num['Na2O'] = [i for i in cpx['Na2O']]
    Oxy_num['K2O'] = [i for i in cpx['K2O']]
    Oxy_num['Cr2O3'] = [i * 3 for i in cpx['Cr2O3']]
    Oxy_num['sum'] = Oxy_num.sum(axis=1)
    Oxy_num['renorm factor'] = [6 / i for i in Oxy_num['sum']]
    
    cpx_new = pd.DataFrame(
        columns=['Si', 'Ti', 'Al_IV', 'Al_VI', 'Al_T', 'Fe', 'Mn', 'Mg', 'Ca', 'Na', 'K', 'Cr'])
    cpx_new['Si'] = cpx['SiO2'] * Oxy_num['renorm factor']
    cpx_new['Ti'] = cpx['TiO2'] * Oxy_num['renorm factor']
    cpx_new['Al_T'] = cpx['Al2O3'] * Oxy_num['renorm factor'] * 2
    cpx_new['Al_IV'] = 2 - cpx_new['Si']
    for i in range(len(cpx_new)):
        cpx_new['Al_VI'][i] = cpx_new['Al_T'][i] - cpx_new['Al_IV'][i] if cpx_new['Al_T'][i] > cpx_new['Al_IV'][i] else 0
    cpx_new['Fe'] = cpx['FeO'] * Oxy_num['renorm factor']
    cpx_new['Mn'] = cpx['MnO'] * Oxy_num['renorm factor']
    cpx_new['Mg'] = cpx['MgO'] * Oxy_num['renorm factor']
    cpx_new['Ca'] = cpx['CaO'] * Oxy_num['renorm factor']
    cpx_new['Na'] = cpx['Na2O'] * Oxy_num['renorm factor'] * 2
    cpx_new['K'] = cpx['K2O'] * Oxy_num['renorm factor'] * 2
    cpx_new['Cr'] = cpx['Cr2O3'] * Oxy_num['renorm factor'] * 2
    
    # calculate cpx species
    cpx_species = pd.DataFrame(
        columns=['Jd', 'CaTs', 'CaTi', 'CrCaTs', 'DiHd', 'EnFs', 'FmCaTs', 'FmTi'],
        index=cpx_new.index)
    
    for i in range(len(cpx_new)):
        cpx_species['Jd'][i] = cpx_new['Na'][i] if cpx_new['Al_VI'][i] > cpx_new['Na'][i] else cpx_new['Al_VI'][i]
        cpx_species['CaTs'][i] = cpx_new['Al_VI'][i] - cpx_new['Na'][i] if cpx_new['Al_VI'][i] > cpx_new['Na'][i] else 0
        cpx_species['CaTi'][i] = (cpx_new['Al_IV'][i] - cpx_species['CaTs'][i]) / 2 if cpx_new['Al_IV'][i] > cpx_species['CaTs'][i] else 0
        cpx_species['CrCaTs'][i] = cpx_new['Cr'][i] / 2
        cpx_species['DiHd'][i] = cpx_new['Ca'][i] - cpx_species['CaTs'][i] - cpx_species['CaTi'][i] - cpx_species['CrCaTs'][i] if cpx_new['Ca'][i] - cpx_species['CaTs'][i] - cpx_species['CaTi'][i] - cpx_species['CrCaTs'][i] > 0 else 0
        cpx_species['EnFs'][i] = (cpx_new['Mg'][i] + cpx_new['Fe'][i] - cpx_species['DiHd'][i]) / 2
        # for Ca-poor pyroxene
        if cpx_new['Ca'][i] < cpx_species['CaTs'][i]:
            cpx_species['CaTs'][i] = cpx_new['Ca'][i]
            cpx_species['CaTi'][i] = cpx_species['CrCaTs'][i] = cpx_species['DiHd'][i] = 0
            cpx_species['FmCaTs'][i] = cpx_new['Al_VI'][i] - cpx_species['CaTs'][i] * 2  ## why 2 times?
            cpx_species['FmTi'][i] = (cpx_new['Al_IV'][i] - cpx_species['CaTs'][i] - cpx_species['FmCaTs'][i]) / 2
            cpx_species['EnFs'] = (cpx_new['Mg'][i] + cpx_new['Fe'][i] - cpx_species['FmCaTs'][i] - cpx_species['FmTi'][i]) / 2
    
    # two equilibrium constant
    K_Jd_liq = cpx_species['Jd'] / (liq_percent['NaO1/2'] * liq_percent['AlO3/2'] * liq_percent['SiO2'] ** 2)
    
    K_Jd_DiHd = (cpx_species['Jd'] * liq_percent['CaO'] * (liq_percent['FeO'] + liq_percent['MgO'])) / (
                cpx_species['DiHd'] * liq_percent['NaO1/2'] * liq_percent['AlO3/2'])
    
    # calculate T and p
    # initial T and p
    T_ini = []
    p_ini = []
    for i in range(len(liq_percent)):
        T_ini.append(10 ** 4 / (6.73 - 0.26 * np.log(K_Jd_DiHd[i]) - 0.86 * np.log(liq_percent['Mg#'][i]) + 0.52 * np.log(liq_percent['CaO'][i])))
    for i in range(len(T_ini)):
        p_ini.append(- 88.3 + 0.00282 * T_ini[i] * np.log(K_Jd_liq[i]) + 0.0219 * T_ini[i] - 25.1 * np.log(liq_percent['CaO'][i] * liq_percent['SiO2'][i]) + 7.03 * liq_percent['Mg#'][i] + 12.4 * np.log(liq_percent['CaO'][i]))
    
    # iterative calculations
    def iter_cal(T0, p0, K_Jd_DiHd, liq_MgN, liq_Na, liq_Si, cpx_Jd, K_Jd_liq, liq_Ca):
        def cal_T(p):
            return 10 ** 4 / (4.60 - 0.437 * np.log(K_Jd_DiHd) - 0.654 * np.log(liq_MgN) - 0.326 * np.log(liq_Na) - 0.00632 * p - 0.92 * np.log(liq_Si) + 0.274 * np.log(cpx_Jd))
    
        def cal_p(T):
            return (-88.3 + 0.00282 * T * np.log(K_Jd_liq) + 0.0219 * T - 25.1 * np.log(liq_Ca * liq_Si) + 7.03 * liq_MgN + 12.4 * np.log(liq_Ca))
        
        for i in range(100):
            p=cal_p(T0)
            T=cal_T(p0)
            if abs(p-p0)<1e-3 and abs(T-T0)<1e-3:
                return p, T
            p0=p
            T0=T
    
    temperature=[]
    pressure=[]
    for i in range(len(cpx_compos)):
        p, T = iter_cal(T_ini[i], p_ini[i], K_Jd_DiHd[i], liq_percent['Mg#'][i], liq_percent['NaO1/2'][i], liq_percent['SiO2'][i], cpx_species['Jd'][i], K_Jd_liq[i], liq_percent['CaO'][i])
        pressure.append(p)
        temperature.append(T)
    
    return temperature, pressure

def cal_overlap_area(x1_array, y1_array, x2_array, y2_array,
					 xerr1, yerr1, xerr2, yerr2,
					 left, right, down, up):
	area_total = 0
	unit_coefficient = 10 if config.H2O_end - config.H2O_begin > 5 else 4
	d_x = config.d_x * unit_coefficient
	d_y = config.d_y * unit_coefficient
	d_area = d_x * d_y
	for i in np.arange(left, right, d_x):
		for j in np.arange(down, up, d_y):
			flag = 0
			for m in range(len(x1_array)):
				for n in range(len(x2_array)):
					if (abs(i - x1_array[m]) < xerr1
						and j < y1_array[m] + yerr1 * np.sqrt(1 - (i - x1_array[m]) ** 2 / xerr1 ** 2)
						and j > y1_array[m] - yerr1 * np.sqrt(1 - (i - x1_array[m]) ** 2 / xerr1 ** 2)
						and abs(i - x2_array[n]) < xerr2
						and j < y2_array[n] + yerr2 * np.sqrt(1 - (i - x2_array[n]) ** 2 / xerr2 ** 2)
						and j > y2_array[n] - yerr2 * np.sqrt(1 - (i - x2_array[n]) ** 2 / xerr2 ** 2)):
						area_total += d_area
						flag = 1 # to make sure each unit will only be counted once
						break
				if flag == 1:
					break
	return area_total

















