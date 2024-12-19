# Pre-eruptive-water-content-calculation

This is the introduction of a discretization method to constrain the pre-eruptive water content of magma.

Di et al. ([2020](https://doi.org/10.2138/am-2020-7137)) proposed a new calculation method to constrain the pre-eruptive water content of magma: by combining a water sensitive thermometer (Lee et al., [2009](https://doi.org/10.1016/j.epsl.2008.12.020)) and a water insensitive thermometer (Putirka et al., [2003](https://doi.org/10.2138/am-2003-1017)) to determine the pre eruption water content of the melt. 
Thermometers sensitive to water are based on the properties of silicon activity in the equilibrium system of olivine clinopyroxene melt, while thermometers insensitive to water are based on the equilibrium process of clinopyroxene melt. 

In theory, for the same volcano, the P-T trajectories given by two sets of thermometers should be similar. Therefore, we can adjust the water content in the water sensitive silicon activity thermometers to achieve maximum overlap between the P-T trajectories calculated by the two sets of thermometers. At this point, the water content is the pre eruption water content of the basalt that we fit.


The input `data` should contain three files, including:

1. Clinopyroxene composition (file: `data/cpx_composition`)

**#SiO<sub>2</sub>	 TiO<sub>2</sub>	Al<sub>2</sub>O<sub>3</sub> FeO<sub>t</sub>	MnO	MgO	CaO Na<sub>2</sub>O K<sub>2</sub>O Cr<sub>2</sub>O<sub>3</sub>**

51.81 	0.14 	1.69 	7.82 	0.39 	14.25 	24.08 	0.30 	0.01 	0.00 ...

2. Melt composition equilibrated with cpx (file: `data/liquid_composition`). **Note**:  1 and 2 should be matched one by one.

**#SiO<sub>2</sub>	 TiO<sub>2</sub>	Al<sub>2</sub>O<sub>3</sub> FeO<sub>t</sub>	MnO	MgO	CaO	Na<sub>2</sub>O	K<sub>2</sub>O	Cr<sub>2</sub>O<sub>3</sub> P<sub>2</sub>O<sub>5</sub>**

44.08	2.43	13.19	17.04	0.27	7.61	11.05	2.16	0.18	0	0.25 ...

3. Primary melt composition (file: `data/primary_melt_composition`):

**# SiO<sub>2</sub> TiO<sub>2</sub> Al<sub>2</sub>O<sub>3</sub> TFe<sub>2</sub>O<sub>3</sub> MgO CaO Na<sub>2</sub>O K<sub>2</sub>O P<sub>2</sub>O<sub>5</sub> TOTAL**

52.74 	2.45 	13.88 	8.49 	0.11 	6.17 	5.50 	3.86 	5.53 	99.69 
...

(The order of the species does not matter, while all the species must be included. The input primary melt compositions should be in equlibrium with olivine + orthopyroxene)


**Manual**:
1. Copy and Paste your primary melt composition and cpx-melt data into the three files mentioned above respectively.
2. Change the parameters in `config.py`. The changeable parameters include a melt Fe</sub>3+ proportion (Fe<sup>3+</sup>/FeT), errors of the thermobarometers, intended H<sub>2</sub>O content range and increment, and the size of unit.
3. Run `demo.py` to get the outcome.
4. The outcome is a figure showing a scatter of H<sub>2</sub>O content (.wt%) and overlap area as well as printing them on the screen.
