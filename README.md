# Overview
This respository provides additional data on the publication "Phenylalanine Modification in Plasma-Driven Biocatalysis Revealed by Solvent Accessibility and Reactive Dynamics in Combination with Protein Mass Spectrometry" by Poggemann et al.[^1].
All simulations were carried out with the LAMMPS simulation package [^2]. The SASA results featured here and in the publication were calculated with the SASA-Analysis Package (https://github.com/hpoggemann/SASA-Analysis).

# Data Structure:

- 00_General_files includes:
	- the force filed used for all calculations
	- the molecule files for all small probe molecues investigated
- 01_SASA includes:
	- The SASA.py files for each enzyme
	- The starting structures for each enzyme
	- SASA result files for CviUPO, AaeUPO and GapA
	- The corresponding plots can be found in the paper SI
- 02_MDs includes:
	- The gerneral input file for the short MDs
	- The gerneral input file for the long MDs
	- Starting structures for all enzymes with and without solvent
	- 3 excel sheets with the details results of long and short MDs

# Citations
[^1]: H.-F. Poggemann et al., “Phenylalanine modification in plasma-driven biocatalysis revealed by solvent accessibility and reactive dynamics in combination with protein mass spectrometry,” The Journal of Physical Chemistry B, 2025, doi: 10.1021/acs.jpcb.5c03518. (https://pubs.acs.org/doi/10.1021/acs.jpcb.5c03518)
[^2]: A. P. Thompson et al., “LAMMPS - a flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales,” Computer Physics Communications, vol. 271, p. 108171, Feb. 2022, doi: 10.1016/j.cpc.2021.108171. (https://doi.org/10.1016/j.cpc.2021.108171)
