# Disclaimer
**This repository contains code and data related to ongoing, non-publishable research. The materials here are part of a work in progress and are not intended for publication or external use.** <br>
For any inquiries or further information, please contact Wook Shin via wook.shin@vanderbilt.edu <br>
The content within includes the analysis of protein entropy based on trajectories from molecular dynamics simulations conducted on Anton supercomputers by D. E. Shaw Research. The use of this [simulation data](https://www.deshawresearch.com/downloads/download_trajectory_sarscov2.cgi/) in **any published work** should be acknowledged by including a citation to: **D. E. Shaw Research, "Molecular Dynamics Simulations Related to SARS-CoV-2," D. E. Shaw Research Technical Data, 2020.**

# DESRES-ANTON-[15235444, 15235449, 15235455, 15256595, 15256598, 15256602]
The dataset released by D. E. Shaw Research contains 6 500-μs simulations of SARS-CoV-2 NiRAN domain targeted peptides. D. E. Shaw Research utilized an in-house peptide-optimization workflow, which combines machine learning with free energy calculations, to predict β-hairpin peptides that would have high affinity for the N-terminal binding site on the nidovirus RdRp-associated nucleotidyltransferase (NiRAN) domain. Their workflow generated a series of amphiphilic peptides that contain both a hydrophobic surface, which in their computational model bound to the NiRAN domain, and a highly positively charged, arginine- and lysine-rich surface that faced away from the NiRAN domain.
![Picture1](https://github.com/rshin1209/desres/assets/25111091/ee9d10e4-506a-4598-af60-1792365ea567)
**Table 1. Amino acid sequences of the peptides are presented.** In the sequence alignment, black font indicates that the residue is present in the wild-type peptide and red indicates that the residue is not present in the wild-type peptide. "Peptide-X-dis" labels denote crosslinked peptides.

# Protein Entropy Analysis
## Python Module Requirements:
* MDtraj
* Networkx
* Numpy
* json
* multiprocessing
## Reproduction
Protein Entropy Analysis was performed on six (peptide-only) NiRAN domain-targeted peptides. (1) A.pdb. (2) B.pdb. (3) C.pdb. (4) A_dis.pdb. (5) B_dis.pdb. (6) C_dis.pdb.<br \>
The main Python script `protein_entropy.py` is used for implementing protein entropy analysis.<br>
Taking the (1) for an example, one can simply run the following two commands to perform protein entropy analysis.<br>
`python protein_entropy.py --reaction A --temperature 298.15 --job 0`
'python protein_entropy.py --reaction A --temperature 298.15 --job 1'
`[reaction] -- Name of pdb file containing simulation trajectory frames`
`[temperature] -- Temperature (Kelvin) for -TS calculation`
`[job] -- 0 for degrees of freedom extraction from Cartesian coordinates and 1 for configurational entropy calculation`


Configurational entropy calculation of a protein begins with extracting Internal Coordinates from Cartesian Coordinates of every frame that exists in simulation trajectories. 

![A_dis_imshow](https://github.com/rshin1209/desres/assets/25111091/8291d79a-d788-486d-85a0-5218dcb194d3)
