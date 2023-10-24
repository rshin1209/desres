# Disclaimer
**This repository contains code and data related to ongoing, non-publishable research. The materials here are part of a work in progress and are not intended for publication or external use.** <br>

For any inquiries or further information, please contact Wook Shin via wook.shin@vanderbilt.edu <br>

The content within includes the analysis of protein entropy based on trajectories from molecular dynamics simulations conducted on Anton supercomputers by D. E. Shaw Research. The use of this [simulation data](https://www.deshawresearch.com/downloads/download_trajectory_sarscov2.cgi/) in **any published work** should be acknowledged by including a citation to: **D. E. Shaw Research, "Molecular Dynamics Simulations Related to SARS-CoV-2," D. E. Shaw Research Technical Data, 2020.**

# DESRES-ANTON-[15235444, 15235449, 15235455, 15256595, 15256598, 15256602]
The dataset released by D. E. Shaw Research contains six 500-μs simulations of SARS-CoV-2 NiRAN domain targeted peptides. D. E. Shaw Research utilized an in-house peptide-optimization workflow, which combines machine learning with free energy calculations, to predict β-hairpin peptides that would have high affinity for the N-terminal binding site on the nidovirus RdRp-associated nucleotidyltransferase (NiRAN) domain. Their workflow generated a series of amphiphilic peptides that contain both a hydrophobic surface, which in their computational model bound to the NiRAN domain, and a highly positively charged, arginine- and lysine-rich surface that faced away from the NiRAN domain.
![Picture1](https://github.com/rshin1209/desres/assets/25111091/ee9d10e4-506a-4598-af60-1792365ea567)
**Table 1. Amino acid sequences of the peptides are presented.** In the sequence alignment, black font indicates that the residue is present in the wild-type peptide and red indicates that the residue is not present in the wild-type peptide. "Peptide-X-dis" labels denote crosslinked peptides.

# Protein Entropy Calculation
## Python Module Requirements:
* MDtraj
* Networkx
* Numpy
* json
* multiprocessing
## Computational Resource
* CPU: Intel(R) Xeon(R) CPU E5-2420 0 @ 1.90GHz
* Nodes: 1
* Job Wall-clock: 01:49:51
* Memory Usage: 9.29 GB

## Reproduction
Protein Entropy Analysis was performed on six (**peptide-only**) NiRAN domain-targeted peptides. (1) A.pdb. (2) B.pdb. (3) C.pdb. (4) A_dis.pdb. (5) B_dis.pdb. (6) C_dis.pdb.<br />

The main Python script `protein_entropy.py` is used for implementing protein entropy analysis.<br>

Taking the (1) for an example, one can simply run the following two commands to perform protein entropy analysis.<br>

```
python protein_entropy.py --reaction A --temperature 298.15 --job 0
python protein_entropy.py --reaction A --temperature 298.15 --job 1
[reaction] -- Name of pdb file containing simulation trajectory frames
[temperature] -- Temperature (Kelvin) for -TS calculation
[job] -- 0 for degrees of freedom extraction from Cartesian coordinates and 1 for configurational entropy calculation
```

After running the two commands above, you will have 5 output files in the following directory: `./A/`. 
```
topology.txt -- Atom indices of bonds, angles, and torsion angles in the target molecule
topology.pdb -- The first frame of the simulation trajectory file
S1D.npy -- Numpy array of 1D-entropies of individual degrees of freedom (The same order presented in topology.txt
entropy.log -- Protein entropy output containing 1D-entropies, mutual information among DoF pairs, MIE entropy, MIST entropy
A_entropy_map.npy -- Entropy matrix containing individual residue entropy and correlation among residue pairs 
```
# Protein Entropy Analysis
## Entropy Matrix
![A_imshow](https://github.com/rshin1209/desres/assets/25111091/f2c8b976-65a9-464f-9684-58c6ef59e47c)
**Figure 1. Entropy Matrix of Peptide A derived from `./A/A_entropy_map.npy`.** Diagonal elements in the matrix signify the entropy of individual residues in -TS (kcal/mol). A lower value indicates higher entropy for the respective residue. On the other hand, off-diagonal elements denote the correlation between pairs of residues. A higher value suggests a stronger correlation between the paired residues.

## Overall Entropy Comparison
![Picture4](https://github.com/rshin1209/desres/assets/25111091/4a99bca9-7f70-45c6-8117-79b5667459fc)
**Table 2. Protein, residue, and backbone entropy comparison.** In this presentation, I conducted a comparison of protein, residue, and backbone entropy. The values are expressed as -TS (kcal/mol), where lower values correspond to higher entropy. The color ${\color{blue}blue}$ indicates the highest entropy, while ${\color{red}red}$ signifies the lowest entropy. According to the comparison, Peptide C-dis exhibited the highest protein and summed residue entropy, whereas Peptide B demonstrated the lowest protein and summed residue entropy.

## Assessment of Structural Dynamics Change by Mutation
![A_B_imshow](https://github.com/rshin1209/desres/assets/25111091/0ac2a61f-db68-4319-9ab5-ae890d4b4b2f)
**Figure 2. Entropy Matrix Comparison between Peptide A and Peptide B, derived from `./A/A_entropy_map.npy` and `./B/B_entropy_map.npy`, respectively.** 
![A_C_imshow](https://github.com/rshin1209/desres/assets/25111091/d80e28fb-7da3-4449-8cf8-4f24d2ed8c39)
**Figure 3. Entropy Matrix Comparison between Peptide A and Peptide C, derived from `./A/A_entropy_map.npy` and `./C/C_entropy_map.npy`, respectively.** 
![B_C_imshow](https://github.com/rshin1209/desres/assets/25111091/b088f388-d8fc-4d81-97e9-03626060e37c)
**Figure 4. Entropy Matrix Comparison between Peptide B and Peptide C, derived from `./B/B_entropy_map.npy` and `./C/C_entropy_map.npy`, respectively.** 

Traditional energy-based models and quantitative techniques like RMSD and RMSF frequently present limitations in providing a thorough evaluation of structural dynamics alterations induced by mutations. By employing the entropy matrix that I have developed, a more comprehensive assessment of the impact of mutations on structural dynamics can be achieved.

As depicted in **Table 1**, **Figure 2** highlights notable entropy changes in M1Y, I7K, and R20K (A#B denotes mutation from A to B at residue number #). Beyond the mutation sites, phenylalanine (F at residue 15) exhibits a substantial alteration in correlations with other residues. This pattern is further evident in **Figure 3**, comparing peptides A and C, where mutation sites include K2Y, R19L, and R20K. Conversely, the comparison between peptides B and C in **Figure 4** reveals anticipated entropy changes solely at the mutation sites Y1M, K2Y, K7I, and R19L.

![123 copy](https://github.com/rshin1209/desres/assets/25111091/f2c17e74-4837-4001-a5ca-1049d748a91d)
**Figure 5. PhenylAlanine (Resi 15) and Lysine (Resi 20) are presented in Peptides B and C.**

These findings suggest a potential influence of lysine (residue 20), present in peptides B and C but not in peptide A, on restricting the conformational space of phenylalanine—closest to the NiRAN domain (**Figure 5**). This observation aligns with the results of viral infection experiments, where peptide A exhibited the highest average inhibition at 100 and 30 µM. It is conceivable that the mutation of residue 20 from arginine to lysine has led to a loss of correlation between phenylalanine and the NiRAN domain for inhibition, resulting instead in an increased correlation within its own peptide.

# License
This repository is licensed under the MIT License - see the LICENSE.md file for details.
