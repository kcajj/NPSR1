# NPSR1 interface and interaction with NPS

Giacomo Castagnetti, Marianna Gordin, Sara Petrini and Anna Salvadori

## Abstract

The neuropeptide S receptor (NPSR) [1] is a G-protein coupled receptor and it is an integral membrane protein. Neuropeptide S (NPS) [2] is the endogenous ligand of NPSR1. The purpose of this study is to analyze the structure of the receptor and to understand its interaction with NPS. To investigate these aspects, Python scripts were implemented to retrieve general information about the amino acids composition and the solvent-accessible surface area; ChimeraX was used to analyze the tridimensional structure of the receptor and its surface of interaction with the neuropeptide; docking webapps were used to study the interaction with NPS. To further improve the reliability of the docking analysis, a comparison with the oxytocin receptor was performed. Our study allowed us to understand how the receptor interacts with its ligand and to visualize their interaction, also with respect to relevant mutations of NPSR1.

Keywords: NPSR1, NPS, interaction, interface.

## 1. Introduction

NPSR1, short for Neuropeptide S receptor, is a G protein-coupled receptor mainly found in the brain. NPSR1 and the associated ligand NPS, a neuropeptide of about 20 amino acids [3], play an essential role in many physiological and behaviour processes, including anxiety and stress responses, arousal and wakefulness, and emotional and sensory processing [4].
Given its role in many important processes, structural studies are fundamental for a better understanding of its function and its potential for pharmacological intervention. Due to the challenges in the crystallisation process of GPCRs, its structure has not been studied much yet. The available studies on the structural relationship between NPSR1 and NPS have highlighted pathogenic mutations: N107I is involved in asthma pathogenesis [5]; Y206H leads to a decreased need for sleep causing side effects as memory deficit [6].
The aim of this study is to analyse the structure of the receptor, focusing on the interface with NPS, in order to better understand the effects of mutations.

## 2. Materials and Methods

### 2.1 AlphaFold

AlphaFold is an artificial intelligence program developed by DeepMind, a subsidiary of Alphabet, which performs predictions of protein structures. The structure of NPSR1 and NPS were retrieved from AlphaFold database [7]. We chose the NPSR1 AlphaFold structure over the SwissModel [8] one because, in SwissModel, the protein lacks terminal regions (from amino acid 1 to 50 and from amino acid 346 to 371). In order to get consistent results with the analysis performed on the sequence of the protein and on its annotation, we preferred to stick with the AlphaFold model, even though the terminal regions have low prediction accuracy.
Furthermore, the structure of NPSR1 present on SwissModel was built starting from the Crystal structure of the human oxytocin receptor, which has 29.59% identity with NPSR1. The crystal structure of oxytocin presents a lot of missing regions, therefore using it for homology modelling does not allow to have a complete and accurate predicted model. We exploited the crystal structure of the Oxytocin receptor, bound with oxytocin and its G-protein, to compare our predicted binding site for NPS with the binding site for Oxytocin.
The NPS model present in Alphafold is 89 amino acids long. Despite recent papers describe NPS as a 20 aa long peptide [9] , we still analysed the 89 amino amino acid model, since it is reported in the official databases. The AlphaFold structure of NPS has a low general pLDDT score as peptide structure prediction is quite challenging. Nevertheless, AlphaFold is still the best choice when lacking an experimental structure [10].

### 2.2 Python

#### Amino acid composition

Different Python scripts for retrieving the amino acid relative abundance were implemented: using the FASTA file to obtain the general distribution of the amino acids in the protein (functionsFASTA.py); using the GFF file to retrieve the transmembrane domain’s sequence in order to compute amino acid distribution in the domain (functionsGFF.py). In all the analysis the relative abundance was computed as: (aa / TOT) * 100, where: aa is the abundance of the amino acid, and TOT is the total number of amino acids of the considered segment (whole protein or TM domain).

#### SASA

A Python script for retrieving the solvent-accessible surface area (SASA) was implemented. The analysis works by inspecting a DSSP file, retrieved from the PDB file of the considered protein or peptide through the DSSP program (see paragraph of DSSP). The script gives as results the total SASA, the SASA of the transmembrane domain of the receptor, which was obtained through the analysis of NPSR1 with ChimeraX, and the SASA of the interface, computed from the area of the receptor and the neuropeptide in their monomeric and in their complex form. 
Another analysis performed with SASA was to predict which residues of NPSR1 are most involved in the interaction with NPS, therefore pointing out the most relevant residues composing the interface. The analysis has been performed starting from the docking model, considering the residues that had a change in accessibility between monomeric and complex form higher than a specific threshold, specifically 50 Å2.

### 2.3 DSSP

The DSSP program [11,12] computes the most likely secondary structure given the tridimensional structure of a protein, by extracting information from the 3D coordinates. Moreover, a web server (http://www.cmbi.umcn.nl/xssp/) is provided to obtain individual DSSP files for private PDB files. In this study, the web server was used to convert the PDB files of NPS, NPSR and the docking model to DSSP files. The output from DSSP contains secondary structure assignments, residue number, one letter amino acid code and accessibility (residue water exposed surface in Å2), one line per residue.

### 2.4 ISPRED4

ISPRED4 is a web-server for predicting protein-protein interaction sites starting from protein structure [13]. We uploaded AlphaFold's PDB file to the website and ran the algorithm with default RSA  (≥ 0.20) and PDB chain “A”. 

### 2.5 PDB

The Protein Data Bank (PDB) [cite_pdb] is a database of 3D structure data for large biological molecules (proteins, DNA, and RNA). A PDB file is a text file consisting of many lines called records, which can be of different types.

### 2.6 ChimeraX

UCSF ChimeraX is a next-generation molecular visualisation and analysis program developed with the support of National Institutes of Health (NIH). 
We used ChimeraX to generate images that could support and visualise our results, and to perform structural analysis of the NPSR1 AlphaFold model.

### 2.7 ClusPro

ClusPro [14, 15, 16] is a fully automated, web-based program for the computational docking of protein structures. It allows to predict the most likely 3D configuration of two proteins that interact together. The docking run was performed with default options and run on the CPU server of ClusPro. 

## 3. Results

### 3.1 Analysis of amino acid composition

### 3.2 Structural analysis




# Citations

1. Reinscheid RK, Ruzza C. Pharmacology, Physiology and Genetics of the Neuropeptide S System. Pharmaceuticals (Basel). 2021;14(5):401.
2. Okamura N, Reinscheid RK. Neuropeptide S: a novel modulator of stress and arousal. Stress. 2007;10(3):221-6.
3. R.K. Reinscheid. Encyclopedia of Neuroscience; Academic Press: University of California at Irvine, CA, USA, 209; pp. 817-819.
4. Zhang Y, Wang Z, et al. Novel neuropeptides as ligands of orphan G protein-coupled receptors. Curr Pharm Des. 2011;17(25):2626-2631.
5. Bernier V, Stocco R, Bogusky MJ, et al. Structure-function relationships in the neuropeptide S receptor: molecular consequences of the asthma-associated mutation N107I. J Biol Chem. 2006;281(34):24704-24712. 
6. Xing L, Shi G, Mostovoy Y, et al. Mutant neuropeptide S receptor reduces sleep duration with preserved memory consolidation. Sci Transl Med. 2019;11(514):eaax2014.
7. AlphaFold. Available online: https://alphafold.ebi.ac.uk/ (accessed on 17th Dec 2023)
8. SwissModel. Available online: https://swissmodel.expasy.org/ (accessed on 17th Dec 2023)
9. Ruzza C, Calò G, Di Maro S, et al. Neuropeptide S receptor ligands: a patent review (2005-2016). Expert Opin Ther Pat. 2017;27(3):347-362.
10. McDonald EF, Jones T, Plate L, Meiler J, Gulsevin A. Benchmarking AlphaFold2 on peptide structure prediction. Structure. 2023;31(1):111-119.e2.
11. Joosten RP, te Beek TA, Krieger E, et al. A series of PDB related databases for everyday needs. Nucleic Acids Res. 2011;39(Database issue):D411-D419.
12. Kabsch W, Sander C. Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers. 1983;22(12):2577-2637.
13. ISPRED4. Available online: https://ispred4.biocomp.unibo.it (Accessed on 19th Dec 2023)
14. ClusPro. Available online: https://cluspro.bu.edu/home.php (Accessed on 27th Dec 2023)
15. Kozakov D, Hall DR, Xia B, et al. The ClusPro web server for protein-protein docking. Nat Protoc. 2017;12(2):255-278. 
16. Kozakov D, Beglov D, Bohnuud T, et al. How good is automated protein docking?. Proteins. 2013;81(12):2159-2166.