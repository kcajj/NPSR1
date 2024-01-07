# Explanation of the scripts

## Fixing the docking output

The docking model produced by ClusPro had some problems. We manipulated the file through a Python script to correct it and to allow its visualisation and subsequent analyses.

In the original file, both chains were labelled as "A" and the IDs of the initial atoms in the two chains were identical. The script allows to rename the second chain and to make its IDs unique. The original file also presented some extra lines that caused troubles during visualisation, these lines were removed.

By running "fix_docking_output.py" we obtain a perfect pdb file named "docking_model.pdb".

## SASA
The Python script "SASA.py" is a tool to retrieve the solvent-accessible surface area (SASA) of a protein (NPSR1) as a whole, of its transmembrane domain and of its interface (considering the specific ligand NPS). It can be also used on a peptide (NPS).

To run this code download the folder "SASA" as it is, then run on the terminal the following command:
```
python3 SASA.py NPS.dssp NPSR.dssp docking_model.dssp
```

Many results will be printed. Among them, there are the SASA values for NPS and NPSR both in their monomeric forms and in the complex, the surface of interaction and the residues that are predicted to be involved in the interaction.

## Amino acids composition
The python scripts "functionsFASTA.py" and "functionsGFF.py" contain functions to analyse the amino acid relative abundance in different conditions: in the whole protein and transmembrane domain, both by counting each amino acid or by property. The output of the functions is used to create barplots of the relative abundances.