# SASA
The Python script "SASA.py" is a tool to retrieve the solvent-accessible surface area (SASA) of a protein (NPSR1) as a whole, of its transmembrane domain and of its interface (considering the specific ligand NPS). It can be also used on a peptide (NPS).

To run this code download the folder "SASA" as it is, then run on the terminal the following command:
python3 SASA.py NPS.dssp NPSR.dssp docking_model.dssp
The option "-v", added at the end of the previous command, activates the verbose modality.

Many results will be printed. Among them, there are the SASA values for NPS and NPSR both in their monomeric forms and in the complex, the surface of interaction and the residues that are predicted to be involved in the interaction.
