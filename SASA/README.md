# SASA
The Python script "SASA.py" is a tool to retrieve the solvent-accessible surface area (SASA) of a protein (NPSR1) as a whole, of its transmembrane domain and of its interface (considering the specific ligand NPS). It can be also used on a peptide (NPS).

To run this code download the folders and the files as they are, then write on the terminal the following commands, in order:
1. python3 SASA_proj_final.py NPSR/NPSR.dssp A
2. python3 SASA_proj_final.py NPS/NPS.dssp A
3. python3 SASA_proj_final.py NPS-NPSR/NPS-NPSR.dssp A
4. python3 SASA_proj_final.py NPS-NPSR/NPS-NPSR.dssp B
The option "-v", added at the end of the previous commands, allows to retrieve additional information.

Using the results of the total SASA of each of the calls, then call the Python script "SASA_interaction.py" to retrieve the SASA of the surface of interaction between NPS and NPSR1. In particular, run:
python3 SASA_interaction.py SASA_NPS SASA_NPSR SASA_NPSR_INTERACTION SASA_NPS_INTERACTION
