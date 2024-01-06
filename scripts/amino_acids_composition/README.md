# Amino acids composition

These python files contain functions to analyse the amino acid relative abundance in different conditions: in the whole protein and transmembrane domain, both by counting each amino acid or by property. The output of the functions is used to create barplots of the relative abundances. 

## Parse FASTA file (_functionsFASTA.py_)
In the function _count_residues_, the FASTA file is traversed line by line to count the occurrence of each amino acid. The relative abundances (count/total number of residues) are saved in a dictionary which is then used in the function _plot_residue_ to plot the relative abundance of each amino acid, colored based on their property. 
The function _count_residues_byproperty_ counts how many amino acids are present in the protein with a specific property. The computed values have been saved in a dictionary then used in _plot_byproperty_ to create the histogram. 
The functions _table_fasta_ and _table_fasta_byproperty_ create a csv file with the outputs of, respectively, _count_residues_ and _count_residues_byproperty_ tab separated. 

## GFF file (_functionsGFF.py_)
The function _retrive_TM_ used the file GFF to retrive the transmembrane domain sequence. Using this information, the function _count_residuesTM_ counts the relative abundance of amino acids in the domain. The output of the function can be plotted using the function _plot_TMvsTOT_. 
The function _plot_comparisonTM_ plots the histogram of the amino acid relative abundance of the TM domain (_count_residuesTM_) compared to the whole sequence (_count_residues_ from _functionsFASTA.py_).
The function _count_residues_bypropertyTM_ counts how many amino acids are present in the transmembrane domain with a specific property. The computed values are  used in _plot_bypropertyTM_ to create the histogram. The function _plot_comparison_bypropertyTM_ creates an histogram of the output of _count_residues_bypropertyTM_ compared to _count_residues_byproperty_ (_functionsFASTA.py_).
The functions _tableTM_ and _table_bypropertyTM_ create a csv file with the outputs of, respectively, _count_residues_ and _count_residues_byproperty_ tab separated. 

## Main 
The functions can be used by giving as parameters the correct file names, 'NPSR1_FASTA.faa' in the functions of _functionsFASTA.py_, 'NPSR1_FASTA.faa' and 'NPSR1_GFF.GFF' for the functions in _functionsGFF.py_. 
