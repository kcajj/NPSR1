import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from functionsFASTA import *

def retrive_TM (GFF, FASTA):
    ''' Retrive from the FASTA file, using the GFF file, the sequence of the transmembrane domain 
        Parameters: name of the FASTA file to be analysed and the name of the corresponding GFF file '''
    
    # Create GFF dataframe from the file 
    df = pd.read_csv(
        GFF,
        sep="\t",
        comment="#",
        names=[
            "seq_id",
            "source",
            "type",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
            ],
        index_col=False
        )
    f = open(FASTA, 'r')
    sequence = ''

    # Save the whole sequence in one variable 
    for line in f:
          if ">" not in line: sequence += line.strip("\n")

    # Retrive the TM sequence using the coordinates from the GFF file 
    TMseq = []
    for i in range(len(df)):
          if df['type'][i] == 'Transmembrane':
                start = df['start'][i]-1
                end = df['end'][i]
                TMseq.append(sequence[start: end])
                # -1 from both 'start' and 'end' because the coordinates in the GFF file are base 1 
                # +1 in the 'end' coordinates because the python splicing ends at f-1 (a[i:f]) 
                # hence just -1 from the start

    return TMseq

def count_residuesTM(GFF, FASTA, print_screen = True): 
    ''' Returns a dictionary containing the relative abundance of each possible amino acid in the transmembrane domain 
        Parameters: name of the FASTA file to be analysed and the name of the corresponding GFF file, and boolean value print_screen to have or not the relative abundances printed on screen'''
    
    # Retrive transmembrane domain sequence 
    TMseq = retrive_TM (GFF, FASTA)

    # Initialization 
    counter = {'G':0,'A':0,'V':0,'P':0, 'I':0, 'L':0,'M':0, 'S':0,  'T':0, 'N':0, 'Q':0, 'C':0, 'F':0, 'W':0, 'Y':0, 'K':0,'R':0,'H':0, 'D':0, 'E':0, "U":0, "O":0, "X":0}
    total = 0

    # Count the residue occurence in the TM sequence 
    for seq in TMseq:
        for r in seq:
            counter[r] += 1
            total += 1

    # Print and frequency computation 
    if print_screen:
        print("The number of each residue is: ")
        for k in counter.keys():
            counter[k] = round(counter[k]/float(total)*100, 2) 
            print(k, counter[k])
    else:
        for k in counter.keys():
            counter[k] = round(counter[k]/float(total)*100, 2)  

    return counter

def plot_TM(GFF, FASTA):
    ''' Barplot displaying the relative abundances of each amino acid in the transmembrane domain 
        Parameter: name of the FASTA file to be analysed and the name of the corresponding GFF file'''
    
    # Retrive amino acid abundance 
    res_count = count_residuesTM(GFF, FASTA, False)

    # Barplot 
    plt.title("Residue relative abundance in the transmembrane domain")
    barlist = plt.bar(res_count.keys(), res_count.values())
    plt.xlabel("Residue")
    plt.ylabel("Relative abundance (%)")

    # Legend 
    patches =[mpatches.Patch(color='orange', label='Non polar'), mpatches.Patch(color='green', label='Polar'), mpatches.Patch(color='pink', label='Aromatic'), mpatches.Patch(color='red', label='Positive charge'), mpatches.Patch(color='blue', label='Negative charge'), mpatches.Patch(color='black', label='Other')]
    plt.legend(handles=patches, prop={'size': 8})

    # Color each bar based on the property of the amino acid 
    dict_values = list(res_count.keys())
    prop_aa = {"NONPOLAR":['G','A','V','P', 'I', 'L','M'], "POLAR":['S',  'T', 'N', 'Q', 'U', 'C'], "AROMATIC":['F', 'W', 'Y'], "+CHARGE":['K','R','H'], "-CHARGE":['D', 'E']}
    for i in range(0, len(dict_values)):
        if dict_values[i] in prop_aa['NONPOLAR']:
            barlist[i].set_color('orange')
        elif dict_values[i] in prop_aa['POLAR']:
            barlist[i].set_color('green')
        elif dict_values[i] in prop_aa['AROMATIC']:
            barlist[i].set_color('pink')
        elif dict_values[i] in prop_aa['+CHARGE']:
            barlist[i].set_color('red')
        elif dict_values[i] in prop_aa['-CHARGE']:
            barlist[i].set_color('blue')
        else:
            barlist[i].set_color('black')
        
    plt.show()

def plot_comparisonTM (GFF, FASTA):
    ''' Barplot displaying the relative abundances of each amino acid in the transmembrane domain compared to the whole sequence 
        Parameter: name of the FASTA file to be analysed and the name of the corresponding GFF file '''
    
    # Retrive the amino acid abundance in the TM domain and whole protein 
    res_countTM = count_residuesTM(GFF, FASTA, False)
    res_countTOT = count_residues(FASTA, False)

    # Set position of bar on X axis 
    barWidth = 0.25
    br1 = np.arange(len(res_countTOT))
    br2 = [x + barWidth for x in br1] 
    br_legend = [x + barWidth/2 for x in br1] 
    
    # Barplots
    barlistTOT = plt.bar(  br1, res_countTOT.values(), width = barWidth)
    barlistTM = plt.bar(br2, res_countTM.values(), width = barWidth)
    plt.title("Whole sequence and transmembrane domain")
    plt.xlabel("Residue")
    plt.ylabel("Relative abundance (%)") 
    plt.xticks(br_legend, res_countTM.keys(), fontsize=10)

    # Legend 
    patches =[mpatches.Patch(color='orange', label='Non polar, TOT'), mpatches.Patch(color='#ff681f', label='Non polar,TM'), 
              mpatches.Patch(color='green', label='Polar, TOT'), mpatches.Patch(color='#78b992', label='Polar, TM'),
              mpatches.Patch(color='pink', label='Aromatic, TOT'),  mpatches.Patch(color='#de70a1', label='Aromatic, TM'),
              mpatches.Patch(color='red', label='Positive charge, TOT'), mpatches.Patch(color='#b22222', label='Positive charge, TM'),
              mpatches.Patch(color='blue', label='Negative charge, TOT'), mpatches.Patch(color='#362a6e', label='Negative charge, TM'),
              mpatches.Patch(color='black', label='Other, TOT'), mpatches.Patch(color='gray', label='Other, TM')]
    plt.legend(handles=patches, prop={'size': 6})

    # Color each bar based on the property of the amino acid  
    dict_values = list(res_countTOT.keys())
    prop_aa = {"NONPOLAR":['G','A','V','P', 'I', 'L','M'], "POLAR":['S',  'T', 'N', 'Q', 'U', 'C'], "AROMATIC":['F', 'W', 'Y'], "+CHARGE":['K','R','H'], "-CHARGE":['D', 'E']}
    for i in range(0, len(dict_values)):
        if dict_values[i][:1] in prop_aa['NONPOLAR']:
            barlistTOT[i].set_color('orange')
            barlistTM[i].set_color('#ff681f')
            
        elif dict_values[i][:1] in prop_aa['POLAR']:
            barlistTOT[i].set_color('green')
            barlistTM[i].set_color('#78b992')
            
        elif dict_values[i][:1] in prop_aa['AROMATIC']:
            barlistTOT[i].set_color('pink')
            barlistTM[i].set_color('#de70a1')

        elif dict_values[i][:1] in prop_aa['+CHARGE']:
            barlistTOT[i].set_color('red')
            barlistTM[i].set_color('#b22222')

        elif dict_values[i][:1] in prop_aa['-CHARGE']:
            barlistTOT[i].set_color('blue')
            barlistTM[i].set_color('#362a6e')

        else:
            barlistTOT[i].set_color('black')
            barlistTM[i].set_color('gray') 
    

    plt.show()

def count_residues_bypropertyTM(GFF, FASTA, print_screen = True): 
    ''' Returns a dictionary containing the relative abundances the residues in the TM domain grouped by their properties  
    Parameter: name of the GFF and FASTA file to be analysed, and boolean value print to have or not the relative abundances printed on screen'''

    # Initialization 
    TMseq = retrive_TM (GFF, FASTA)
    counter = {"NONPOLAR":0, "POLAR":0, "AROMATIC":0, "+CHARGE":0, "-CHARGE":0}
    prop_aa = {"NONPOLAR":['G','A','V','P', 'I', 'L','M'], "POLAR":['S',  'T', 'N', 'Q', 'U', 'C'], "AROMATIC":['F', 'W', 'Y'], "+CHARGE":['K','R','H'], "-CHARGE":['D', 'E']}
    total = 0
    
    # Residue count by property 
    for seq in TMseq:
        for r in seq:
            if r in prop_aa["NONPOLAR"]:
                counter['NONPOLAR'] += 1
            elif r in prop_aa["POLAR"]:
                counter['POLAR'] += 1
            elif r in prop_aa["AROMATIC"]:
                counter['AROMATIC'] += 1
            elif r in prop_aa["+CHARGE"]:
                counter['+CHARGE'] += 1
            elif r in prop_aa["-CHARGE"]:
                counter['-CHARGE'] += 1
            total += 1

    # Print and frequency calculation 
    if print_screen: 
        print("The number of each residue is: ")
        for k in counter.keys():
            counter[k] = round(counter[k]/float(total)*100, 2) 
            print(k, counter[k])
    else: 
        for k in counter.keys():
            counter[k] = round(counter[k]/float(total)*100, 2)
    return counter

def plot_bypropertyTM(GFF, FASTA):
    ''' Barplot displaying the relative abundances of each amino acid in the transmembrane domain collapsed by the property (obtained from the function count_residues()) 
        Parameter: name of the GFF and FASTA file to be analysed '''
    
    # Retrive dictionary 
    res_count = count_residues_bypropertyTM(GFF, FASTA, False)

    # Barplot 
    plt.title("Residue relative abundance transmembrane domain by property")
    plt.bar(res_count.keys(), res_count.values(), color = ['orange', 'green', 'pink', 'red', 'blue'])
    plt.xlabel("Properties")
    plt.ylabel("Relative abundance (%)")

    # Legend 
    patches =[mpatches.Patch(color='orange', label='Non polar'), mpatches.Patch(color='green', label='Polar'), mpatches.Patch(color='pink', label='Aromatic'), mpatches.Patch(color='red', label='Positive charge'), mpatches.Patch(color='blue', label='Negative charge'), mpatches.Patch(color='black', label='Other')]
    plt.legend(handles=patches, prop={'size': 8})
        
    plt.show()

def plot_comparison_bypropertyTM (GFF, FASTA):
    ''' Barplot displaying the relative abundances of each amino acid in the transmembrane domain compared to the whole sequence 
        Parameter: name of the FASTA file to be analysed and the name of the corresponding GFF file '''
    
    # Retrive the amino acid abundance in the TM domain and whole protein 
    res_countTM = count_residues_bypropertyTM(GFF, FASTA, False)
    res_countTOT = count_residues_byproperty(FASTA, False)

    # Set position of bar on X axis 
    barWidth = 0.25
    br1 = np.arange(len(res_countTOT))
    br2 = [x + barWidth for x in br1] 
    br_legend = [x + barWidth/2 for x in br1] 
    
    # Barplots
    plt.bar( br1, res_countTOT.values(), width = barWidth, color = ['orange', 'green', 'pink', 'red', 'blue', 'black'])
    plt.bar(br2, res_countTM.values(), width = barWidth, color = ['#ff681f', '#78b992', '#de70a1', '#b22222', '#362a6e', 'gray'])
    plt.title("Whole sequence and transmembrane domain by property")
    plt.xlabel("Residue")
    plt.ylabel("Relative abundance (%)") 
    plt.xticks(br_legend, res_countTM.keys(), fontsize=10)

    # Legend 
    patches =[mpatches.Patch(color='orange', label='Non polar, TOT'), mpatches.Patch(color='#ff681f', label='Non polar,TM'), 
              mpatches.Patch(color='green', label='Polar, TOT'), mpatches.Patch(color='#78b992', label='Polar, TM'),
              mpatches.Patch(color='pink', label='Aromatic, TOT'),  mpatches.Patch(color='#de70a1', label='Aromatic, TM'),
              mpatches.Patch(color='red', label='Positive charge, TOT'), mpatches.Patch(color='#b22222', label='Positive charge, TM'),
              mpatches.Patch(color='blue', label='Negative charge, TOT'), mpatches.Patch(color='#362a6e', label='Negative charge, TM'),
              mpatches.Patch(color='black', label='Other, TOT'), mpatches.Patch(color='gray', label='Other, TM')]
    plt.legend(handles=patches, prop={'size': 6})

    
    plt.show()

def tableTM(GFF, FASTA): 
    ''' Print in the file "TablesGFF1.csv" the residues count tab separated
    Parameters: FASTA file to be analysed and corresponding GFF file'''

    # Initialization 
    res_countTM = count_residuesTM(GFF, FASTA, False)
    res_countTOT = count_residues(FASTA, False)
    d = {}

    # Create a dictionary compatible with the pandas dataframe 
    for k in res_countTM.keys():
        d[k] = [ res_countTOT[k], res_countTM[k], abs(round(res_countTOT[k]-res_countTM[k], 2))]

    # Create the pandas dataframe
    df = pd.DataFrame.from_dict(d, orient='index', columns = ['Whole protein', 'Transmembrane domain', 'Difference (absolute value)'])

    # Transform the dataframe in CSV file
    df.to_csv('TablesGFF1.csv', sep ='\t') 

def table_bypropertyTM(GFF, FASTA): 
    ''' Print in the file "TablesGFF1.csv" the residues count by property tab separated
    Parameters: FASTA file to be analysed and corresponding GFF file'''

    # Initialization 
    res_countTM = count_residues_bypropertyTM(GFF, FASTA, False)
    res_countTOT = count_residues_byproperty(FASTA, False)
    d = {}

    # Create a dictionary compatible with the pandas dataframe 
    for k in res_countTM.keys():
        d[k] = [res_countTOT[k], res_countTM[k], abs(round(res_countTOT[k]-res_countTM[k], 2))]

    # Create the pandas dataframe
    df = pd.DataFrame.from_dict(d, orient='index', columns = ['Whole protein', 'Transmembrane domain', 'Difference (absolute value)'])

    # Transform the dataframe in CSV file
    df.to_csv('TablesGFF2.csv', sep ='\t') 
