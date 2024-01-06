from pandas import *
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd


def count_residues(FASTA, print_screen = True): 
    ''' Returns a dictionary containing the relative abundance of each possible amino acid
        Parameters: name of the FASTA file to be analysed, and boolean value print_screen to have or not the relative abundances printed on screen'''
    
    # Variable initialization 
    f = open(FASTA, 'r')
    counter = {'G':0,'A':0,'V':0,'P':0, 'I':0, 'L':0,'M':0, 'S':0,  'T':0, 'N':0, 'Q':0, 'C':0, 'F':0, 'W':0, 'Y':0, 'K':0,'R':0,'H':0, 'D':0, 'E':0, "U":0, "O":0, "X":0}
    total = 0

    # Cycle the FASTA file counting each residue and the total number of residues 
    for line in f:
        if ">" not in line:
            for r in line.strip("\n"):
                counter[r] += 1
                total += 1

    # Print of the result and simultanous computation of the frequency
    if print_screen: 
        print("The number of each residue is: ")
        for k in counter.keys():
            counter[k] = round(counter[k]/float(total)*100, 2)
            print(k, counter[k])
    else:
        for k in counter.keys():
            counter[k] = round(counter[k]/float(total)*100, 2)
    return counter

def plot_residue(FASTA):
    ''' Barplot displaying the relative abundances of each amino acid (obtained from the function count_residues()) 
        Parameter: name of the FASTA file to be analysed '''

    # Retrive relative abundances 
    res_count = count_residues(FASTA, print_screen = False)

    # Barplot 
    plt.title("Residue relative abundance")
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

def count_residues_byproperty(FASTA, print_screen=True): 
    ''' Returns a dictionary containing the relative abundances the residues grouped by their properties  
    Parameter: name of the GFF and FASTA file to be analysed, and boolean value print_screen to have or not the relative abundances printed on screen'''

    # Initialization 
    f = open(FASTA, 'r')
    counter = {"NONPOLAR":0, "POLAR":0, "AROMATIC":0, "+CHARGE":0, "-CHARGE":0}
    prop_aa = {"NONPOLAR":['G','A','V','P', 'I', 'L','M'], "POLAR":['S',  'T', 'N', 'Q', 'U', 'C'], "AROMATIC":['F', 'W', 'Y'], "+CHARGE":['K','R','H'], "-CHARGE":['D', 'E']}
    total = 0
    
    # Residue count by property 
    for line in f:
        if ">" not in line:
            for r in line.strip("\n"):
                if r in prop_aa["NONPOLAR"]:
                    counter['NONPOLAR'] += 1
                elif r in prop_aa["POLAR"]:
                    counter['POLAR'] += 1
                elif r in prop_aa["AROMATIC"]:
                    counter['AROMATIC'] += 1
                elif r in prop_aa["+CHARGE"]:
                    counter['+CHARGE'] += 1
                if r in prop_aa["-CHARGE"]:
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

def plot_byproperty(FASTA):
    ''' Barplot displaying the relative abundances of each amino acid collapsed by the property (obtained from the function count_residues()) 
        Parameter: name of the FASTA file to be analysed '''
    
    # Retrive dictionary 
    res_count = count_residues_byproperty(FASTA, False)

    # Barplot 
    plt.title("Residue relative abundance")
    plt.bar(res_count.keys(), res_count.values(), color = ['orange', 'green', 'pink', 'red', 'blue'])
    plt.xlabel("Properties")
    plt.ylabel("Relative abundance (%)")

    # Legend 
    patches =[mpatches.Patch(color='orange', label='Non polar'), mpatches.Patch(color='green', label='Polar'), mpatches.Patch(color='pink', label='Aromatic'), mpatches.Patch(color='red', label='Positive charge'), mpatches.Patch(color='blue', label='Negative charge'), mpatches.Patch(color='black', label='Other')]
    plt.legend(handles=patches, prop={'size': 8})
        
    plt.show()

def table_fasta(FASTA): 
    ''' Print in the file "TablesFASTA1.csv" the residues count tab separated
    Parameters: FASTA file to be analysed'''

    # Initialization 
    res_count = count_residues(FASTA, False)
    d = {}

    # Create a dictionary compatible with the pandas dataframe 
    for k in res_count.keys():
        d[k] = [res_count[k]]

    # Create the pandas dataframe 
    df = pd.DataFrame.from_dict(d, orient = 'index', columns=['Relative abundances(%)'])

    # Transform the dataframe in CSV file 
    df.to_csv('TablesFASTA1.csv', sep ='\t') 

def table_fasta_byproperty(FASTA): 
    ''' Print in the file "TablesFASTA2.csv" the residues count by property tab separated
    Parameters: FASTA file to be analysed'''

    # Initialization 
    res_count = count_residues_byproperty(FASTA, False)
    d = {}

    # Create a dictionary compatible with the pandas dataframe 
    for k in res_count.keys():
        d[k] = [res_count[k]]

    # Create the pandas dataframe
    df = pd.DataFrame.from_dict(d, orient = 'index', columns=['Relative abundances(%)'])

    # Transform the dataframe in CSV file 
    df.to_csv('TablesFASTA2.csv', sep ='\t') 