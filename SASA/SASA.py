import sys

def parse_dssp(file, chain):
    ''' Parse the DSSP file and save amino acid, accessibility and residue number columns for each residue.
        Parameters: name of the DSSP file, chain on which SASA is computed. '''
    
    # Initialise list, state variable and amino acids count 
    dssp_list = []
    s = 0
    SA_count = 0 # solvent accessible aa count

    # Read the DSSP file and read each line
    f = open(file)
    for line in f:
        line = line.rstrip()

        # If reached the line "  #  RESIDUE"
        if line.find("  #  RESIDUE") > -1:
            s = 1
            continue

        # Check if reached the line "  #  RESIDUE" and if the chain is correct 
        if s == 0 or line[11] != chain:
            continue   

        # Read the variables (aa, acc, rn)
        # Amino acid
        aa = line[13]  
        if aa.islower(): # SS-bridge CYS is indicated with lower case letters   
            aa = 'C' # convert to unbridged cystine   
        if aa == '!':
            continue
        # Accessbility
        acc = float(line[34:38])
        if acc > 35: # threshold to find the amino acids that are exposed with an area > 35 angstrom^2
            SA_count += 1
        # Residue number
        rn = int(line[5:10])
        
        # Add the variables to the list
        aa_dssp = [aa, acc, rn]
        dssp_list.append(aa_dssp)

    return dssp_list, SA_count


def retrieve_interface(filename): 
    ''' Reads a rtf or txt file containing the H-bonds involved in the interaction between the ligand (NPS, chain B) and the protein (NPSR, chain A)
        Parameter: name of the rtf or txt file '''
    
    # Initialise the list
    interface_list = []

    # Read the rtf or txt file and read each line
    file = open(filename)
    for line in file:
        line = line.rstrip()
        if line == "": # skip empty lines
            continue

        # Read lines containing the H-bonds records 
        if line[0] == "/":
            donor = (line[1], int(line[7:10].strip())) # donor chain and coordinate
            acceptor = (line[17], int(line[23:26].strip())) # acceptor chain and coordinate
            
            # Flag to easily retrieve the coordinate in the function "get_interface_acc"
            flag = 0
            if "A" in acceptor:
                flag = 1
            
            # Add the variables to the list
            interface_list.append([donor, acceptor, flag])

    return interface_list


def get_total_SASA(dssp_list):
    ''' Returns the total solvent-accessible surface area (SASA) of the protein/peptide.
        Parameter: the list containing amino acid, accessiblity and residue number (derived from the function "parse_dssp) '''
    
    # Initialisation of the accessibilities list
    accs = []

    # Add the accessibility of each residue to the accessibilities list
    for aa in dssp_list:
        accs.append(aa[1])
    
    return sum(accs) # return the sum of accessibilities


def get_transmembrane_SASA(dssp_list, tm_list):
    ''' Returns the total solvent-accessible surface area (SASA) of the transmembrane domain of the protein.
        Parameter: the list containing amino acid, accessiblity and residue number (derived from the function "parse_dssp), and the list of the TM residues (derived from ChimeraX) '''
    
    # Initialisation of the accessibilities list
    accs = []

    # Create a new TM list, containing all the coordinates of the TM residues
    new_tm_list = []
    for element in tm_list:
        for i in range(element[0], element[1]):
            new_tm_list.append(i)
    
    # Add the accessibility only of residues belonging to the TM list, to the accessibilities list
    for aa in dssp_list:
        if aa[2] in new_tm_list:
            accs.append(aa[1])

    return sum(accs) # return the sum of accessibilities

 
def get_interface_SASA(dssp_list, interface_list):
    ''' Returns the total solvent-accessible surface area (SASA) of the interface of the protein, given its binding with the ligand.
        Parameter: the list containing amino acid, accessiblity and residue number (derived from the function "parse_dssp), and the list of the interface residues (derived from ChimeraX using the docking model) '''

    # Initialisation of the accessibilities list
    accs = []

    # Create a new interaction list, containing the unique coordinates of the residues involved in binding
    aa_inter_list = []
    for interaction in interface_list:
        aa_number = interaction[interaction[2]][1] # choose chain A (receptor) interface residue coordinates according to the flag added in the function "retrieve_interface"
        aa_inter_list.append(aa_number)
    original_list = aa_inter_list.copy() # copy the original list
    print("Total number of H bonds:", len(aa_inter_list))
    aa_inter_list = list(set(aa_inter_list)) # remove duplicates
    aa_inter_list.sort() # sort list 

    if verbose == "-v":
        print("Coordinates of residues involved in the interaction:", aa_inter_list)
        print("Number of residues involved in the interaction:", len(aa_inter_list))
    
    # Add the accessibility only of residues belonging to the interface list, to the accessibilities list
    for coordinate in aa_inter_list:
        aa = dssp_list[coordinate-1]
        accs.append(aa[1])
        
        if verbose == "-v":
            print(f"\tThe residue {aa[0]} with coordinate {coordinate} is involved in interaction through {original_list.count(coordinate)} H-bond(s)")
    
    return sum(accs) # return the sum of accessibilities


def interaction_surface(s_nps, s_npsr, s_npsr_inter, s_nps_inter):
    ''' Computes the surface of interaction between NPS and NPSR1
        Parameters: the values of SASA retrieved with the four specificed calls of "SASA.py" '''
    return ( float(s_nps) + float(s_npsr) - ( float(s_npsr_inter) + float(s_nps_inter) ) ) / 2


def find_interface_residues(dsspfile_monomeric, dsspfile_complex):
    ''' Finds the residues at the interface by comparing the dssp files of chains in the monomeric and complex forms and keeping only residues with difference in accessibility greater than 50.
        Parameters: dssp file of the monomeric protein, dssp file of the complex. '''
    # Parse dssp files
    dssp_monomeric = parse_dssp(dsspfile_monomeric, "A")[0]
    dssp_complex = parse_dssp(dsspfile_complex, "A")[0]

    # Keep only residues with a difference in accessibility greater than 50
    interface_residues = []
    for i in range(len(dssp_monomeric)):
        diff = dssp_monomeric[i][1] - dssp_complex[i][1]
        diff = diff if diff >= 0 else -diff
        if diff > 50:
            interface_residues.append( ( dssp_monomeric[i][0], dssp_monomeric[i][2] ) )
    return interface_residues


if __name__ == '__main__':

    # Read command line arguments

    # Input: 3 files
    file_NPS = sys.argv[1]
    file_NPSR = sys.argv[2]
    file_docking = sys.argv[3]
    if len(sys. argv) == 5:
        verbose = sys.argv[4]
    else:
        verbose = ""

    # Create the inputs
    inputs = [
        (file_NPS, "A"),
        (file_NPSR, "A"),
        (file_docking, "A"), # receptor
        (file_docking, "B") # NPS
        ]

    # Initialise the SASA list (dict)
    SASA_list = {}

    # Run the analysis on each input
    for inp in inputs:
        file, chain = inp[0], inp[1]
        print("\n• Analysis on:", file, chain)

        dssp_list, aa_count = parse_dssp(file, chain)[0], parse_dssp(file, chain)[1]

        # 1. Amino acids counts
        print("\nAMINO ACIDS COUNTS")
        print("Total amino acids:", len(dssp_list))
        print("Exposed amino acids:", aa_count)

        # 2. SASA analysis
        print("\nSASA ANALYSIS")

        # 2.1. Total SASA
        total_SASA = get_total_SASA(dssp_list)
        SASA_list[(file, chain)] = total_SASA
        print(f"Total SASA: {total_SASA} Å^2")
        
        # 2.2. Transmembrane domain SASA
        if file == "NPSR.dssp":
            tm_coordinates = [
                (53,73),
                (83,103),
                (124,144),
                (165,185),
                (213,233),
                (276,296),
                (313,333)    ]
                # the coordinates were obtained with ChimeraX
            tm_SASA = get_transmembrane_SASA(dssp_list, tm_coordinates)
            print(f"Transmembrane SASA for NPSR: {tm_SASA} Å^2")

        # 3. Interface and interaction analysis
        if file == "NPSR.dssp":
            print("\nINTERFACE ANALYSIS")
            interface_list = retrieve_interface("interface.rtf")
                # the file "interface.rtf" was obtained with ChimeraX using the docking model
            inter_SASA = get_interface_SASA(dssp_list, interface_list)
            print("Interface SASA for NPSR: ", inter_SASA)
    
    print("")

    # Compute the surface of interaction
    print("SURFACE OF INTERACTION ANALYSIS")
    s_nps, s_npsr, s_npsr_inter, s_nps_inter = "","","",""
    SASA_results = [s_nps, s_npsr, s_npsr_inter, s_nps_inter]
    for i in range(len(SASA_list.values())):
        SASA_results[i] = list(SASA_list.values())[i]
    s_nps, s_npsr, s_npsr_inter, s_nps_inter = SASA_results[0], SASA_results[1], SASA_results[2], SASA_results[3]

    surface_interaction = interaction_surface(s_nps, s_npsr, s_npsr_inter, s_nps_inter)
    print("Surface of interaction:", surface_interaction, "Å^2\n")

    # Interface residues
    inter_res = find_interface_residues(file_NPSR, file_docking)
    print("The residues involved in the interaction are:")
    for res in inter_res:
        print(res)
    print("")

