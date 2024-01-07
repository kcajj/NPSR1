def fix_pdb(input_filename, output_filename):
    '''
    takes in input the old file and the new one. It finds where the first chain stops,
    then over writes all the successive lines fixing the chain label and the ID.
    '''

    b=0
    pdb=open(input_filename,"r")
    new_pdb=open(output_filename,"w")

    for line in pdb:

        #if it is the second chain
        if line=="HEADER lig.000.09.pdb\n":
            b=1 #stop skipping the lines
            new_pdb.write(line)
            continue
        
        #skip some extra lines that are present in the file
        if line[:3]=="REM" or line[:3]=="TER":
            continue

        #if it is the first chain
        if b==0:
            #just copy each line
            new_pdb.write(line)

            #take the id of the atom
            piece=line[7:11].strip()
            if piece != '':
                if piece[0] in "1234567890":
                    last_id=int(piece)
            continue

        #if it is the second chain
        else:
            last_id+=1
            new_line=line[:7]+str(last_id)+line[11:21]+'B'+line[22:] #rebuild the line
            new_pdb.write(new_line)

if __name__=="__main__":
    input_filename="structure_models/docking_output.pdb"
    output_filename="structure_models/docking_model.pdb"

    fix_pdb(input_filename, output_filename)
