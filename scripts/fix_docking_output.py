if __name__=="__main__":
    input_filename="model.000.09.pdb"
    output_filename="docking_model.pdb"

    b=0
    pdb=open(input_filename,"r")
    new_pdb=open(output_filename,"w")
    for line in pdb:
        #if it is the second chain
        if line=="HEADER lig.000.09.pdb\n":
            b=1
            new_pdb.write(line)
            continue
        if b==0:
            # for the first chain just copy each line
            new_pdb.write(line)

            #take the id:
            piece=line[7:11].strip()
            if piece != '':
                if piece[0] in "1234567890":
                    last_id=int(piece)

            continue
        else:
            last_id+=1
            new_line=line[:7]+str(last_id)+line[11:21]+'B'+line[22:]
            new_pdb.write(new_line)
