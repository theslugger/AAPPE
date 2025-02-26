import joblib

f = open('../data/UniProt_Info2.tab','r')
lines1 = f.readlines()
f.close()
lines1.pop(0)

f = open('../data/INDEX_general_PL_name.2016','r')
lines2 = f.readlines()
f.close()
lines2 = lines2[6:]

f = open('../outfiles/Pocket_Reduced_List.txt','r')
lines3 = f.readlines()
f.close()
pdbs = [line[:-1] for line in lines3]


uniprot_info_dict = {}
for line in lines1:
    line = line[:-1].split('\t')
    if (len(line[1])>0) and (line[1][-1] == ';'):
        line[1] =  line[1][:-1]
    uniprot_info_dict[line[0]] = line

pdb_uniprot_dict = {}
for line in lines2:
    pdb_uniprot_dict[line[:4]] = [line[12:18],line[20:-1]]

outlines = ['PDB ID'+'\t'+'UniProt ID'+'\t'+'Gene name'+'\t'+'Organism'+'\t'+'Protein name'+'\t'+'Entry name'+'\n']
outlines2 = ['PDB ID'+'\t'+'Residues of The Pocket'+'\n']
for pdb in pdbs:
    if pdb in pdb_uniprot_dict:
        uniprot = pdb_uniprot_dict[pdb][0]
        if uniprot in uniprot_info_dict.keys():
            outlines.append(pdb+'\t'+'\t'.join(uniprot_info_dict[uniprot])+'\n')
        else:
            outlines.append(pdb+'\t'+pdb_uniprot_dict[pdb][0]+'\t'+''+'\t'+''+pdb_uniprot_dict[pdb][1]+'\n')

        with open('/data1/liuwei/targetPre/targetPre/part1/1_DataPreprocess/outfiles1/'+pdb+'/'+pdb+'_Pocket_AAlist.txt','r') as f:
            lines = f.readlines()

        outlines2.append(pdb+'\t'+';'.join([i[:-1] for i in lines])+'\n')
    else:
        print(f"PDB: {pdb} not found in pdb_uniprot_dict, skipping.")


f = open('../outfiles/Pocket_Reduced_Info.txt','w')
f.writelines(outlines)
f.close()

joblib.dump(outlines,'../outfiles/Pocket_Reduced_Info.txt.z')


f = open('../outfiles/Pocket_Reduced_Residues.txt','w')
f.writelines(outlines2)
f.close()

joblib.dump(outlines2,'../outfiles/Pocket_Reduced_Residues.txt.z')
