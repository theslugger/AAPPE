import random

f = open('../data/Pocket_Reduced_Info.txt','r')
lines1 = f.readlines()
f.close()
lines1.pop(0)

uniprot_pdb_dict = {}
for line in lines1:
    line = line.split('\t')
    if line[1] in uniprot_pdb_dict.keys():
        uniprot_pdb_dict[line[1]].append(line[0])
    else:
        uniprot_pdb_dict[line[1]] = [line[0]]

#print len(uniprot_pdb_dict)
#print uniprot_pdb_dict.items()[:3]


f = open('../outfiles1/General_TargetUniprot_ChEMBL_InavtiveCompound.txt','r')
lines2 = f.readlines()
f.close()

outlines = []
for line in lines2:
    tmp = line.split('\t',1)
    u = tmp[0].split('|')[1]
    outlines.append(random.sample(uniprot_pdb_dict[u],1)[0]+'\t'+tmp[1])


f = open('../outfiles2/PDB_InactiveCompound.txt','w')
f.writelines(random.sample(outlines,3743))
f.close()
