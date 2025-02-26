import os
import itertools
from tqdm import tqdm

f = open('../data/pdblist1.txt','r')
lines = f.readlines()
f.close()
pdblist = [fn.strip() for fn in lines]
#print 'pdblist: ',len(pdblist)

f = open('../data/INDEX_general_PL_name.2020','r')
lines = f.readlines()
f.close()
lines.sort(key=lambda x:x[6:10],reverse=True)


Uniprot_PDB_Dict = {}
UniprotID = []
for line in lines:
    line = line[:-1].split('  ')
    if line[0] in pdblist:
        tmp = line[2]
        UniprotID.append(line[2])
        if tmp in Uniprot_PDB_Dict.keys():
            Uniprot_PDB_Dict[tmp].append(line[0])
        else:
            Uniprot_PDB_Dict[tmp] = [line[0]]
    else:
        continue
#print(len(Uniprot_PDB_Dict))


Uniprot_PDB_ReduceDict = {}
for key in tqdm(Uniprot_PDB_Dict.keys()):
    PDB_Pocket_Dict = {}
    for pdb in Uniprot_PDB_Dict[key]:
        f = open('/data1/liuwei/targetPre/targetPre/part1/1_DataPreprocess/outfiles1/'+pdb+'/'+pdb+'_Pocket_AAlist.txt')    ##outfiles1 in step1
        AAlines = f.readlines()
        f.close()

        if len(AAlines) >= 5:
            PDB_Pocket_Dict[pdb] = [i[:-1] for i in AAlines]

    if len(PDB_Pocket_Dict) == 0:
        continue

    Uniprot_PDB_ReduceDict[key] = list(PDB_Pocket_Dict.keys())


#print 'Uniprot_PDB_ReduceDict: ',len(Uniprot_PDB_ReduceDict)



outlines1 = []
pdb_num = 0
outlines2 = []

for key in Uniprot_PDB_ReduceDict.keys():
    outlines1.append(key+'\t'+str(len(Uniprot_PDB_ReduceDict[key]))+'\t'+';'.join(Uniprot_PDB_ReduceDict[key])+'\n')
    outlines2.append('\n'.join(Uniprot_PDB_ReduceDict[key])+'\n')
    pdb_num += len(Uniprot_PDB_ReduceDict[key])

#print 'Available PDB Num: ',pdb_num




f = open('../outfiles/Pocket_Reduced_BindingDBInfo.txt','w')
f.writelines(sorted(outlines1))
f.close()


f = open('../outfiles/Pocket_Reduced_List.txt','w')
f.writelines(sorted(outlines2))
f.close()

f = open('../outfiles/UniprotID_List.txt','w')
f.writelines('\n'.join(list(set(UniprotID)))+'\n')
f.close()
