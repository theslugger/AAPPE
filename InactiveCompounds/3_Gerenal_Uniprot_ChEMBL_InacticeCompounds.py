import random

f = open('../outfiles1/General_Uniprot_ChEMBL_ID.txt','r')
lines1 = f.readlines()
f.close()

f = open('../data/chembl_target_CHEMBLID_inactive.tab','r')
lines2 = f.readlines()
f.close()

f =open('../data/Pocket_Reduced_Info.txt','r')
lines3 = f.readlines()
f.close()



###Targets_Componud_Dict={'target':'compounds'} ChEMBL Database
Targets_Componud_Dict = {}
ChEMBL_Uniprot_Dict = {}
for line in lines1:
    line = line.split('\t')
    Targets_Componud_Dict[line[1][:-1]] = []
    ChEMBL_Uniprot_Dict[line[1][:-1]] = line[0]


###
for line in lines2:
    tmp = line.split('\t')
    if tmp[-1][:-1] in Targets_Componud_Dict.keys():
        Targets_Componud_Dict[tmp[-1][:-1]].append(tmp[1]+'\t'+tmp[2])

#print len(Targets_Componud_Dict)
#print Targets_Componud_Dict.items()[:3]




outlines = []

Targets_With_Compounds = ['TargetName','\t','Compound Number','\n']


for key in Targets_Componud_Dict.keys():
    if Targets_Componud_Dict[key] != []:
        Targets_With_Compounds.append(key+'|'+ChEMBL_Uniprot_Dict[key]+'\t'+str(len(Targets_Componud_Dict[key]))+'\n')
        if len(Targets_Componud_Dict[key]) < 7:
            for value in Targets_Componud_Dict[key]:
                outlines.append(key+'|'+ChEMBL_Uniprot_Dict[key]+'\t'+value+'\n')
        else:
            for value in random.sample(Targets_Componud_Dict[key],7):
                outlines.append(key+'|'+ChEMBL_Uniprot_Dict[key]+'\t'+value+'\n')



###


#print 'The number of targets with inactive compounds: ',len(Targets_With_Compounds)-4




f = open('../outfiles1/General_TargetChEMBL_With_InactiveCompounds_NumList.txt','w')
f.writelines(Targets_With_Compounds)
f.close()


f = open('../outfiles1/General_TargetUniprot_ChEMBL_InavtiveCompound.txt','w')
f.writelines(outlines)
f.close()
