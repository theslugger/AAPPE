import os
f = open('../data/Pocket_Reduced_Info.txt','r')
lines = f.readlines()
f.close()
lines.pop(0)


UniprotID_list = []
for line in lines:
    line = line[:-1].split('\t')
    UniprotID_list.append(line[1])



UniprotID_WithoutRepeat= list(set(UniprotID_list))
#print len(UniprotID_WithoutRepeat)

f = open('../outfiles1/INDEX_General_UniprotName_List.txt','w')
f.writelines('\n'.join(UniprotID_WithoutRepeat)+'\n')
f.close()
