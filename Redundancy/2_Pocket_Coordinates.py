#The time of a file:0.025s
import os
from tqdm import tqdm

dirs = os.listdir('../refined-set')


for name in tqdm(dirs):
    f = open('../refined-set/'+name+'/'+name+'_protein.pdb','r')
    lines1 = f.readlines()
    f.close()

    f = open('../outfile1/'+name+'/'+name+'_Pocket_AAlist.txt','r')
    lines2 = f.readlines()
    f.close()

    AAlist = [line[:-1] for line in lines2]
    outlines = []
    for line in lines1:
        if line[17:26] in AAlist:
            outlines.append(line)

    f = open('../outfile1/'+name+'/'+name+'_pocket.pdb','w')
    f.writelines(outlines)
    f.close()
