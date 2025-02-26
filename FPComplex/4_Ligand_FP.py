import pybel
import os
from tqdm import tqdm

dirs = os.listdir('../refined-set')

f = open('../outfile2/Pocket_without_5AA.txt','r')
del_pockets = [i[:-1] for i in f.readlines()]
f.close()

#print 'del_pockets:', del_pockets
dirs = [i for i in dirs if i not in del_pockets]

#print len(dirs)


Ligand_FP = []
for name in tqdm(dirs):
    mol = pybel.readfile('sdf','../refined-set/'+name+'/'+name+'_ligand.sdf').__next__()

    if os.path.exists('../outfile1/'+name+'/') == False:
        os.mkdir('../outfile1/'+name+'/')


    mol.write(format='smi',filename='../outfile1/'+name+'/'+name+'_ligand.smi',overwrite=True)


    fplist = ['0']*1024
    fp = mol.calcfp()
    bits = fp.bits
    for bit in bits:
        fplist[bit-1] = '1'

    Ligand_FP.append(name+'_ligand'+'|'+' '.join(fplist)+'\n')


f = open('../outfile2/Ligand_FP.txt','w')
f.writelines(Ligand_FP)
f.close()
