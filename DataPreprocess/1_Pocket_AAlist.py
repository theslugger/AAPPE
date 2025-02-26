import os
import re
from tqdm import tqdm


dirs = os.listdir('../refined-set/')

for name in tqdm(dirs):
    f = open('../refined-set/' + name + '/' + name + '_protein.pdb', 'r')
    lines1 = f.readlines()
    f.close()

    f = open('../refined-set/' + name + '/' + name + '_ligand.mol2', 'r')
    lines2 = f.readlines()
    f.close()

    outlines = []

    m = lines2.index('@<TRIPOS>ATOM\n')
    n = lines2.index('@<TRIPOS>BOND\n')

    for i in range(m + 1, n):
        mol_line = re.sub(r'\s{2,}', ' ', lines2[i])
        tmp = mol_line.split()
        if tmp[1].startswith('H'):
            continue
        else:
            mol_xyz = tmp[2:5]

            for line in lines1:
                if line.startswith('ATOM'):
                    pro_xyz = [line[30:38].strip(), line[38:46].strip(), line[46:54].strip()]
                    x = float(mol_xyz[0]) - float(pro_xyz[0])
                    y = float(mol_xyz[1]) - float(pro_xyz[1])
                    z = float(mol_xyz[2]) - float(pro_xyz[2])
                    power = pow(x, 2) + pow(y, 2) + pow(z, 2)
                    if (power < 25) & (line[17:26] + '\n' not in outlines):
                        outlines.append(line[17:26] + '\n')
                else:
                    continue

    if not os.path.exists('../outfile1/' + name + '/'):
        os.mkdir('../outfile1/' + name + '/')

    f = open('../outfile1/' + name + '/' + name + '_Pocket_AAlist.txt', 'w')
    f.writelines(outlines)
    f.close()
