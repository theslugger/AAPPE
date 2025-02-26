import os
import itertools
from tqdm import tqdm


Pocket_5A_FP = ['ProName|']
AA_Without_AllAtoms = []
AAList = ['ARG','LYS','HIS','SER','GLN','THR','CYS','ASN','TYR','ASP',
                 'GLU','PHE','GLY','ALA','LEU','ILE','VAL','PRO','MET','TRP']

AAPairs = []
for i in range(len(AAList)):
    for n in range(i,len(AAList)):
        AAPairs.append(AAList[i]+'_'+AAList[n])

#print len(AAPairs)
#print AAPairs

AAPairs_Index = {}
for i in range(len(AAPairs)):
    for n in range(10):
        AAPairs_Index[AAPairs[i]+'_'+str(n)] = i*10+n
        Pocket_5A_FP.append(AAPairs[i]+'_'+str(n+1)+' ')


#print AAPairs_Index['ASP_ALA_3']
#print Pocket_5A_FP.index('ASP_ALA_3 ')


#print len(AAPairs_Index)
Pocket_5A_FP[-1] = Pocket_5A_FP[-1][:-1]+'\n'

dirs = os.listdir('../refined-set/')
#dirs = ['1pau']

#f = open('FileNameList.txt','w')
#f.writelines('\n'.join(dirs))
#f.close()

n = 0


Pocket_without_5AA = []
for name in tqdm(dirs):
    f = open('../outfile1/'+name+'/'+name+'_pocket.pdb','r')
    lines = f.readlines()
    f.close()

    f = open('../outfile1/'+name+'/'+name+'_Pocket_AAlist.txt','r')
    lines2 = f.readlines()
    f.close()

    if len(lines2) < 5:
        Pocket_without_5AA.append(name)
        continue

    Pocket_5A_FP.append(name+'|')

    AAxyz = {}   ### The non-hydrogen atom which is farthest from the C
    Tmp_His = {} ###center of aromatic ring
    Tmp_Phe = {} ###center of aromatic ring
    Tmp_Pro = {} ###center of aromatic ring
    Tmp_Trp = {} ###center of aromatic ring
    for line in lines:
        if ('CZ  ARG' in line) or ('NZ  LYS' in line) or ('OG  SER' in line) or ('CD GLN' in line) or ('OG1 THR' in line) or ('SG  CYS' in line) or ('CG  ASN' in line) or ('OH  TYR' in line) or ('CA  GLY' in line) or ('CB  ALA' in line) or ('CG  LEU' in line) or ('CD1 ILE' in line) or('CB  VAL' in line) or ('CG  ASP' in line) or ('CD  GLU' in line) or ('CE  MET' in line):
            AAxyz[line[17:26]] = [float(line[30:38].strip()),float(line[38:46].strip()),float(line[46:54].strip())]
        if ('CG  HIS' in line) or ('CD2 HIS' in line) or ('NE2 HIS' in line) or ('CE1 HIS' in line) or ('ND1 HIS' in line):
           if line[17:26] not in Tmp_His.keys():
                Tmp_His[line[17:26]] = [(float(line[30:38].strip()),float(line[38:46].strip()),float(line[46:54].strip()))]
           else:
                Tmp_His[line[17:26]].append((float(line[30:38].strip()),float(line[38:46].strip()),float(line[46:54].strip())))
        if ('CG  PHE' in line) or ('CD1 PHE' in line) or ('CD2 PHE' in line) or ('CE1 PHE' in line) or ('CE2 PHE' in line) or ('CZ  PHE' in line):
            if line[17:26] not in Tmp_Phe.keys():
                Tmp_Phe[line[17:26]] = [(float(line[30:38].strip()),float(line[38:46].strip()),float(line[46:54].strip()))]
            else:
                Tmp_Phe[line[17:26]].append((float(line[30:38].strip()),float(line[38:46].strip()),float(line[46:54].strip())))
        if ('N   PRO' in line) or ('CA  PRO' in line) or ('CB  PRO' in line) or ('CG  PRO' in line) or ('CD  PRO' in line):
            if line[17:26] not in Tmp_Pro.keys():
                Tmp_Pro[line[17:26]] = [(float(line[30:38].strip()),float(line[38:46].strip()),float(line[46:54].strip()))]
            else:
                 Tmp_Pro[line[17:26]].append((float(line[30:38].strip()),float(line[38:46].strip()),float(line[46:54].strip())))
        if ('CD2 TRP' in line) or ('CE2 TRP' in line) or ('CE3 TRP' in line) or ('CZ2 TRP' in line) or ('CZ3 TRP' in line) or ('CH2 TRP' in line):
            if line[17:26] not in Tmp_Trp.keys():
                 Tmp_Trp[line[17:26]] = [(float(line[30:38].strip()),float(line[38:46].strip()),float(line[46:54].strip()))]
            else:
                 Tmp_Trp[line[17:26]].append((float(line[30:38].strip()),float(line[38:46].strip()),float(line[46:54].strip())))

    if Tmp_His != {} :
        for key in Tmp_His:
            if len(Tmp_His[key]) == 5:
                AAxyz[key] = [(Tmp_His[key][0][0] + Tmp_His[key][1][0] +Tmp_His[key][2][0] +Tmp_His[key][3][0] +Tmp_His[key][4][0]) / 5,
                               (Tmp_His[key][0][1] + Tmp_His[key][1][1] +Tmp_His[key][2][1] +Tmp_His[key][3][1] +Tmp_His[key][4][1]) / 5,
                               (Tmp_His[key][0][2] + Tmp_His[key][1][2] +Tmp_His[key][2][2] +Tmp_His[key][3][2] +Tmp_His[key][4][2]) / 5]
            else:
                AA_Without_AllAtoms.append(name+'\t'+key+'\n')
    if Tmp_Pro != {}:
        for key in Tmp_Pro:
            if len(Tmp_Pro[key]) == 5:
                AAxyz[key] = [(Tmp_Pro[key][0][0] + Tmp_Pro[key][1][0] +Tmp_Pro[key][2][0] +Tmp_Pro[key][3][0] +Tmp_Pro[key][4][0]) / 5,
                               (Tmp_Pro[key][0][1] + Tmp_Pro[key][1][1] +Tmp_Pro[key][2][1] +Tmp_Pro[key][3][1] +Tmp_Pro[key][4][1]) / 5,
                               (Tmp_Pro[key][0][2] + Tmp_Pro[key][1][2] +Tmp_Pro[key][2][2] +Tmp_Pro[key][3][2] +Tmp_Pro[key][4][2]) / 5]
            else:
                AA_Without_AllAtoms.append(name+'\t'+key+'\n')
    if Tmp_Phe != {}:
        for key in Tmp_Phe:
            if len(Tmp_Phe[key]) == 6:
                AAxyz[key] = [(Tmp_Phe[key][0][0] + Tmp_Phe[key][1][0] +Tmp_Phe[key][2][0] +Tmp_Phe[key][3][0] +Tmp_Phe[key][4][0] + Tmp_Phe[key][5][0]) / 6,
                               (Tmp_Phe[key][0][1] + Tmp_Phe[key][1][1] +Tmp_Phe[key][2][1] +Tmp_Phe[key][3][1] +Tmp_Phe[key][4][1] + Tmp_Phe[key][5][1]) / 6,
                               (Tmp_Phe[key][0][2] + Tmp_Phe[key][1][2] +Tmp_Phe[key][2][2] +Tmp_Phe[key][3][2] +Tmp_Phe[key][4][2] + Tmp_Phe[key][5][2]) / 6]
            else:
                AA_Without_AllAtoms.append(name+'\t'+key+'\n')
    if Tmp_Trp != {}:
        for key in Tmp_Trp:
            if len(Tmp_Trp[key]) == 6:
                AAxyz[key] = [(Tmp_Trp[key][0][0] + Tmp_Trp[key][1][0] +Tmp_Trp[key][2][0] +Tmp_Trp[key][3][0] +Tmp_Trp[key][4][0] + Tmp_Trp[key][5][0]) / 6,
                               (Tmp_Trp[key][0][1] + Tmp_Trp[key][1][1] +Tmp_Trp[key][2][1] +Tmp_Trp[key][3][1] +Tmp_Trp[key][4][1] + Tmp_Trp[key][5][1]) / 6,
                               (Tmp_Trp[key][0][2] + Tmp_Trp[key][1][2] +Tmp_Trp[key][2][2] +Tmp_Trp[key][3][2] +Tmp_Trp[key][4][2] + Tmp_Trp[key][5][2]) / 6]
            else:
                AA_Without_AllAtoms.append(name+'\t'+key+'\n')
#    print 'len(AAxyz)=',len(AAxyz)
#    print AAxyz

    outlines = []
    for i in range(len(AAPairs_Index)):
        outlines.append(0)

    for (m,n) in list(itertools.combinations(AAxyz.keys(),2)):
        tmp_AAPairs1 = m[:3]+'_'+n[:3]
        tmp_AAPairs2 = n[:3]+'_'+m[:3]
        x = AAxyz[m][0] - AAxyz[n][0]
        y = AAxyz[m][1] - AAxyz[n][1]
        z = AAxyz[m][2] - AAxyz[n][2]
        power = pow(x,2) + pow(y,2) + pow(z,2) ### square of distance


        if tmp_AAPairs1 in AAPairs:
            tmp_AAPairs = tmp_AAPairs1
        if tmp_AAPairs2 in AAPairs:
            tmp_AAPairs = tmp_AAPairs2


        for i in range(10):
            if (power >= pow(2*i,2)) & (power < pow(2*(i+1),2)): ### square of distance
                index = AAPairs_Index[tmp_AAPairs+'_'+str(i)]
                outlines[index] = outlines[index] + 1
#                print tmp_AAPairs,power,tmp_AAPairs+'_'+str(i),index

#    print len(outlines)

#    for i in range(210):
#        print outlines[i*20:20*(i+1)]
    #print name,' The num of %d file has been run.'%(dirs.index(name)+1)

    Pocket_5A_FP.append(' '.join([str(i) for i in outlines]))
    Pocket_5A_FP.append('\n')


f = open('../outfile2/Pocket_5A_FP.txt','w')
f.writelines(Pocket_5A_FP)
f.close()

f = open('../outfile2/AA_Without_AllAtoms.txt','w')
f.writelines(AA_Without_AllAtoms)
f.close()

f = open('../outfile2/Pocket_without_5AA.txt','w')
f.writelines('\n'.join(Pocket_without_5AA)+'\n')
f.close()
