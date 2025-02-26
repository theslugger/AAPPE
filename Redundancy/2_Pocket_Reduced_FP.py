import joblib


f = open('../data/Pocket_5A_FP.txt','r')
lines1 = f.readlines()
f.close()
lines1.pop(0)

f = open('../outfiles/pdblist.txt','r')
lines2 = f.readlines()
f.close()
lines2 = [line[:-1] for line in lines2]


pocket_fp_dict = {}
for line in lines1:
    line = line.split('|')
    pocket_fp_dict[line[0]] = line[1]


outlines = []
for pdb in lines2:
    if pdb in pocket_fp_dict:
        outlines.append(pdb + '|' + pocket_fp_dict[pdb])
    else:
        print(f"Key {pdb} not found in pocket_fp_dict, skipping.")

f = open('../outfiles/Pocket_Reduced_FP.txt','w')
f.writelines(outlines)
f.close()

joblib.dump(outlines,'../outfiles/Pocket_Reduced_FP.txt.z')
