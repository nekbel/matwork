dic={}
flag=1
with open("C:/Users/User/Downloads/mart_export.txt") as f:
    for line in f:
        if flag:
            flag=0
            continue
        k,v = line.strip().split("\t")
        dic[v] = k

#%%
nf = open("C:/Users/User/Downloads/quant_mod.txt", "wt")
flag=1
with open("C:/Users/User/Downloads/quant.sf", "r+") as q:
    for line in q:
        if flag: 
            flag=0
            nf.write(line)
            continue
        data = line.split("\t")
        ens = data[0].split('.')[0]
        data = [dic[ens]] + data
        nf.write("\t".join(data))
nf.close()
#%%
f = open("C:/Users/User/Downloads/quant_mod.txt", "r")

for i in range(0,9):
    print(f.readline())
f.close()
