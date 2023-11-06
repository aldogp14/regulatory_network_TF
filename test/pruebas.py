import pandas as pd
l = ['r1', 'r2', 'r3', 'r4', 'r2', 'r2']

c = pd.Series(l).value_counts()
c_s = c[c==1].index
print(c_s)
ins = []
for item_s in c_s:
    print(item_s)
    ins.extend(index for index, item in enumerate(l) if item_s == item and index%2 == 0)
print(ins)

m = []

if m: print('algo')
else: print('vacio')
print(l[l=='r3'].index)
l.insert(0, 4)
print(l)

if not 'hola' in m: print('hey')