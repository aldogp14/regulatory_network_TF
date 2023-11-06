import pandas as pd
l = ['r1', 'r2', 'r3', 'r4', 'r2', 'r2']

print(type(l[-1]))

if type(l) == list: print(1)
if type(l[-1]) == str: print(0)