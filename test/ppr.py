l  = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

m = l[:2] + l[2:][::-1]



def isSublist(sublist, larger_list):
    len_sublist = len(sublist)

    for i in range(len(larger_list) - len_sublist + 1):
        if larger_list[i:i + len_sublist] == sublist:
            return True

    return False

n = ['d','e', 'f']

if isSublist(n, l): print('yes')
else: print('no')