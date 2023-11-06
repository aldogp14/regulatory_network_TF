import test as t
import pandas as pd  

lines_pathway_rxns = t.readTxts('test/test.txt')
df_pathway_rxns = t.makeDF(lines_pathway_rxns, 0)
print(df_pathway_rxns)

row = 0

list_all_pathways = []
pathway_names = [df_pathway_rxns.iloc[0,0]]
p_memory = []

def getListPathways():
    # el for sirve para recorrer las filas de 'df_pathways_rxns', va otorgando los indices
    for row in range(df_pathway_rxns.shape[0]):
        # este if sirve para conocer cuando existe un cambio de pathway, cuando ya no es la misma con la que se trabajo en la iteracion pasada
        if df_pathway_rxns.iloc[row,0] != pathway_names[-1]:
            list_all_pathways.append(p_memory) # guardar la via ya terminada en 'pathways'
            p_memory = [] # resetear la lista donde se guarda la via actual/con la que se trabaja
            pathway_names.append(df_pathway_rxns.iloc[row,0]) # actualizar el nombre de la via
        # esto se hace para cada iteracion, guardar las reacciones de la fila en cuestion
        p_memory.append(df_pathway_rxns.iloc[row, 1])
        p_memory.append(df_pathway_rxns.iloc[row, 2])
        # saber si es la ultima iteracion de todo el for, la ultima fila. No puede ir junto con el primer if porque sino no se guardan las dos ultimas reactions de esta ultima via
        if row == df_pathway_rxns.shape[0]-1: list_all_pathways.append(p_memory)
    return(list_all_pathways)


def getTypePath(current_p):
    begins = []
    ends = []
    if max(pd.Series(current_p).value_counts()) == 2 and min(pd.Series(current_p).value_counts()) == 2: 
        return('cycle', ends, begins)
    # sino fue un ciclo simple entonces hay que discernir si es 'lineal' o empieza como ciclo y luego como linea
    counts = pd.Series(current_p).value_counts()
    uniques = counts[counts==1].index
    for item_s in uniques:
        begins.extend(index for index, item in enumerate(current_p) if item_s == item and index%2 == 0)
        ends.extend(index for index, item in enumerate(current_p) if item_s == item and index%2 == 1)
    
    not_uniques = [item for item in current_p if current_p.count(item)>1]
    for item_s in list(set(not_uniques)):
        c_begin = 0
        c_end = 0
        indexes = [index for index, item in enumerate(current_p) if item_s == item]
        for index in indexes:
            if index%2 == 1 and c_begin == 0: c_end += 1
            elif index%2 ==0 and c_end == 0: c_begin += 1
            else: break
            if c_end == len(indexes): ends.extend(indexes)
            elif c_begin == len(indexes): begins.extend(indexes)        

    # si tiene inicios entonces es una via normal o lineal, sino los tiene entonces inicia como ciclo y despues se integre
    if begins: return('normal', ends, begins)
    else: return('cycleBegin', ends, begins)

def getRouteLinear(current_p, i, k, breakps=[]):
    path.append(current_p[i+1])
    # guardar el indice que se acaba de usar para despues no tomarlo en cuenta
    used_index = i+1
    k+=1
    # guardar los indices que corresponden a la siguiente reaccion. Estos tienen que ser pares porque sino serian sucesoras y no predecesoras
    prev_len = len(match_indexs)
    match_indexs.extend(index for index, item in enumerate(current_p) if path[k] == item and index != used_index and index%2==0 and not index in breakps)
    dif = len(match_indexs) - prev_len
    return(k, dif)


def getRouteCycleB(current_p, i, breakps=[]):
    path.insert(0,current_p[i-1])
    # guardar el indice que se acaba de usar para despues no tomarlo en cuenta
    used_index = i-1
    # guardar los indices que corresponden a la siguiente reaccion. Estos tienen que ser pares porque sino serian sucesoras y no predecesoras
    prev_len = len(match_indexs)
    match_indexs.extend(index for index, item in enumerate(current_p) if path[0] == item and index != used_index and index%2==1 and not current_p[index-1] in path and not index in breakps)
    dif = len(match_indexs) - prev_len
    return(dif)

list_all_pathways = getListPathways()

all_paths_routes = []
for current_p in list_all_pathways:
    routes_path = []
    type_p, ends, begins = getTypePath(current_p)
    if type_p == 'cycle':
        # guardar los primeros dos reacciones/items que esten en current_p. Al menos se puede estar seguro que estan unidas
        path = [current_p[0], current_p[1]]
        k = 1
        # obtener los indices de donde vuelve a aparcer la ultima de reaccion de 'path' en 'current_p'
        match_indexs = [index for index, item in enumerate(current_p) if path[k] == item and index != 1]
        # recorrer los indices obtenidos
        for i in match_indexs:
            # ver que nuestro indice no sea el ultimo
            if i+1 < len(current_p):
                path.append(current_p[i+1])
                # guardar el indice que se acaba de usar para despues no tomarlo en cuenta
                used_index = i+1
                k+=1
                # guardar los indices que corresponden a la siguiente reaccion. Estos tienen que ser pares porque sino serian sucesoras y no predecesoras
                match_indexs.extend(index for index, item in enumerate(current_p) if path[k] == item and index != used_index and index%2==0)
        routes_path.append(path)

    elif type_p == 'normal':
        for begin in begins:
            breakps = []
            while begin != None:
                # guardar los primeros dos reacciones/items que esten en current_p. Al menos se puede estar seguro que estan unidas
                path = [current_p[begin], current_p[begin+1]]
                k = 1
                # obtener los indices de donde vuelve a aparcer la ultima de reaccion de 'path' en 'current_p'
                match_indexs = [index for index, item in enumerate(current_p) if path[k] == item and index != begin+1 and not index in breakps and index%2 == 0]
                dif = 0
                # recorrer los indices obtenidos
                for i in match_indexs:
                    if dif <= 1:
                        # ver que nuestro indice no sea el ultimo
                        if i+1 < len(current_p):
                            # ver que no se trate de un ciclo, para no representarlo
                            if not current_p[i+1] in path:
                                k, dif = getRouteLinear(current_p, i, k)
                    else:
                        breakps = [i, i+1]
                        k, dif = getRouteLinear(current_p, i, k, breakps)

                routes_path.append(path)
                if breakps == []: break

    elif type_p == 'cycleBegin':
        for end in ends:
            breaker = 0
            while end != None:
                # guardar los primeros dos reacciones/items que esten en current_p. Al menos se puede estar seguro que estan unidas
                path = [current_p[end-1], current_p[end]]
                k = 0
                breakps = []
                # obtener los indices de donde vuelve a aparcer la ultima de reaccion de 'path' en 'current_p'
                match_indexs = [index for index, item in enumerate(current_p) if path[k] == item and index != end-1 and not index in breakps]
                dif = 0
                # recorrer los indices obtenidos
                for i in match_indexs:
                    if dif <= 1:
                        # ver que no se trate de un ciclo, para no representarlo
                        if not current_p[i-1] in path:
                            dif = getRouteCycleB(current_p, i)
                            if dif == 0: breaker = 1 ; break
                    else:
                        breakps = [i, i-1]
                        dif = getRouteCycleB(current_p, i, breakps)
                
                routes_path.append(path)
                if breaker: break

    all_paths_routes.append(routes_path)

for i in all_paths_routes:
    for j in i:
        print(j)

