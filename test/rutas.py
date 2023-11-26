import test as t
import pandas as pd  

lines_pathway_rxns = t.readTxts('test/test.txt')
df_pathway_rxns = t.makeDF(lines_pathway_rxns, 0)
print(df_pathway_rxns)

def getListPathways():
    row = 0
    pathway_names = [df_pathway_rxns.iloc[0,0]]
    list_all_pathways = []
    p_memory = []
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

def getBeginsEnds(current_p):
    begins, ends = [],[]
    # sacar las rxns que sean inicios y finales de la via
    # son aquellas que aparecen una vez o que siempre aparecen en posiciones iniciales(pares) o siempre en posiciones finales (impares)
    counts = pd.Series(current_p).value_counts()
    # obtener las rxns que aparecen una vez
    uniques = counts[counts==1].index
    # discernir si son inicios o finales
    for item_s in uniques:
        begins.extend(index for index, item in enumerate(current_p) if item_s == item and index%2 == 0)
        ends.extend(index for index, item in enumerate(current_p) if item_s == item and index%2 == 1)
    
    # obtener rxns que aparecen mas de una vez
    not_uniques = [item for item in current_p if current_p.count(item)>1]
    # ir recorriendo las rxns que aparecen mas de una vez
    for item_s in list(set(not_uniques)):
        c_begin, c_end = 0, 0
        # obtener los indices de donde aparece la rxn en cuestion
        indexes = [index for index, item in enumerate(current_p) if item_s == item]
        # el siguiente for sirve para asegurar que todos los indices que corresponden a una rxn sean solo inicios o solo finales
        # en caso de que asi lo fuera entonces se agregan dichos indices como posiciones de inicios
        for index in indexes:
            if index%2 == 1 and c_begin == 0: c_end += 1
            elif index%2 ==0 and c_end == 0: c_begin += 1
            else: break
            if c_end == len(indexes): ends.extend(indexes)
            elif c_begin == len(indexes): begins.extend(indexes)
    return(begins, ends)

def getBranches(current_p):
    branches = []
    branch_rxns = list(set([item for index, item in enumerate(current_p) if current_p.count(item)>2 and not index in begins and not index in ends]))
    branch_rxns.extend(list(set([item for index, item in enumerate(current_p) if current_p.count(item)>1 and (index in begins or index in ends)])))
    for item_s in branch_rxns:
        indexes = [index for index, item in enumerate(current_p) if item_s == item]
        c_pre, c_post = 0, 0
        pre, post = [], []
        for index in indexes:
            if index%2 == 0:
                c_pre += 1
                pre.append(current_p[index+1])
                pre.append(current_p[index])
            else:
                c_post += 1
                post.append(current_p[index])
                post.append(current_p[index-1])
        if c_pre > c_post:
            branches.extend(pre)
        else:
            branches.extend(post)
    return(branches, branch_rxns)

def getIndexes(current_p, path, k, used_index, avoid, indexes, direction, dif):
    def getIndexesDown(current_p, path, k, used_index, avoid, indexes):
        prev_len = len(indexes)
        for index, item  in enumerate(current_p):
            if item == path[k] and index != used_index and index%2 == 0 and not index in avoid and not current_p[index+1] in path:
                indexes.append(index)
        dif = len(indexes) - prev_len
        return(dif)

    def getIndexesUp(current_p, path, k, used_index, avoid, indexes):
        prev_len = len(indexes)
        for index, item  in enumerate(current_p):
            if item == path[k] and index != used_index and index%2 == 1 and not index in avoid and not current_p[index-1] in path:
                indexes.append(index)
        dif = len(indexes) - prev_len
        return(dif)

    if direction == 'down':
        dif = getIndexesDown(current_p, path, k, used_index, avoid, indexes)
        if dif == 0 and not path[-1] in [current_p[i] for i in ends] and not path[-1] in [current_p[j] for j in begins]:
            direction = 'up'
            dif = getIndexesUp(current_p, path, k, used_index, avoid, indexes)
    else:
        dif = getIndexesUp(current_p, path, k, used_index, avoid, indexes)
        if dif == 0 and not path[-1] in [current_p[i] for i in ends] and not path[-1] in [current_p[j] for j in begins]:
            direction = 'down'
            dif = getIndexesDown(current_p, path, k, used_index, avoid, indexes) 
    return(dif, direction)

list_all_pathways = getListPathways()
all_paths_all_routes = []
routes_path = []

for current_p in list_all_pathways:
    discriminator = 0
    if max(pd.Series(current_p).value_counts()) == 2 and min(pd.Series(current_p).value_counts()) == 2:
        discriminator = 1
        for rxn in list(set(current_p)):
            disc = [index%2 for index,item in enumerate(current_p) if item == rxn]
            if sum(disc)!=1 or len(disc)!=2: discriminator = 0
    if discriminator:
        path = [current_p[0], current_p[1]]
        k = 1
        indexes = [index for index, item in enumerate(current_p) if path[k] == item and index%2 == 0 and index != begin+1 and not index in avoid and not current_p[index+1] in path]
        for i in indexes:
            # ver que nuestro indice no sea el ultimo
            if i+1 < len(current_p):
                path.append(current_p[i+1])
                # guardar el indice que se acaba de usar para despues no tomarlo en cuenta
                used_index = i+1
                k+=1
                # guardar los indices que corresponden a la siguiente reaccion. Estos tienen que ser pares porque sino serian sucesoras y no predecesoras
                indexes.extend(index for index, item in enumerate(current_p) if path[k] == item and index != used_index and index%2==0 and item != path[0])
        routes_path.append(path)                 
    else:
        # obtener inicios y finales
        begins, ends = getBeginsEnds(current_p)
        # obtener las rxns donde hay ramificaciones
        branches, branch_rxns = getBranches(current_p)
        # agregar a la via los pares de rxns donde hay ramificaciones duplicadas e invertidas
        current_p.extend(branches)

        for begin in begins:
            avoid = []
            while begin != None:
                direction = 'down'
                avoid_temp = []
                n_branches_route = 0
                # guardar los primeros dos reacciones/items que esten en current_p. Al menos se puede estar seguro que estan unidas
                path = [current_p[begin], current_p[begin+1]]
                k = 1
                indexes = [index for index, item in enumerate(current_p) if path[k] == item and index%2 == 0 and index != begin+1 and not index in avoid and not current_p[index+1] in path]
                dif = len(indexes)
                if indexes:
                    for index in indexes:
                        k += 1
                        if dif == 1:
                            if direction == 'down': used_index = index+1
                            else: used_index = index-1
                            path.append(current_p[used_index])
                            dif, direction = getIndexes(current_p, path, k, used_index, avoid, indexes, direction, dif)

                        elif dif > 1:
                            n_branches_route += 1
                            avoid.insert(0, index)
                            avoid_temp.append(index)

                            if direction == 'down': used_index = index+1 
                            else: used_index = index-1 
                            path.append(current_p[used_index])

                            for c in range(dif-1):
                                indexes.pop()
                            dif, direction = getIndexes(current_p, path, k, used_index, avoid, indexes, direction, dif)

                            if path[-1] in ends or path[-1] in begins:
                                routes_path.append(path)
                
                if n_branches_route > 1:
                    avoid = avoid[:-(n_branches_route-1)]

                if not (path[::-1] in routes_path or path in routes_path):
                    routes_path.append(path)
                if avoid_temp == []: break
        
        for end in ends:
                avoid = []
                while end != None:
                    direction = 'up'
                    avoid_temp = []
                    n_branches_route = 0
                    # guardar los primeros dos reacciones/items que esten en current_p. Al menos se puede estar seguro que estan unidas
                    path = [current_p[end], current_p[end-1]]
                    k = 1
                    indexes = [index for index, item in enumerate(current_p) if path[k] == item and index%2 == 1 and index != end-1 and not index in avoid and not current_p[index-1] in path]
                    dif = len(indexes)
                    if indexes:
                        for index in indexes:
                            k += 1
                            if dif == 1:
                                if direction == 'down': used_index = index+1
                                else: used_index = index-1
                                path.append(current_p[used_index])
                                dif, direction = getIndexes(current_p, path, k, used_index, avoid, indexes, direction, dif)

                            elif dif > 1:
                                n_branches_route += 1
                                avoid.insert(0, index)
                                avoid_temp.append(index)

                                if direction == 'down': used_index = index+1 
                                else: used_index = index-1 
                                path.append(current_p[used_index])

                                for c in range(dif-1):
                                    indexes.pop()
                                dif, direction = getIndexes(current_p, path, k, used_index, avoid, indexes, direction, dif)
                                    
                    if n_branches_route > 1:
                        avoid = avoid[:-(n_branches_route-1)]

                    if not (path[::-1] in routes_path or path in routes_path):
                        routes_path.append(path)
                    if avoid_temp == []: break

all_paths_all_routes.append(routes_path)

'''def definitiveDecision(query, contrast):
    for rxn_q, rxn_c in zip(query, contrast[::-1]):
        if rxn_c != rxn_q: return 0
        if rxn_c == query[-1]: 
            return 1

for i, path_routes in enumerate(all_paths_all_routes):
    deletions = []
    delete = 0
    print(len(path_routes))
    for j, rout_query in enumerate(path_routes):
        len_query = len(rout_query)
        for k, rout_contrast in enumerate(path_routes):
            if j == k or k in deletions: continue
            if len_query == len(rout_contrast) and rout_query[0] == rout_contrast[-1] and rout_query[-1] == rout_contrast[0]:
                delete = definitiveDecision(rout_query, rout_contrast)
            elif rout_query == rout_contrast: 
                delete = 1
            if delete:
                path_routes.pop(k)
                delete = 0
        deletions.append(j)
    print(len(path_routes))
 
'''
for i in all_paths_all_routes:
    print(len(i))
    for j in i:
        print(j)
        print('\n')