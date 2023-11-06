import test as t
import pandas as pd  

lines_pathway_rxns = t.readTxts('data/pathway_rxns_test.txt')
df_pathway_rxns = t.makeDF(lines_pathway_rxns, 0)
print(df_pathway_rxns)

row = 0

list_all_pathways = []
pathway_name = df_pathway_rxns.iloc[0,0]
p_memory = []

# el for sirve para recorrer las filas de 'df_pathways_rxns', va otorgando los indices
for row in range(df_pathway_rxns.shape[0]):
    # este if sirve para conocer cuando existe un cambio de pathway, cuando ya no es la misma con la que se trabajo en la iteracion pasada
    if df_pathway_rxns.iloc[row,0] != pathway_name:
        list_all_pathways.append(p_memory) # guardar la via ya terminada en 'pathways'
        p_memory = [] # resetear la lista donde se guarda la via actual/con la que se trabaja
        pathway_name = df_pathway_rxns.iloc[row,0] # actualizar el nombre de la via
    # esto se hace para cada iteracion, guardar las reacciones de la fila en cuestion
    p_memory.append(df_pathway_rxns.iloc[row, 1])
    p_memory.append(df_pathway_rxns.iloc[row, 2])
    # saber si es la ultima iteracion de todo el for, la ultima fila. No puede ir junto con el primer if porque sino no se guardan las dos ultimas reactiosn de esta ultima via
    if row == df_pathway_rxns.shape[0]-1: list_all_pathways.append(p_memory)


pathways = []
'''
for current_p in list_all_pathways:
    path = [current_p[0], current_p[1]]
    # aqui es mejor iterar por indices y no por los valores/items de la 'current_p' porque asi se puede hacer sumas y comparaciones lejanas
    k = 1   
    match = [index for index, item in enumerate(current_p) if path[k] == item and index != 1]
    for i in match:
        print(match, i)
        k += 1  
        for j in range(i+1, len(current_p), 2):
            # if  not current_p[j] in path and current_p[j+1] == path[0]: path.insert(0, current_p[j])
            if not current_p[j] in path[:-1]:
                # este if es para ver si el final de la 'current_p' no regresa al inicio de la via
                if not current_p[j+1] in path: 
                    path.append(current_p[j+1])
                    used_index = (j+1)
                    appended = 1
        if appended:         
            match = [index for index, item in enumerate(current_p) if path[k] == item and index != used_index]
    print(path)
    ''' 

# for que recorra cada una de las vias
for current_p in list_all_pathways:
    # guardar los primeros dos reacciones/items que esten en current_p. Al menos se puede estar seguro que estan unidas
    path = [current_p[0], current_p[1]]
    k = 1
    # obtener los indices de donde vuelve a aparcer la ultima de reaccion de 'path' en 'current_p'
    match_indexs = [index for index, item in enumerate(current_p) if path[k] == item and index != 1]
    # recorrer los indices obtenidos
    for i in match_indexs:
        # ver que nuestro indice no sea el ultimo
        if i+1 < len(current_p):
            # ver que no se trate de un ciclo, para no representarlo
            if not current_p[i+1] in path:
                path.append(current_p[i+1])
                # guardar el indice que se acaba de usar para despues no tomarlo en cuenta
                used_index = i+1
                k+=1
                # guardar los indices que corresponden a la siguiente reaccion. Estos tienen que ser pares porque sino serian sucesoras y no predecesoras
                match_indexs.extend(index for index, item in enumerate(current_p) if path[k] == item and index != used_index and index%2==0)
    
    # guardar una lista con las reacciones que no se usaron de current_p
    left_rxns = [rxn for rxn in set(current_p) if not rxn in set(path)]
    #print(left_rxns)
    #print(path)
    # ir sacando las reacciones que no se usaron
    for left in left_rxns:
        path_var = path
        # indice de la reaccion 'left' en 'current_p'
        index_left = current_p.index(left)
        # indice de la reaccion/item siguiente de 'left' en 'current_p'
        index_next = current_p.index(left) + 1
        if current_p.count(left) == 0:
            # hay que preguntar si left se integra a la via, ya sea al incio o en medio, o por su parte es el inicio de una integracion a la via
            if current_p[index_next] in path and index_left%2 == 0:
                # hay que ver si la reaccion que le seguia a left es la que teniamos al inicio en 'path'
                if current_p[index_next] == path[0]: 
                    path.insert(0, left)
                    print('entre inicio')
                    left_rxns.remove(left)
                # sino fuera el caso, entonces left se mete a la via, por el medio
                else:
                    index_cut = path_var.index(current_p[index_next])
                    path_var = path_var[index_cut:]
                    path_var.insert(0, left)
                    print('entre mitad')
                    left_rxns.remove(left)
        else:
            path_var = [left]
            k = 0   
            match_indexs = [index for index, item in enumerate(current_p) if path_var[k] == item and index != 1]
            for i in match_indexs:
                if i+1 < len(current_p):
                    if not current_p[i+1] in path_var:
                        path_var.append(current_p[i+1])
                        used_index = i+1
                        k+=1
                        match_indexs.extend(index for index, item in enumerate(current_p) if path_var[k] == item and index != used_index and index%2==0)
            for rxn in path_var:
                if rxn in left_rxns: left_rxns.remove(rxn)
        print(left_rxns)
        print(f'path var: {path_var}')