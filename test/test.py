# importar librerias
import pandas as pd 
from io import StringIO
import re
from operator import *
from itertools import chain

# funcion para leer arhivos de texto tabulares con texto de encabezado
def readTxts(location):
    with open(location , 'r') as f:
        # leer por lineas para poder operarlo como lista de lineas
        output = f.readlines()
        # distingue las lineas del encabezado porque empiezan con '#' y las omite
        output = [l for l in output if not l.startswith('#')]
        return(output)
    

# funcion para sacar un patron (bnumbers) de cada linea   
def getPatterns(lines, pattern):
    patterns = []
    # buscar en cada linea el patron que concuerda con la ocurrencia de un bnumber y guardarlos en una lista
    for line in lines:
        match = re.search(pattern, line)
        # hay que verificar que exista match, sino lo hay match guarda None
        if match:
            extraction = match.group()
            patterns.append(extraction)
        # sino hubo match entonces hay que guardar vacio en esa posicion
        else: patterns.append('')
    return(patterns)


# funcion para pasar de lista de lineas a un string de todo el texto
def makeDF(lines, header_status=None):
    table = ''.join(lines)
    df = pd.read_csv(StringIO(table), delimiter='\t', header=header_status)
    return(df)


# funcion para sustituir una columna con muchos valores (gene synonyms) por uno de interes (bnumber) 
def subsPatternToDF(df, subs, c_receptor):
    df.iloc[:, c_receptor] = subs
    return(df)

lines_pathway_rxns = readTxts('test/test.txt')
df_pathway_rxns = makeDF(lines_pathway_rxns, 0)

# hacer el dataframe para el archivo reaction_to_genes.txt
lines_reaction_bn = readTxts('data/reaction_to_genes.txt')
df_reaction_bn = makeDF(lines_reaction_bn, 0)

# hacer el dataframe para el archivo RegulonDB_geneidentifiers.txt
lines_gene_IDs = readTxts('data/RegulonDB_geneidentifiers.txt')
b_numbers = getPatterns(lines_gene_IDs, r'b\d{4}')
df_gene_IDs = makeDF(lines_gene_IDs)
df_gene_IDs = subsPatternToDF(df_gene_IDs, b_numbers, 5)

# hacer el dataframe para el archivo RegulonDB_NetworkTFGene.txt
lines_TF_gene = readTxts('data/RegulonDB_NetworkTFGene.txt')
df_TF_gene = makeDF(lines_TF_gene, 0)

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
    
    return list_all_pathways, pathway_names

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

def getBranches(current_p, begins, ends):
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
    return branches

def getIndexes(current_p, path, k, used_index, avoid, indexes, direction, dif, begins, ends):
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

list_all_pathways, pathways_names = getListPathways()

def getRoutes(current_p):
    routes_path = []
    discriminator = 0
    if max(pd.Series(current_p).value_counts()) == 2 and min(pd.Series(current_p).value_counts()) == 2:
        discriminator = 1
        for rxn in list(set(current_p)):
            disc = [index%2 for index,item in enumerate(current_p) if item == rxn]
            if sum(disc)!=1 or len(disc)!=2: discriminator = 0
    if discriminator:
        path = [current_p[0], current_p[1]]
        k = 1
        indexes = [index for index, item in enumerate(current_p) if path[k] == item and index%2 == 0 and index != 1]
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
        branches = getBranches(current_p, begins, ends)
        # agregar a la via los pares de rxns donde hay ramificaciones duplicadas e invertidas
        current_p.extend(branches)

        for begin in begins:
            avoid = []
            entered = 0
            counter = 1
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
                            dif, direction = getIndexes(current_p, path, k, used_index, avoid, indexes, direction, dif, begins, ends)

                        elif dif > 1:
                            n_branches_route += 1
                            avoid.append(index)
                            #avoid.insert(0, index)
                            avoid_temp.append(index)

                            if direction == 'down': used_index = index+1 
                            else: used_index = index-1 
                            path.append(current_p[used_index])

                            for c in range(dif-1):
                                indexes.pop()
                            dif, direction = getIndexes(current_p, path, k, used_index, avoid, indexes, direction, dif, begins, ends)

                            if path[-1] in ends or path[-1] in begins:
                                routes_path.append(path)
                
                #if n_branches_route > 1:
                    #avoid = avoid[:-(n_branches_route-1)]

                if n_branches_route > 1:
                    #print(n_branches_route)
                    #print(f'avoid before before {avoid}')
                    if entered:
                        avoid.extend(avoid[:counter-1])
                    #print(f'avoid before {avoid}')
                    avoid = avoid[-counter:]
                    memory_avoid = avoid
                    entered = 1
                    counter += 1
                    #print(f'avoid after {memory_avoid}')

                if not (path[::-1] in routes_path or path in routes_path):
                    routes_path.append(path)
                if avoid_temp == []: break
        
        for end in ends:
                avoid = []
                entered = 0
                counter = 1
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
                                dif, direction = getIndexes(current_p, path, k, used_index, avoid, indexes, direction, dif, begins, ends)

                            elif dif > 1:
                                n_branches_route += 1
                                avoid.append(index)
                                #avoid.insert(0, index)
                                avoid_temp.append(index)

                                if direction == 'down': used_index = index+1 
                                else: used_index = index-1 
                                path.append(current_p[used_index])

                                for c in range(dif-1):
                                    indexes.pop()
                                dif, direction = getIndexes(current_p, path, k, used_index, avoid, indexes, direction, dif, begins, ends)
                                    
                    #if n_branches_route > 1:
                        #avoid = avoid[:-(n_branches_route-1)]

                    if n_branches_route > 1:
                        #print(n_branches_route)
                        #print(f'avoid before before {avoid}')
                        if entered:
                            avoid.extend(avoid[:counter-1])
                            #print(f'avoid before {avoid}')
                        avoid = avoid[-counter:]
                        memory_avoid = avoid
                        entered = 1
                        counter += 1
                        #print(f'avoid after {memory_avoid}')

                    if not (path[::-1] in routes_path or path in routes_path):
                        routes_path.append(path)
                    if avoid_temp == []: break

    return routes_path

def doShit(rout, pathway_name):
    def getInterestBNumbers(df, rout):
        bnumbers_rout = []
        for item in rout:
            match = []
            match = [index for index, item_col_search in enumerate(df.iloc[:, 0]) if item_col_search == item]
            if len(match)>0:
                bnumbers_per_rxn = []
                # el for sirve para recorrer los indices guardados en 'match' 
                for i in match:
                    # guardar en una lista temporal los FTs por cada gen
                    if not df.iloc[i, 1] in bnumbers_per_rxn: # saegurarse de no repetir un bnumber para cada rxn
                        bnumbers_per_rxn.append(df.iloc[i, 1])
                bnumbers_rout.append(bnumbers_per_rxn)
            else: bnumbers_rout.append('unknown')

        return bnumbers_rout

    # funcion que toma la lista de bnumbers de la via y extrae su correspondiente gen (bnumber >> gen)
    def getInterestGenes(df, bnumbers_rxn, c_to_search=5, c_to_subs=1):# la columna 5 contiene los bNumbers
        genes_rout = []
        for bnumbers in bnumbers_rxn:
            genes_per_rxn = []
            for bnumber in bnumbers:
                # obtener el indice donde existe el match, sino hubiera match entonces se guarda vacio
                match = df.loc[df.iloc[:, c_to_search] == bnumber] # en df.iloc se regresa una lista con un True en la posicion en la que coincida la busqueda
                # verificar que si exista un match, sino lo hay poner 'unknown' en 'output_list'
                if not match.empty: 
                    genes_per_rxn.append(df.iloc[match.index[0], c_to_subs])
                else: genes_per_rxn.append('unknown')
            genes_rout.append(genes_per_rxn)

        return genes_rout

    # funcion que toma la lista de genes de la via y extrae los correspondientes TFs (gen >> TFs)
    def getInterestTFs(df, genes_route, c_to_search=4, c_to_subs=1): # la columna 4 contiene el nombre de los genes
        tfs_rout = []
        for genes in genes_route:
            unknowns_counter = 0
            tfs_rxn = []
            for gen in genes:
                # como puede haber mas de un FT por cada gen entonces se crea una lista de match's
                match = []
                # se usa una list comprehension para sacar los indices donde exista un match entre el gen buscado y la columna con los genes
                match = [index for index, item_col_search in enumerate(df.iloc[:, c_to_search]) if item_col_search == gen]
                # verificar que si haya match's, sino los hubiera entonces match estaria vacio
                if len(match) > 0:
                    tfs_gene = []
                    # el for sirve para recorrer los indices guardados en 'match' 
                    for i in match:
                        # guardar en una lista temporal los FTs por cada gen
                        if not df.iloc[i, c_to_subs] in tfs_gene and not df.iloc[i, c_to_subs] in tfs_rxn: # saegurarse de no repetir un TF para cada gen
                            tfs_gene.append(df.iloc[i, c_to_subs])
                    tfs_rxn.extend(tfs_gene)
                else: unknowns_counter += 1

            if unknowns_counter == len(genes):
                # hay que guardarlo como lista, aunque solo sea un elemento, sino luego hay problemas al aplanar en callerFraction 
                tfs_rout.append(['unknown']) 
            else: tfs_rout.append(tfs_rxn)

        # al final se regresa una lista de listas, en la que las listas interiores contienen los FTs asociados a cada gen
        return tfs_rout

    bns = getInterestBNumbers(df_reaction_bn, rout)
    genes = getInterestGenes(df_gene_IDs, bns)
    TFs = getInterestTFs(df_TF_gene, genes)
    # hacer un DataFrame con los bnumbers, el nombre de los genes y los FTs que los regulan
    bN_gene_TF = pd.DataFrame({'bNumber': bns, 'gene': genes, 'TF': TFs})
        
    fractions = []
    type_subpath = []
    tf_most_ocurred = []
    occurrences = []
    unknowns = []
    other_TFs = []

    # funcion que calcula la fraccion de una subvia controlada por un mismo TF, el mas ocurrente
    def getFraction(l, counts, case, uk, current_path):
        most_ocurrence = counts[case]
        fractions.append(most_ocurrence/l)
        # guardar todos los factores de transcripcion mas representados de la via, para los casos en los que no solo es uno
        current_path = current_path.tolist()
        if counts.tolist().count(most_ocurrence) != 1:
            # evitar guardar el mismo factor dos veces
            temp = list(set([item for item in current_path if current_path.count(item) == most_ocurrence and item != 'unknown']))
            # si solo se encontro un FT hay que evitar guardarlo como lista
            if len(temp) == 1: tf_most_ocurred.append(temp[0])
            else: tf_most_ocurred.append(temp)
        else: 
            tf_most_ocurred.append(counts.index[case])
        occurrences.append(most_ocurrence)
        # contar cuantos otros factores de transcripcion regulan la subvia pero no fueron los mas recurrentes
        # si es 'unknown' el mas recurrente, entonces no hay otros factores
        if tf_most_ocurred[-1] == 'unknown': other_TFs.append(0)
        # si se apendeo un solo FT y no una lista de FTs entonces puede que si haya otros FTs
        elif tf_most_ocurred[-1] == str: other_TFs.append(l - most_ocurrence - uk)
        # si todos los factores de la subvia tambien son los mas recurrentes entonces no quedan mas TFs
        elif all(item in current_path for item in tf_most_ocurred[-1]): other_TFs.append(0)
        # si no se cumplio nada de lo anterior entonces tambien puede ser que queden algunos FTs
        else: other_TFs.append(l - most_ocurrence - uk)

    # funcion que prepara todo lo que necesita getFraction para despues llamarla
    def callerFraction(current_path, len_sub):
        # tenemos una lista de lista, hay que 'aplanarla'. El metodo chain hace esto, despues hay que pasarlo a tipo de dato pd.Series
        current_path_flat = pd.Series(list(chain(*current_path)))
        # .value_counts() saca, de una lista, cuantas veces ocurre cada item unico
        counts = current_path_flat.value_counts()
        # sacar la longitud de la via 
        len_current_path = len(current_path)
        # preguntar de que caso se trata, por lo general no queremos que se tomen en cuenta los 'unkwown's.
        # en este caso estamos viendo sobre 'counts' que esta aplanado, por eso 'unknown' no va entre parentesis
        case = 0 if counts.idxmax() != 'unknown' or len(counts) == 1 else 1
        # sacar el numero de 'unknown's por cada subvia
        # en este caso estamos viendo sobre la lista de listas 'current_path', hay que preguntar por ['unknown'] 
        uk = countOf(current_path, ['unknown']) if countOf(current_path,  ['unknown']) else 0
        getFraction(len_current_path, counts, case, uk, current_path_flat)
        unknowns.append(uk)
        type_subpath.append(len_sub)

    # funcion que hace la particion de la via en subvias y pregunta si sus FTs son el mismo
    def getSubpathways():
        # guardar en una variable unicamente los TFs, ya en orden, de la via
        TFs = bN_gene_TF['TF']
        len_path = len(TFs)
        subroutes = []
        pathways = []
        # en este primer for se van sacando los tamanos de subvia que puede haber. En este se recorren las sublongitudes
        for len_sub in range(len_path, 1, -1):
            when_stop = len_path-len_sub+1
            final_pos = len_sub
            # en este segundo for se van haciendo las subvias con el tamano acorde al primer for. En este se recorren las 'combinaciones' de cada sublongitud
            for initial_pos in range(0, when_stop):
                subpath = TFs[initial_pos:final_pos]
                callerFraction(subpath, len_sub)
                final_pos+=1
                subroutes.append(rout[initial_pos:final_pos])
                pathways.append(pathway_name)
        # crear el dataframe de salida
        output = pd.DataFrame({'pathway': pathways, 'length': type_subpath, 'TF': tf_most_ocurred,  'fraction': fractions, 'occurrences': occurrences, 'unknowns': unknowns, 'other TFs': other_TFs, 'subpath': subroutes})   
        
        return output
    
    final_df = getSubpathways()

    # dar formato a las variables a imprimir para el output final
    final_ouput_1 = final_df.to_csv(index=False, header=0, sep='\t').replace('\n', '') # pasar el df a string para quitar los indices de python (van del 0 a nfilas y no representa nada)

    # guardar el ouput en el archivo de salida
    with open('output.txt', 'a') as out_file:
        out_file.write(f'{final_ouput_1}')

with open('output.txt', 'w') as out_file:
    out_file.write('# Output file from one.py\n')
    out_file.write('# pathway\tlength\tTF\tfraction\tocurrences\tunknowns\tother_TFs\tsubroute\n')
  
for current_p, pathway_name in zip(list_all_pathways, pathways_names):
    routes_path = getRoutes(current_p)
    print(len(routes_path))
    for rout in routes_path:
        print(f'{rout}\n')
        doShit(rout, pathway_name)