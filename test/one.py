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


# funcion que toma la lista de bnumbers de la via y extrae su correspondiente gen (bnumber >> gen)
def getInterestGenes(df, list_search, c_to_search=5, c_to_subs=1):# la columna 5 contiene los bNumbers
    output_list = []
    for item in list_search:
        # obtener el indice donde existe el match, sino hubiera match entonces se guarda vacio
        match = df.loc[df.iloc[:, c_to_search] == item] # en df.iloc se regresa una lista con un True en la posicion en la que coincida la busqueda
        # verificar que si exista un match, sino lo hay poner 'unknown' en 'output_list'
        if not match.empty: 
            output_list.append(df.iloc[match.index[0], c_to_subs])
        else: output_list.append('unknown')
    return(output_list)

# funcion que toma la lista de genes de la via y extrae los correspondientes TFs (gen >> TFs)
def getInterestTFs(df, list_search, c_to_search=4, c_to_subs=1): # la columna 4 contiene el nombre de los genes
    output_list = []
    for item in list_search:
        # como puede haber mas de un FT por cada gen entonces se crea una lista de match's
        match = []
        # se usa una list comprehension para sacar los indices donde exista un match entre el gen buscado y la columna con los genes
        match = [index for index, item_col_to_search in enumerate(df.iloc[:, c_to_search]) if item_col_to_search == item]
        # verificar que si haya match's, sino los hubiera entonces match estaria vacio
        if len(match) > 0:
            temp_list = []
            # el for sirve para recorrer los indices guardados en 'match' 
            for i in match:
                # guardar en una lista temporal los FTs por cada gen
                if not df.iloc[i, c_to_subs] in temp_list: # saegurarse de no repetir un TF para cada gen
                    temp_list.append(df.iloc[i, c_to_subs])
            output_list.append(temp_list)
        else: output_list.append(['unknown'])
    # al final se regresa una lista de listas, en la que las listas interiores contienen los FTs asociados a cada gen
    return(output_list)
    
fractions = []
type_subpath = []
tf_most_ocurred = []
occurrences = []
unknowns = []
other_TFs = []

# funcion que calcula la fraccion de una subvia controlada por un mismo TF, el mas ocurrente
def getFraction(l, counts, case, uk):
    most_ocurrence = counts[case]
    fractions.append(most_ocurrence/l)
    tf_most_ocurred.append(counts.index[case])
    occurrences.append(most_ocurrence)
    # contar cuantos otros factores de transcripcion regulan la subvia pero no fueron los mas recurrentes
    if tf_most_ocurred[-1] == 'unknown': other_TFs.append(0)
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
    # en este caso estamos viendo sobre la lista de listas 'curren_path', por eso hay que preguntar sobre ['unknown']
    uk = countOf(current_path, ['unknown']) if countOf(current_path, ['unknown']) else 0
    unknowns.append(uk)
    getFraction(len_current_path, counts, case, uk)
    type_subpath.append(len_sub)

# funcion que hace la particion de la via en subvias y pregunta si sus FTs son el mismo
def getSubpathways():
    # guardar en una variable unicamente los TFs, ya en orden, de la via
    TFs = bN_gene_TF['TF']
    len_path = len(TFs)
    # en este primer for se van sacando los tamanos de subvia que puede haber. En este se recorren las sublongitudes
    for len_sub in range(len_path, 1, -1):
        when_stop = len_path-len_sub+1
        final_pos = len_sub
        # en este segundo for se van haciendo las subvias con el tamano acorde al primer for. En este se recorren las 'combinaciones' de cada sublongitud
        for initial_pos in range(0, when_stop):
            subpath = TFs[initial_pos:final_pos]
            callerFraction(subpath, len_sub)
            final_pos+=1
    # crear el dataframe de salida
    output = pd.DataFrame({'subpath': type_subpath, 'TF': tf_most_ocurred,  'fraction': fractions, 'occurrences': occurrences, 'unknowns': unknowns, 'other TFs': other_TFs})   
    return(output)

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

def getRouteLinear(current_p, i, k, breakps=[], path):
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

def getPathroutes():
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
    return(all_paths_routes)


# hacer el dataframe para el archivo RegulonDB_geneidentifiers.txt
lines_gene_IDs = readTxts('data/RegulonDB_geneidentifiers.txt')
b_numbers = getPatterns(lines_gene_IDs, r'b\d{4}')
df_gene_IDs = makeDF(lines_gene_IDs)
df_gene_IDs = subsPatternToDF(df_gene_IDs, b_numbers, 5)

# hacer el dataframe para el archivo RegulonDB_NetworkTFGene.txt
lines_TF_gene = readTxts('data/RegulonDB_NetworkTFGene.txt')
df_TF_gene = makeDF(lines_TF_gene, 0)

# leer el archivo de la via y guardar, en otra variable, sus bnumbers asociados
pathway = pd.read_csv('data/pathway_gns_test.txt', delimiter='\t', header=None)
bNumbers_pathway = pathway.iloc[:,1]

genes = getInterestGenes(df_gene_IDs, bNumbers_pathway)
TFs = getInterestTFs(df_TF_gene, genes)
# hacer un DataFrame con los bnumbers, el nombre de los genes y los FTs que los regulan
bN_gene_TF = pd.DataFrame({'bNumber': bNumbers_pathway, 'gene': genes, 'TF': TFs})
print(bN_gene_TF)
  
final_df = getSubpathways()

# dar formato a las variables a imprimir para el output final
final_ouput_1 =final_df.to_string(index=False) # pasar el df a string para quitar los indices de python (van del 0 a nfilas y no representa nada)
final_ouput_2_1 =final_df.iloc[0,1] # tomar el FT mas representado de la via entera
final_ouput_2_2 = round(final_df.iloc[0,2]*100, 2) # redondear a dos decimales y pasarlo a porcentaje la fraccion del FT mas represenatado

# imprimir el output del script
print(f'Lista de subvias con el factor de transcripcion que mas la regula y la fraccion de la (sub)via regulada por este mismo:\n{final_ouput_1}')
print(F'\nEl factor de transcripcion mas representado de la via es \'{final_ouput_2_1}\' con un {final_ouput_2_2}%')

lines_pathway_rxns = readTxts('test/test.txt')
df_pathway_rxns = makeDF(lines_pathway_rxns, 0)
print(df_pathway_rxns)

list_all_pathways = getListPathways()

for i in all_paths_routes:
    for j in i:
        print(j)
