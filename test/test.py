# importar librerias
import pandas as pd 
from io import StringIO
import re
from operator import *

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

# hacer el dataframe para el archivo RegulonDB_geneidentifiers.txt
lines_gene_IDs = readTxts('data/RegulonDB_geneidentifiers.txt')
b_numbers = getPatterns(lines_gene_IDs, r'b\d{4}')
df_gene_IDs = makeDF(lines_gene_IDs)
df_gene_IDs = subsPatternToDF(df_gene_IDs, b_numbers, 5)

# hacer el dataframe para el archivo RegulonDB_NetworkTFGene.txt
lines_TF_gene = readTxts('data/RegulonDB_NetworkTFGene.txt')
df_TF_gene = makeDF(lines_TF_gene, 0)

# leer el archivo de la via y guardar, en otra variable, sus bnumbers asociados
pathway = pd.read_csv('test/test.txt', delimiter='\t', header=None)
bNumbers_pathway = pathway.iloc[:,1]

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
                temp_list.append(df.iloc[i, c_to_subs])
            output_list.append(temp_list)
        else: output_list.append('unknown')
    # al final se regresa una lista de listas, en la que las listas interiores contienen los FTs asociados a cada gen
    return(output_list)

genes = getInterestGenes(df_gene_IDs, bNumbers_pathway)
TFs = getInterestTFs(df_TF_gene, genes)
# hacer un DataFrame con los bnumbers, el nombre de los genes y los FTs que los regulan
bN_gene_TF = pd.DataFrame({'bNumber': bNumbers_pathway, 'gene': genes, 'TF': TFs})
print(bN_gene_TF)
    
fractions = []
type_subpath = []
tf_most_ocurred = []
occurrences = []
unknowns = []

# funcion que calcula la fraccion de una subvia controlada por un mismo TF, el mas ocurrente
def getFraction(l, counts, case):
    most_ocurrence = counts[case]
    fractions.append(most_ocurrence/l)
    tf_most_ocurred.append(counts.index[case])
    occurrences.append(most_ocurrence)

# funcion que prepara todo lo que necesita getFraction para despues llamarla
def callerFraction(current_path, len_sub):
    # .value_counts() saca, de una lista, cuantas veces ocurre cada item unico
    counts = current_path.value_counts()
    # sacar la longitud de la via 
    len_current_path = len(current_path)
    # preguntar de que caso se trata, por lo general no queremos que se tomen en cuenta los 'unkwown's.
    case = 0 if counts.idxmax() != 'unknown' or len(counts) == 1 else 1
    # sacar el numero de 'unknown's por cada subvia
    uk = countOf(current_path, 'unknown') if countOf(current_path, 'unknown') > 0 else 0
    unknowns.append(uk)
    getFraction(len_current_path, counts, case)
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
    output = pd.DataFrame({'subpath': type_subpath, 'TF': tf_most_ocurred,  'fraction': fractions, 'occurrences': occurrences, 'unknowns': unknowns})   
    return(output)
  
final_df = getSubpathways()

# dar formato a las variables a imprimir para el output final
final_ouput_1 =final_df.to_string(index=False) # pasar el df a string para quitar los indices de python (van del 0 a nfilas y no representa nada)
final_ouput_2_1 =final_df.iloc[0,1] # tomar el FT mas representado de la via entera
final_ouput_2_2 = round(final_df.iloc[0,2]*100, 2) # redondear a dos decimales y pasarlo a porcentaje la fraccion del FT mas represenatado

# imprimir el output del script
print(f'Lista de subvias con el factor de transcripcion que mas la regula y la fraccion de la (sub)via regulada por este mismo:\n{final_ouput_1}')
print(F'\nEl factor de transcripcion mas representado de la via es \'{final_ouput_2_1}\' con un {final_ouput_2_2}%') 