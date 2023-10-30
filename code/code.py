# importar librerias
import pandas as pd 
from io import StringIO
import re

# funcion para leer arhivos de texto tabulares con texto de encabezado
def readTxts(location, begin):
    with open(location , 'r') as f:
        # leer por lineas para poder operarlo como lista de lineas
        # begin tiene el numero de linea en donde ya empieza la tabla
        output = f.readlines()[begin:]
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
lines_gene_IDs = readTxts('data/RegulonDB_geneidentifiers.txt', 16)
b_numbers = getPatterns(lines_gene_IDs, r'b\d{4}')
df_gene_IDs = makeDF(lines_gene_IDs)
df_gene_IDs = subsPatternToDF(df_gene_IDs, b_numbers, 5)

# hacer el dataframe para el archivo RegulonDB_NetworkTFGene.txt
lines_TF_gene = readTxts('data/RegulonDB_NetworkTFGene.txt', 37)
df_TF_gene = makeDF(lines_TF_gene, 0)

# leer el archivo de la via y guardar, en otra variable, sus bnumbers asociados
pathway = pd.read_csv('data/pathway_gns_test.txt', delimiter='\t', header=None)
bNumbers_pathway = pathway.iloc[:,1]

# funcion, toma items de una lista y los busca en una columna para entregar el valor de otra columna de la misma fila
def getInterest(df, list_search, c_to_search, c_to_subs=1):
    output_list = []
    for item in list_search:
        # regresa una lista con un True en la posicion en la que coincida la busqueda
        match = df.loc[df.iloc[:, c_to_search] == item]
        # sacar los indices de fila
        if not match.empty: # verificar que si exista un match, sino lo hay poner None en 'index_in_df'
            index_in_df = match.index[0]
        else: index_in_df = None
        # guardar en output el valor de la columna de interes
        # en base a si index_in_df tiene indice o vacio llenar segun sea el caso
        if index_in_df != None:
            output_list.append(df.iloc[index_in_df, c_to_subs])
        else: output_list.append('unknown')
    return(output_list)

# funcion que transforma de bNumber a TF
def getTFsPathway(initial_geneIDs):
    genes = getInterest(df_gene_IDs, initial_geneIDs, 5) # la columna 5 contiene los bNumbers
    TFs = getInterest(df_TF_gene, genes, 4) # la columna 4 contiene el nombre de los genes
    output_df = pd.DataFrame({'bNumber': initial_geneIDs, 'gene': genes, 'TF': TFs})
    return(output_df)

bN_gene_TF = getTFsPathway(bNumbers_pathway)
    
fractions = []
type_subpath = []
tf_most_ocurred = []

# funcion que calcula la fraccion de una subvia controlada por un mismo TF, el mas ocurrente
def getFraction(l, counts, case):
    most_ocurrence = counts[case]
    fractions.append(most_ocurrence/l)
    tf_most_ocurred.append(counts.index[case])

# funcion que prepara todo lo que necesita getFraction para despues llamarla
def callerFraction(current_path, len_sub):
    # .value_counts() saca, de una lista, cuantas veces ocurre cada item unico
    counts = current_path.value_counts()
    # sacar la longitud de la via 
    len_current_path = len(current_path)
    # preguntar de que caso se trata, por lo general no queremos que se tomen en cuenta los 'unkwown's.
    case = 0 if counts.idxmax() != 'unknown' or len(counts) == 1 else 1
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
    output = pd.DataFrame({'subpath': type_subpath, 'TF': tf_most_ocurred,  'fraction': fractions})   
    return(output)  
subpath_tf_fraction = getSubpathways()

# dar formato a las variables a imprimir para el output final
final_ouput_1 =subpath_tf_fraction.to_string(index=False) # pasar el df a string para quitar los indices de python (van del 0 a nfilas y no representa nada)
final_ouput_2_1 =subpath_tf_fraction.iloc[0,1] # tomar el FT mas representado de la via entera
final_ouput_2_2 = round(subpath_tf_fraction.iloc[0,2]*100, 2) # redondear a dos decimales y pasarlo a porcentaje la fraccion del FT mas represenatado

# imprimir el output del script
print(f'Lista de subvias con el factor de transcripcion que mas la regula y la fraccion de la (sub)via regulada por este mismo:\n{final_ouput_1}')
print(F'\nEl factor de transcripcion mas representado de la via es \'{final_ouput_2_1}\' con un {final_ouput_2_2}%') 