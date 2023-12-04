import pandas as pd
from io import StringIO
from ast import literal_eval


# funcion para leer arhivos de texto tabulares con texto de encabezado
def readTxts(location):
    with open(location , 'r') as f:
        # leer por lineas para poder operarlo como lista de lineas
        output = f.readlines()
        # distingue las lineas del encabezado porque empiezan con '#' y las omite
        output = [l for l in output if not l.startswith('#')]
        return(output)

# funcion para pasar de lista de lineas a un string de todo el texto
def makeDF(lines, header_status=None):
    table = ''.join(lines)
    df = pd.read_csv(StringIO(table), delimiter='\t', header=header_status)
    return(df)

# crear el data frame y quedarse unicamente con la columna que contiene los FTs 
lines = readTxts('./with_globals.txt')
df = makeDF(lines, 0)
tfs = df.iloc[:, 2]

tfs_t = []

# funcion para ocnvertir representaciones de listas en listas
def convert_string_to_list(item):
    try:
        return literal_eval(item)
    except (ValueError, SyntaxError):
        # regresar el item origianl sino puede ser convertido
        return item

# ir haciendo una lista que contenga los FTs ya como items y no como listas
for item in tfs:
    i = convert_string_to_list(item)
    if type(i) is str:
        tfs_t.append(i)
    else:
        tfs_t.extend(i)

# escribir todo en un archivo
with open('tf_with_globals.txt', 'w') as file:
    file.write('TF\n')

with open('tf_with_globals.txt', 'a') as file:
    for item in tfs_t:
        file.write(item + '\n')

