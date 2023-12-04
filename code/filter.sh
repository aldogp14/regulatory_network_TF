#!/bin/bash

# Establecer los TF que no queremos tener en cuenta
tf='CRP\|H-NS\|ArcA\|Lrp\|FNR\|IHF\|Fis'

# Filtrar las líneas basándose en la columna de TF
awk -v tf="$tf" '$2 !~ tf' ./data/RegulonDB_NetworkTFGene.txt > ./data/RegulonDB_NotGlobal_TFs.txt

# Contar las líneas del archivo filtrado
wc -l ./data/RegulonDB_NotGlobal_TFs.txt
