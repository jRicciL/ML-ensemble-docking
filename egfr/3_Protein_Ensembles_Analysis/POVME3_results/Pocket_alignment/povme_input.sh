
# Nombre de la trayectoria en formato pdb
PDBFileName                   ../../PDB_175_egfr_POCKET_ALL.pdb 

# Resolucion (tamaño de las sondas en A)
GridSpacing                    1.0 
# Datos de la esfera de inclusión (centor X Y Z y radio en A)
# The same center used for docking
InclusionSphere        -43	-6	2	14

# CRITERIOS DE DEFINICIÓN DE LA CAVIDAD
DistanceCutoff                          1.09
# Método de exclusión automático que determina donde acaba el pocket
ConvexHullExclusion                     first
# Se dfine una esfera menor, que indica cual es la cavidad principal
# en caso de que un frame muestre cavidades discontinuas, así POVME
# sólo mide la principal
SeedSphere      -43	-6	2	4.0	
# Número de puntos en común que determina si dos cavidades son continuas (una sola)
ContiguousPointsCriteria        3

# EJECUCIÓN
NumProcessors               2
OutputFilenamePrefix        EGFR-POCKET-RESIDUES_align/res_

