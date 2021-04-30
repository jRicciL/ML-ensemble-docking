
# Nombre de la trayectoria en formato pdb
PDBFileName                   PDB_136_fxa_SECONDARY-STRUCT-RESIDUES.pdb 

# Resolucion (tamaño de las sondas en A)
GridSpacing                    1.0 
# Datos de la esfera de inclusión (centor X Y Z y radio en A)
InclusionSphere        7	61	70	15 

# CRITERIOS DE DEFINICIÓN DE LA CAVIDAD
DistanceCutoff                          1.09
# Método de exclusión automático que determina donde acaba el pocket
ConvexHullExclusion                     first
# Se dfine una esfera menor, que indica cual es la cavidad principal
# en caso de que un frame muestre cavidades discontinuas, así POVME
# sólo mide la principal
SeedSphere      7	61	70	4.0	
# Número de puntos en común que determina si dos cavidades son continuas (una sola)
ContiguousPointsCriteria        3

# EJECUCIÓN
NumProcessors               4 
OutputFilenamePrefix        FXA_SEC-STRUCT-RESIDUES_align/res_

