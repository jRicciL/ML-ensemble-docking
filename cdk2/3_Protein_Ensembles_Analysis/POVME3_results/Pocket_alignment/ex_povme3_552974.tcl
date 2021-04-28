# Nombre de la trayectoria en formato pdb
PDBFileName                     PDB_402_cdk2_PISANI_ALL.pdb

# Resolucion (tamaño de las sondas en A)
GridSpacing                     1.0
# Datos de la esfera de inclusión (centor X Y Z y radio en A)
InclusionSphere         -12.5     206.5     113.8     12.0

# CRITERIOS DE DEFINICIÓN DE LA CAVIDAD
DistanceCutoff                          1.09
# Método de exclusión automático que determina donde acaba el pocket
ConvexHullExclusion                     first
# Se dfine una esfera menor, que indica cual es la cavidad principal
# en caso de que un frame muestre cavidades discontinuas, así POVME
# sólo mide la principal
SeedSphere      -12.5     206.5     113.8     4.0
# Número de puntos en común que determina si dos cavidades son continuas (una sola)
ContiguousPointsCriteria        3

# EJECUCIÓN
NumProcessors               128
OutputFilenamePrefix        ./CDK2_VOL_PISANI_402_2/res_
