
# Instalaciones
En principio se puede hacer todo con [NetCoMi](https://github.com/stefpeschel/NetCoMi). Está instalado en el módulo `R/4.3.1` del cluster. Copio y pego las instrucciones de instalación del paquete:

> ```r
> # Required packages
> install.packages("devtools")
> install.packages("BiocManager")
> 
> # Install NetCoMi
> devtools::install_github("stefpeschel/NetCoMi", 
>                          dependencies = c("Depends", "Imports", "LinkingTo"),
>                          repos = c("https://cloud.r-project.org/",
>                                    BiocManager::repositories()))
> ```
> If there are any errors during installation, please install the missing dependencies manually.
> 
> In particular the automatic installation of SPRING and SpiecEasi (only available on GitHub) does sometimes not work. These packages can be installed as follows (the order is important because SPRING depends on SpiecEasi):
> 
> ```r
> devtools::install_github("zdk123/SpiecEasi")
> devtools::install_github("GraceYoon/SPRING")
> ```


# Preprocesado

Algunas veces me he encontrado con redes a las que les salen unas estructuras un poco raras.

Esa estructura de abajo no es muy normal. Después de muchas pruebas, creo que esa especie de cogollo, piña :pineapple: o lo que sea viene de nodos con abundancias muy pequeñas, y suele desaparecer si se hace un **filtrado**. A mí incluso me han salido bichos con abundancia = 0 en todas las muestras, que no sé muy bien qué hacían en mi objeto `phyloseq` en primer lugar. 
 

# Generar las redes
Hay distintos métodos para generar las redes. Nosotras generalmente utilizamos SPIEC-EASI desde NetCoMi, un paquete de R que permite generar las redes y hacer distintos análisis con ellas.

La información del [GitHub de NetCoMi](https://github.com/stefpeschel/NetCoMi) es bastante completa, pero dejo por aquí algunas cosas.

> Ahora la documentación está en [netcomi.de](https://netcomi.de)

Se usa la función `netConstruct`. **No hace falta hacer la transformación CLR**, forma parte de los pasos que aplica SPIEC-EASI. El resultado es un objeto `microNet`.

```r
net.mho <- netConstruct(data = physeq.mho,
                        dataType = "counts",
                        measure = "spieceasi",
                        sparsMethod = "none",
                        nboot = 1000,
                        cores = 4L) # numero de cores
```

> SpiecEasi tiene dos modos a la hora de ejecutarlo, "glasso" (por defecto) y "mb". Muchos artículos especifican que se usó MB, aunque no tengo muy clara la explicación matemática. El caso es que para correrlo en ese modo se puede añadir el parámetro `measurePar = list(method='mb', icov.select.params = list(rep.num = 100))` a `netConstruct`. 

Después se puede hacer un objeto `microNetProps`, que saca algunas características de la red y permite además visualizarla.

```r
props.mho <- netAnalyze(net.mho,
                        centrLCC = TRUE,
                        clustMethod = "cluster_fast_greedy",
                        hubPar = c("degree", "betweenness"),
                        hubQuant = 0.9,
                        weightDeg = FALSE, normDeg = FALSE)
```

Los parámetros `hubPar` y `hubQuant` definen las propiedades y el cuantil utilizados para definir los **keystone taxa**. En este caso, NetCoMi identificará como keystone taxa aquellos nodos (bichos) con un grado y _betweenness_ mayor al 90% de los nodos. 


Los **keystone taxa** son microorganismos que en teoría son necesarios para que el ecosistema pueda mantener sus funciones. Las dos propiedades que se utilizan para identificarlos son características que ayudan a identificar nodos que, a nivel matemático, son relevantes en la red:
- **Grado:** es el número de vecinos que tiene un nodo, es decir, el número de conexiones que tiene.
- **Betweenness:** es una medida de centralidad. Aquí entra otro concepto, que es el camino más corto entre dos nodos de una red. Si calculamos todos los caminos más cortos entre todos los nodos de la red, y vemos cuántos pasan por cada nodo, tendríamos su betweenness. En redes de transporte es un concepto bastante intuitivo, Atocha sería un nodo con un betweenness muy alto en la red de Cercanías :nerd_face: 


# Visualizar las redes

La visualización se puede hacer desde R con NetCoMi, pero reconozco que yo no me apaño y me gusta más usar [Cytoscape](https://cytoscape.org/), que es un programa con interfaz gráfica que permite hacer muchas cosas y personalizar bastante las redes de forma muy sencilla. Para pasar de R a Cytoscape tengo esta función de R:

```r

exportNet <- function(net, physeq, props, filename){
  # EDGES
  # Create edge object from the edge list exported by netConstruct()
  edges <- net$edgelist1 %>%
    select(v1, v2, asso) %>%
    rename(Source = v1,
           Target = v2,
           Weight = asso) %>%
    mutate(Type = "Undirected",
           Sign = if_else(Weight < 0, "Negative", "Positive"))
  
  # NODES
  # get taxonomic info:
  taxa <- data.frame(physeq@tax_table) %>%
    rownames_to_column(var = "Label")
  
  # get mean abundances:
  mean_ab <- data.frame(
    "Abundance" = apply(data.frame(physeq@otu_table), 1, mean)) %>%
    rownames_to_column(var = "Label")
  
  # get clusters:
  clusters <- data.frame("Cluster" =  props$clustering$clust1) %>%
    rownames_to_column(var = "Label")
  
  # get hubs:
  hubs <- data.frame(isHub = rep(1, length(props$hubs$hubs1)),
                     "Label" = props$hubs$hubs1)
  
  # join all of them together:
  metadata <- taxa %>%
    left_join(mean_ab, by = "Label") %>%
    left_join(clusters, by = "Label") %>%
    left_join(hubs, by = "Label") %>%
    mutate(isHub = if_else(is.na(isHub), 0, 1))
  
  
  # WRITE CSV FILES:
  write_csv(edges,
            file = paste0(filename, "_edges.csv"),
            quote = "needed")
  write_csv(metadata,
            paste0(filename, "_metadata.csv"),
            quote = "needed")
}
```

Esta función usa como parámetros el objeto microNet (`net`), el objeto phyloseq que se utilizó para generarlo(`physeq`), el objeto microNetProps (`props`), y el nombre que queramos darle al archivo (`filename`). Por ejemplo:

```r
exportNet(net.mho, physeq.mho, props.mho, "mho")
```

La salida de esto son dos archivos que podemos importar a Cytoscape:
- `mho_edges.csv`: matriz con todas las ramas (edges) que unen unos nodos con otros, junto con el peso asociado a cada una. Se importa a Cytoscape con File > Import > Network from File.
- `mho_metadata.csv`: contiene las características de cada nodo: su abundancia relativa, si son un keystone taxa o no (columna `isHub`), y todos los niveles taxonómicos. Una vez importado el fichero anterior, se importa a Cytoscape con File > Import > Table from File.

Una vez importados los dos ficheros, se puede empezar a personalizar los colores, tamaños de los nodos, su disposición (_Layout_)... También tiene herramientas de análisis (Tools > Analyze Network)

> Estos ficheros vienen bien también si se quiere mover la red a otros paquetes, librerías, etc., como NetworkX


# Analizar las redes

Para sacar una tabla con algunas propiedades básicas de la red, en realidad nos sirve el objeto `microNetProps`:

```r
summary(props.mho)
```

```
Component sizes
    
size: 164
   #:   1
______________________________
Global network properties

Whole network:
                                 
Number of components      1.00000
Clustering coefficient    0.14338
Modularity                0.67388
Positive edge percentage 90.39146
Edge density              0.02102
Natural connectivity      0.00759
Vertex connectivity       1.00000
Edge connectivity         1.00000
Average dissimilarity*    0.99288
Average path length**     3.57500
-----
*: Dissimilarity = 1 - edge weight
**: Path length = Units with average dissimilarity
```

> También se podría hacer en Cytoscape, que tiene una opción "Analyze Network" en el menú Tools. Saldrá un cuadro preguntando si queremos analizarlo como un grafo dirigido, no hay que seleccionar esa opción porque nuestros edges no tienen direccionalidad.

![image](https://hackmd.io/_uploads/ryCigt-6C.png)

- **Component sizes / Number of components:** si la red tiene una sola componente conexa o varias. Si se genera una sola red, es una sola componente, si se forman varias redes chiquititas, se dice que hay varias componentes conexas. Si fuera el caso, además de sacar las propiedades de toda la red, saca las de la subred más grande (Largest connected component (LCC)), que son las que yo miraría.

Un par de propiedades importantes son: 

- **Clustering coefficient:** mide el nivel de agrupamiento de los nodos, es un valor entre 0-1. 
- **Average path length:** también se llama a veces shortest path length, mide cómo de separados están los nodos. Mide la distancia de un nodo a los demás nodos del grafo. En español se llama camino característico. En teoría, las redes con caminos característicos cortos son más robustas a ataques. 

Estas propiedades se calculan para cada nodo y se hace la media (es lo que se devuelve cuando da propiedades de toda la red). Son propiedades relevantes porque se usan para definir el "tipo" de estructura que tiene la red.

Del resto de propiedades, quizá destacaría la densidad de edges, que da una idea de la proporción de nodos y edges que hay, y la modularidad.

> Hace poco me encontré [este artículo](https://academic.oup.com/bib/article/22/2/1639/5734243), que no me he leído completo pero en el que me gustó la tabla de propiedades de redes, y explica bastante bien.



Más abajo da información sobre los _hubs_, que son los keystone taxa:

```
 k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Rikenellaceae|g__Alistipes|s__Alistipes_onderdonkii                       
 k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Eubacteriaceae|g__Eubacterium|s__Eubacterium_sp_CAG_274                       
 k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Blautia|s__Ruminococcus_torques                            
 k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Roseburia|s__Roseburia_sp_CAG_303                          
 k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Oscillospiraceae|g__Oscillibacter|s__Oscillibacter_sp_PC13                    
 k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Ruminococcus|s__Ruminococcus_callidus                      
 k__Bacteria|p__Firmicutes|c__Erysipelotrichia|o__Erysipelotrichales|f__Erysipelotrichaceae|g__Catenibacterium|s__Catenibacterium_mitsuokai
 k__Bacteria|p__Firmicutes|c__Negativicutes|o__Veillonellales|f__Veillonellaceae|g__Veillonella|s__Veillonella_infantium 
 ```
 
 Y da información sobre las medidas de centralidad de cada uno de ellos. 
 
 
# Comparar redes
NetCoMi permite también comparar las redes entre ellas. Esto computacionalmente es un poco costoso y a veces tarda... Yo estoy intentando optimizar la forma de hacerlo en el cluster, pero esta semana ya sabes que se está complicando. En cualquier caso el código sería algo así:

```r
comp_muo_mhno.net <- netConstruct(data = physeq.mhno,
                                  data2 = physeq.muo,
                                  dataType = "counts",
                                  measure = "spieceasi",
                                  sparsMethod = "none",
                                  nboot = 1000,
                                  cores = 4L)
saveRDS(comp_muo_mhno.net, "MHNOvsMUO_net_a4f.rds")

comp_muo_mhno.props <- netAnalyze(comp_muo_mhno.net,
                                  centrLCC = TRUE,
                                  clustMethod = "cluster_fast_greedy",
                                  hubPar = c("degree", "betweenness"),
                                  weightDeg = FALSE, normDeg = FALSE)
saveRDS(comp_muo_mhno.props, "MHNOvsMUO_props_a4f.rds")

# 
# featVec <- data.frame(tax_table(physeq.mhno))$Phylum
# names(featVec) <- rownames(tax_table(physeq.mhno))
# paletaClus <- MetBrewer::met.brewer("Renoir",n = 12)
# cluslcol <- as.character(paletaClus) 
# 
# plot(comp_muo_mhno.props, 
#      nodeTransp = 15,
#      colorVec=cluslcol,
#      nodeColor = "feature",
#      featVecCol = featVec,
#      nodeSize = "fix",
#      hubTransp = 10,
#      posCol = "darkgray", 
#      edgeTranspLow = 20,
#      cexLabels=1,
#      highlightHubs = T,
#      mar = c(2, 3, 3, 5),
#      sameLayout = TRUE, layoutGroup = "union",
#      groupNames = c("MHNO", "MUO"))

comp_muo_mhno <- netCompare(comp_muo_mhno.props,
                            permTest = TRUE,
                            nPerm = 100L,
                            storeCountsPerm = TRUE,
                            fileStoreCountsPerm = c("countsPermMHNOvsMUO",
                                                    "countsPermMUO"),
                            storeAssoPerm = TRUE,
                            fileStoreAssoPerm = "assoPerm",
                            seed = 123456,
                            cores = 4L)
```

Hay muchos ejemplos en el GitHub de NetCoMi, de todas formas, para que busques la opción que más te guste. Yo esto todavía no lo tengo claro.


# Otros paquetes

Si haciendo todo esto te quedas con ganas de más, hay paquetes que permiten hacer un montón de cosas con redes, `igraph` en R y `NetworkX` en Python. Yo personalmente prefiero NetworkX porque es donde aprendí a analizar redes y ya tengo algunos scripts :nerd_face: 
