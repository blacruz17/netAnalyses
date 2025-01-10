des microbioma

Hay distintos métodos para generar las redes. Nosotras generalmente utilizamos SPIEC-EASI desde NetCoMi, un paquete de R que permite generar las redes y hacer distintos análisis con ellas.

La información del [GitHub de NetCoMi](https://github.com/stefpeschel/NetCoMi) es bastante completa, pero dejo por aquí algunas cosas.


# Generar las redes

Se usa la función `netConstruct`. **No hace falta hacer la transformación CLR**, forma parte de los pasos que aplica SPIEC-EASI. El resultado es un objeto `microNet`.

```r
net.mho <- netConstruct(data = physeq.mho,
                        dataType = "counts",
                        measure = "spieceasi",
                        sparsMethod = "none",
                        nboot = 1000,
                        cores = 4L)
```

Después se puede hacer un objeto `microNetProps`, que saca algunas características de la red y permite además visualizarla.

```r
props.mho <- netAnalyze(net.mho,
                        centrLCC = TRUE,
                        clustMethod = "cluster_fast_greedy",
                        hubPar = c("degree", "betweenness"),
                        hubQuant = 0.9,
                        weightDeg = FALSE, normDeg = FALSE)
```

Los parámetros `hubPar` y `hubQuant` definen las propiedades y el cuantil utilizados para definir los **keystone taxa**. En este caso, NetCoMi identificará como keystone taxa aquellos nodos (bichos) con un grado y _betweenness_ mayor