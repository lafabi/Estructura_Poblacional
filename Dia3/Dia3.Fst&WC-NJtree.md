# Indices de Diferencición *Fst*

El índice Fst (originalmente de Fisher) de Weir y Cockerham es un parámetro de genética de poblaciones utlizado para cuantificar la diferenciación genética entre subpoblaciones. Este índice nos permite evaluar la proporción de la variación genética total presente entre las subpoblaciones en comparación con la variación genética total en toda la población.

El método de Weir y Cockerham descrito en 1984, es particularmente importante en genómica pues proporciona una estimación robusta pues es capaz de corregir la varianza y covarianza de los alelos por el tamaño de la muestra en relación al tamaño de la población total.

Calculamos indices de Fst weir and Cockherman por sitio

```
vcftools --vcf archivo.vcf --weir-fst-pop pop1.txt --weir-fst-pop pop2.txt --out pop1-pop2
```

```
rm(list=ls())
graphics.off()
set.seed(1)

setwd("path/")

library(tidyverse)
library(ggplot2)

# Read .weir.fst file
fst_data <- read.table("path/pop1-pop2.weir.fst", header = TRUE)

# Plot using ggplot2
ggplot(fst_data, aes(x = CHROM, y = WEIR_AND_COCKERHAM_FST)) +
  geom_point(size = 1, color = "red") + 
  theme_classic() +   
  labs(title = "Distribución de Fst",
       x = "CHROM",
       y = "WEIR_AND_COCKERHAM_FST")

Histograma por Cromosomas

ggplot(fst_data, aes(x = WEIR_AND_COCKERHAM_FST)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  labs(title = "Distribución de Fst",
       x = "Fst",
       y = "Frecuencia")
```
### Con la finalidad de hacer los cálculos estadístico del promedio y la desviación estandar de los índices, debemos tomar en consideración la normalización de los datos, pues pueden verse afectado por el sesgo del tamaño de cada cromosoma o scaffolds. Para ello debemos conocer el nombre y tamaño de cada cromosoma (esto lo obtenemos a partir de un archivo .fai que es un índice que se obtiene con *samtools faidx* del genoma de referencia).

```

fst_data <- fread("pop1-pop2.weir.fst")
tamaño_cromosoma <- fread("tamaño_cromosoma_scaffolds.txt")
merged_data <- merge(fst_data , tamaño_cromosoma, by="CHR")
merged_data$normalized_FST <- merged_data$WEIR_AND_COCKERHAM_FST / merged_data$size
mean_normalized_fst <- mean(merged_data$normalized_FST)
std_dev_normalized_fst <- sd(merged_data$normalized_FST)
t_test <- t.test(merged_data$normalized_FST)

print(paste("X_FST_normalizado:", mean_normalized_fst))
print(paste("DS_FST_Norm:", std_dev_normalized_fst))
print(paste("p_value:", t_test$p.value))

```

## Árboles Neighbor Joining 
Los métodos basados en  Neighbor Joining (NJ) son árboles de distancias y que pretenden complementar los análisis exploratorios de clusterización. Para ellos se crea una matriz de distancia  a partir de una matriz de distancias o similitud entre especies o secuencias. Existen muchos métodos en linea para construir árboles de distancias, en esta ocasión utilizaremos vcf2pop, para ello clonaremos el repositorio donde vive el programa y lo invocamos desde nuestra máquina local: 



```
https://github.com/sansubs/vcf2pop

git clone https://github.com/sansubs/vcf2pop.git

```
Abrimos el recurso en línea html haciendo point and click en el ícono del programa que resposa dentro del directiorio.

Seguidamente, cargamos el archivo *vcf* y escogemos el método de obtención de la matriz de distancia y graficamos en la misma página o descargamos el archivo newick. En este punto puedes copiar y pegar el archivo newick y graficar un árbol más coqueto con la siguiente herramienta en línea 

```
https://itol.embl.de/ 
```
