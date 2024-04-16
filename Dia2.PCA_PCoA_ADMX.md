## Análisis de Estructura PCA y PCoA
Los PCAs se obtienen como forma de análisis exploratorio de los datos, su principal función es la de reducir la dimensionalidad de un conjunto de datos. En genómica los PCAs permiten observar la  direccionalidad y patrones de agrupamiento, ordenándolos en función de cuánta varianza explican.

Instalemos lo programas que utilizaremos en esta sección

```bash
conda install -c bioconda plink
conda install -c bioconda admixture
```

# El primer paso es convertir el archivo VCF a formato Plink. 
```
plink --vcf archivo1.vcf --recode --out archivo2.vcf --allow-extra-chr
```
> + En este script la función --allow-extra-chr permite el uso de cromosomas adicionales en comparación a el número de cromosomas de humanos
> + 

# Pruning o purga primer paso, se identifican los SNPs que muestran desequilibrio de ligamiento. Cualquier SNP que presente una correlación con r igual o mayor a 0.1 será eliminado.
```
plink --file archivo2 --indep-pairwise 50 5 0.1 --allow-extra-chr
```
> + En este script la función --indep-pairwise 50 5 0.1 especifica la forma en que se realizará el pruning, se filtran los SNPs que tienen una correlación (r) igual o mayor a 0.1 pq se consideran que están  ligados. Los números 50 y 5 indican la ventana de tamaño y el paso entre SNPs considerados.

# Pruning o purga segundo paso, se excluyen los snps que no presentan desequilibrio de ligamiento.

```
plink --file archivo2 --extract plink.prune.in --recode vcf --out archivo3.vcf --allow-extra-chr
```
> + La función --extract plink.prune.in especifica el archivo que contiene los SNPs pruneados.
> + Mientras que la función --recode vcf indica que los archivos de salida se recodifican a formato VCF.

# Se realiza el PCA
```
plink --vcf archivo3.vcf --pca --out archivo3-LD-0.1 --allow-extra-chr
```


Luego de correr el script. obtendremos los archivos eigenvalue y eigenvector.Estos archivos son ampliamente utilizados en análisis multivariados. 

    Los eigenvalues son valores que caracterizan la escala o la magnitud de la varianza explicada por cada componente principal en un PCA. Representan la varianza explicada por cada componente principal y se muestran en orden descendente, lo que significa que el primer eigenvalue es el más grande y explica la mayor parte de la varianza en los datos.

    Los eigenvectors son vectores que representan las direcciones o ejes de máxima variabilidad en los datos. Cada eigenvector está asociado con un eigenvalue y define una componente principal en el espacio de características de los datos.
    Rrepresentan las combinaciones lineales de las variables originales que definen las nuevas variables ortogonales (componentes principales) que explican la mayor parte de la varianza en los datos, éstosse utilizan para transformar los datos originales en un nuevo conjunto de datos proyectados en el espacio de las componentes principales.


Una vez obtenido los eigenvalues y los eigenvectors procedemos a graficarlos en R

```
# Limpiamos el ambiente de R
rm(list=ls())
graphics.off()
set.seed(1)

# Cargamos las librerías y dependencia
library(ggplot2)
library(stringr)

### Observemos los Eigenvalues

# Cargamos los datos de eigenvalues
Val <- read.table("/home/fabiola/Documentos/Spheniscus/External-memory/Spheniscus/vcf-raw/Workshop/PCA3.eigenval") 
Val$PC <- c(1:20)  # concatenar valores del 1 a 17 como valores de la variable PC
colnames(Val) <- c("percent","PC") # nombrar columnas
print(Val)

# Graficamos  los eigenvalues en Barplot o Histogramas 
ggplot(Val[1:5,], aes(x=Val[1:5,]$PC, y=Val[1:5,]$percent), xlab = "PC") +
  geom_bar(stat = "identity", width = 0.5) +
  labs(y= "% variance", x = "PC")+
  lims(y=c(0,10)) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10))


### Observemos los Eigenvectors

## Cargamos los datos de eigenvector
Vec<-read.table("/home/fabiola/Documentos/Spheniscus/External-memory/Spheniscus/vcf-raw/Workshop/PCA3.eigenvec")  
Vec <- Vec[,(1:5), drop=FALSE] 
colnames(Vec)<-c("ID","Species","PCA1","PCA2","PCA3") 
print(Vec$Species)  
print(Vec$ID) 


## Graficamos los Componentes 1 vs 2.  PC1 vs PC2

p1 <- ggplot (Vec, aes(x= PCA1, y= PCA2, color = Species))+  
  geom_point(size=4)+  
  theme_classic()+  
  labs(x="PC1 (4.44%)", y="PC2 (2.84%)") +  
  theme(axis.title.x = element_text(face="bold", vjust=0, size=rel(1.5))) +  
  theme(axis.title.y = element_text(face="bold", vjust=1.5, size=rel(1.5)))+  
  theme(panel.border = element_rect(colour="black", fill=NA, size=1))  
print(p1)

p2 <- p1 + scale_colour_manual(values = c("#CA3F3F","#64AD3F","#E8B547","red"))
print(p2)


## Graficamos los Componentes 2 vs 3. PC2 vs PC3

p3 <- ggplot (Vec, aes(x= PCA2, y= PCA3, color = Species))+
  geom_point(size=4)+
  theme_gray()+
  labs(x="PC2 (1.99%)", y="PC3 (1.8%)") +
  theme(axis.title.x = element_text(face="bold", vjust=0, size=rel(1.5))) +
  theme(axis.title.y = element_text(face="bold", vjust=1.5, size=rel(1.5))) +
  theme(panel.border = element_rect(colour="black", fill=NA, size=1))
print(p3)
p4 <- p3 + scale_colour_manual(values = c("#CA3F3F","#64AD3F","#E8B547"))
print(p4)

```

## Análisis de Estructura hallar número de clústers más probables y el grado de mezcla


# Convertir a formato plink admisible reconocible para admixture (recode 12)
```
plink --vcf archivo1.vcf --recode 12 --out archivo2-ADMX --allow-extra-chr
```


# Pruning, primer paso. Se estima cuales están bajo desequilibrio de ligamiento. Aquí cualquier SNP que se correlacione con un r de 0.1 o mayor será eliminado
```
plink --file archivo2-ADMX --indep-pairwise 50 5 0.1 --allow-extra-chr
```

# Pruning, segundo paso. Se extraen los snps que no presentan desequilibrio de ligamiento.
```
plink --file archivo2-ADMX --extract plink.prune.in --recode --out archivo3-ADMX --allow-extra-chr 0
```

# Com este loop probaremos  distintos valores de K (cantidad de poblaciones ancestrales), con K={1 2 3 4 5 6 7 8 9 10}.
#El flag '--cv' permite evaluar varios K con una 'validación cruzada'. Calculan el error estandar de la validacion cruzada para cada K, y el de menor
#valor indica que tiene mayor sensibilidad.
for K in 1 2 3 4 5 6 7 8 9 10
do
admixture --cv archivo3-ADMX.ped $K -j8 | tee log${K}.out
done

# Para observar los distintos K ocupamos el comando *grep*
grep -h CV *.out > ID_log_admix.txt

Para graficar los resultados ocuparemos una herramienta en línea llamada pophelper:

http://pophelper.com/



    
