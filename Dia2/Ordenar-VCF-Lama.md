# Cómo cambiamos el orden de las muestras en un VCF según su ID sample?

Tomaremos como ejemplo el archivo VCF de *Lama*. Este VCF posee un orden de ID que nos dificulta hacer los análisis posteriores de PCA y Admixture, este tutorial explicará cómo solucionarlo y cómo hacer si tenemos que ordenar otros archivos VCFs propios en el futuro. Lo primero que haremos es conocer el orden de las muestras en el VCF de *Lama*, para ello ocupamos el programa bcftools:

```
bcftools query -l guanaco-3chr.vcf
```

Idealmente ordenaremos el vcf tal cual aparece en el archivo lama.pop. Para ello,extraemos con grep la primera columna del archivo lama.pop en un nuevo archivo de texto que se llamará sort-lama-ID.txt. 

```
awk '{print $1} lama.pop > sort-lama-ID.txt

```
Ahora, ordenamos el VCF de guanaco tal como aparece en el archivo de texto sort-lama-ID.txt con el programa bcftools y generaremos un archivo VCF con un nuevo orden de esta forma: 

```
bcftools view -S sort-lama-ID.txt guanaco-3chr.vcf -o guanaco-3chr-sort.vcf

```
### Consideraciones generales

Para ordenar un archivo VCF deben coincidir tanto el número como los ID de las muestras en ambos archivos, el VCF y el .txt, pudiendo variar el orden pero conservando el mismo número y el ID de las muestras.

Para facilitar los análisis al trabajar con archivos VCFs conviene ordenar los ID de las muestras según nos acomode:  orden latitudinal ascendente o descendente, regiones simpátricas, puntos cardinales etc. 
