## Edición y aplicación de filtros en archivos VCF

Antes de comenzar, asegúrate de tener instalado Conda. Luego, sigue estos pasos para instalar los programas necesarios:


```bash
conda install -c bioconda bcftools
conda install -c bioconda vcftools
```

Esto te permitirá instalar los programas requeridos utilizando Conda.


# Explorar y editar información dentro de los archivos VCFs

Por lo general, la edición de archivos VCFs se realiza mediante scripts o programas especializados que permiten visualizar, agregar, modificar o eliminar elementos de manera controlada y documentada.

Conocer el nombre de las muestras presentes en el vcf
```
bcftools query -l archivo.vcf
```
Cambiar el nombre de las muestras
```
bcftools reheader --samples ID.pop.txt -o archivo2.vcf archivo1.vcf
```
 Anotar información relevante en el archivo (cromosoma y posición)
```
bcftools annotate --set-id +'%CHROM\_%POS\' -Ov -o archivo2.vcf archivo1.vcf
```
Conocer el número de sitios variante e invariantes en VCF
```
grep  -c -v "^#" input.vcf
```
Conocer sólo el número  de los sitios variantes
```
egrep -v "^#" input.vcf | wc -l
```

## Filtros Básicos para archivos VCFs

Los filtros básicos en un archivo VCF (Variant Call Format) identifican y seleccionan variantes genéticas que cumplen con ciertos criterios de calidad y confiabilidad en el proceso de llamada de variantes.Por ejemplo: 

valor mínimo de Q de 30 (en escala phred)
```
vcftools --gzvcf archivo1.vcf.gz --minQ 30  --recode --recode-INFO-all --out archivo2
```
Fitro con base en la profundidad media de cobertura de las variantes genéticas. Aquí se filtra con base en un umbral mínimo (--min-meanDP) o máximo (--max-meanDP) para la profundidad media de cobertura. En este ejemplo, se conservarán solo las variantes que tienen una profundidad media de cobertura igual o superior a 1.7 e igual o inferior a 50

```
vcftools --vcf archivo1.vcf --min-meanDP 1.7 --max-meanDP 10 --recode  --recode-INFO-all --out archivo2
```
Cantidad mínima de reads para considerar un genotipo. Se establecer el número mínimo de reads que deben cubrir una variante para que el genotipo sea considerado válido.
```
vcftools --vcf archivo1.vcf --minDP 3 --recode --recode-INFO-all --out archivo2
```
Máximo de missing information del 0.95 o  'Proporción de individuos en el estudio para los cuales la información correspondiente de SNP NO está ausente.
```
vcftools --vcf archivo1.vcf --max-missing 0.95 --recode  --recode-INFO-all --out archivo2
```
Enmascarar los indels y snps alrededor de indels. En este paso se oculta o filtrar las variantes de inserción/deleción (indels) y los SNPs que están cerca de las indels en un archivo VCF. Este filtro se hace en general  para evitar los errores o SNPs falsos alrededor de los indels.
```
bcftools filter --SnpGap 5 --IndelGap 5 -Ov -o archivo2.vcf archivo1.vcf
```
Eliminar todo los missing data. Muchos programas son sensibles a los "missing information", en este caso conviene eliminar todos los datos faltantes con el script a continuación.
```
 bcftools view -e 'GT[*] = "mis"' archivo1.vcf > archivo2.vcf
```

#  Filtros específicos 

Minor Allele Frequency (MAF). Filtra con base en la frecuencia con la que ocurre el alelo menos común en una población. Dependiendo de la pregunta de investigación, la naturaleza de los datos conviene ajustar el parámetro para filtrar bajo este criterio.
```
vcftools --vcf archivo1.vcf --maf 0.01 --recode --recode-INFO-all --out archivo2
```

Existen otros filtros como los que se mencionan a continuación:
> + Minor Allele Count (MAC)
> + Heterocigociddad
> + Hardy & Weinberg
> + Desequilibrio de Ligamiento (Lo veremos en el segundo día)


Eliminar cromosomas de un vcf. Muchas veces, queremos eliminar cromosomas de nuestro análisis, por ejemplo los cromosomas sexuales o mitocondriales, en este caso aplicamos este filtro:

```
vcftools --vcf archivo1.vcf --not-chr nombre_del_cromosoma --recode --recode-INFO-all --out archivo2
```

Crear un vcf sólo con ciertos cromosomas. En ocasiones, puedes explorar la estructura, selección o mezcla por cromosomas entre poblaciones, para ello aplicamos este filtro:

```
vcftools --vcf archivo1.vcf --chr nombre_del_cromosoma --recode --recode-INFO-all --out archivo2
```

