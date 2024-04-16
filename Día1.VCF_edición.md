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
vcftools --gzvcf archivo1.vcf.gz --minQ 30  --max-missing 0.95 --recode --recode-INFO-all --out archivo2
```
Fitro con base en la profundidad media de cobertura de las variantes genéticas
```
vcftools --vcf archivo1.vcf --min-meanDP 1.7 --max-meanDP 10 --recode  --recode-INFO-all --out archivo2
```
Cantidad mínima de reads para considerar un genotipo
```
vcftools --vcf archivo1.vcf --minDP 3 --recode --recode-INFO-all --out archivo2
```
Máximo de missing information del 0.95 o  'Proporción de individuos en el estudio para los cuales la información correspondiente de SNP NO está ausente.
```
vcftools --vcf archivo1.vcf --max-missing 0.95 --recode  --recode-INFO-all --out archivo2
```
Enmascarar los indels y snps alrededor de indels
```
bcftools filter --SnpGap 5 --IndelGap 5 -Ov -o archivo2.vcf archivo1.vcf
```
Eliminar todo los missing data
```
 bcftools view -e 'GT[*] = "mis"' archivo1.vcf > archivo2.vcf
```

#  Filtros específicos 

Minor Allele Frequency (MAF)
```
vcftools --vcf archivo1.vcf --maf 0.01 --recode --recode-INFO-all --out archivo2.vcf
```

Existen otros filtros como los que se mencionan a continuación:
> + Minor Allele Count (MAC)
> + Heterocigociddad
> + Hardy & Weinberg
> + Desequilibrio de Ligamiento (Lo veremos en el segundo día)



