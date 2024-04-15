## Ahora, interactuemos con los Archivos VCF

Antes de comenzar, asegúrate de tener instalado Conda. Luego, sigue estos pasos para instalar los programas necesarios:

1. **bcftools**:  
```bash
conda install -c bioconda bcftools
conda install -c bioconda vcftools
conda install -c bioconda admixture
```

Esto te permitirá instalar los programas requeridos utilizando Conda.


# Conocer el número de sitios variante e invariantes en VCF

grep  -c -v "^#" input.vcf

# Sólo de los sitios variantes

egrep -v "^#" input.vcf | wc -l

#Conocer el nombre de las muestras presentes en el vcf

bcftools query -l archivo.vcf

#cambio el nombre de las muestras
~/bcftools/bcftools reheader --samples ID.poo.txt -o archivo2.vcf archivo1.vcf

#Anotar información relevante en el archivo (cromosoma y posición)
~/bcftools/bcftools annotate --set-id +'%CHROM\_%POS\' -Ov -o archivo2.vcf archivo1.vcf

##Filtros Básicos

#valor mínimo de Q de 30 (en escala phred)
~/vcftools/src/cpp/vcftools --gzvcf archivo1.vcf.gz --minQ 30  --minDP 3 --min-meanDP 1.7 --max-meanDP 10 --max-missing 0.95 --recode --recode-INFO-all --out archivo2.vcf

#Máximo de missing information del 0.95 o  'Proporción de individuos en el estudio para los cuales la información correspondiente de SNP NO está ausente.

~/vcftools/src/cpp/vcftools --vcf archivo1.vcf --max-missing 0.95 --recode --recode-INFO-all --out archivo2.vcf

#--minDP indica la cantidad mínima de reads para considerar un genotipo
~/vcftools/src/cpp/vcftools --vcf archivo1.vcf --minDP 3 --recode --recode-INFO-all --out archivo2.vcf

#Enmascarar los indels y snps alrededor de indels
~/bcftools/bcftools filter --SnpGap 5 --IndelGap 5 -Ov -o archivo2.vcf archivo1.vcf


#Eliminar todo el missing data

 bcftools view -e 'GT[*] = "mis"' archivo1.vcf > archivo2.vcf


##Filtros específicos (MAF)

~/vcftools/src/cpp/vcftools --vcf archivo1.vcf --maf 0.01 --recode --recode-INFO-all --out archivo2.vcf

#Enmascarar los indels y snps alrededor de indels
~/bcftools/bcftools filter --SnpGap 5 --IndelGap 5 -Ov -o archivo2.vcf archivo1.vcf


Edición

# Cambiar el nombre de las muestras
~/bcftools/bcftools reheader --samples ID.pop.txt -o archivo1.vcf archivo2.vcf

# Anotar ID con cromosoma y posición
~/bcftools/bcftools annotate --set-id +'%CHROM\_%POS\' -Ov -o archivo2.vcf archivo2.vcf

