# Qué es un archivio VCF?

Es un archivo de texto de llamado de variantes que tiene una extensión .vcf (archivo VCF) y que es  es el resultado de un proceso de bioinformático donde se congrega información de los sitios genéticos variantes y/o invariantes una o más muestras. Por lo general, una muestra de ADN es secuenciada a través de un sistema de secuenciación, produciendo archivos de secuencia cruda. En el proceso de resecuenciación, los datos de secuencia cruda experimentan una serie de ediciones, filtros, anotaciones, luego se alinean contra un genoma de referencia, creando archivos BAM/SAM como resultado. A partir de ahí se realiza el llamado de variantes que no es más que la  identificación de cambios en un genoma particular en comparación con el genoma de referencia. Esa salida se almacena en un formato de llamada de variantes, VCF que es un acrónimo de *Variant Calling File*

En el **Formato de Llamada de Variantes (VCF)**, hay 3 secciones principales en cada archivo:

1. **Líneas de Información Meta:** Múltiples líneas precedidas por símbolos de doble almohadilla (`##`).

2. **Línea de Encabezado:** Una sola línea precedida por un símbolo de almohadilla (`#`).

3. **Líneas de Datos:** El resto del archivo con 1 posición por línea.

##Ejemplo:

##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##bcftoolsVersion=1.10.2-90-gf27f849+htslib-1.10.2-58-g2e36fa6
##bcftoolsCommand=mpileup -Ou -f /data6/testacc/fabiola/GCA_010078495.1_BGI_Enov.V1_genomic.fna -b /data6/testacc/fabiola/Fabi.bamlist -d 1000 -Q 20 -a DP,SP -r VULB01000001.1 –-skip-indels
##reference=file:///data6/testacc/fabiola/GCA_010078495.1_BGI_Enov.V1_genomic.fna
##contig=<ID=VULB01000001.1,length=1067>
##contig=<ID=VULB01000002.1,length=2639>
##contig=<ID=VULB01000003.1,length=4521>
##contig=<ID=VULB01000004.1,length=3054>
##contig=<ID=VULB01000005.1,length=3064>
