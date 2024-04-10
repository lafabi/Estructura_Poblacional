# Qué es un archivio VCF?

Es un archivo de texto de llamado de variantes que tiene una extensión .vcf (archivo VCF) y que es  es el resultado de un proceso de bioinformático donde se congrega información de los sitios genéticos variantes y/o invariantes una o más muestras. Por lo general, una muestra de ADN es secuenciada a través de un sistema de secuenciación, produciendo archivos de secuencia cruda. En el proceso de resecuenciación, los datos de secuencia cruda experimentan una serie de ediciones, filtros, anotaciones, luego se alinean contra un genoma de referencia, creando archivos BAM/SAM como resultado. A partir de ahí se realiza el llamado de variantes que no es más que la  identificación de cambios en un genoma particular en comparación con el genoma de referencia. Esa salida se almacena en un formato de llamada de variantes, VCF que es un acrónimo de *Variant Calling File*

En el **Formato de Llamada de Variantes (VCF)**, hay 3 secciones principales en cada archivo:

1. **Líneas de Información Meta:** Múltiples líneas precedidas por símbolos de doble almohadilla (`##`).

2. **Línea de Encabezado:** Una sola línea precedida por un símbolo de almohadilla (`#`).

3. **Líneas de Datos:** El resto del archivo con 1 posición por línea.

##Ejemplo:


![Encabezado de VCFs](https://github.com/lafabi/Figuras/Fig1.vcf-head.png)

**Figura 1.** Encabezado de archivo VCF

Línea de cabecera

Cada archivo VCF tiene una única línea de encabezado que tiene 8 campos obligatorios separados por pestañas que representan columnas para cada línea de datos:

#CHROM POS ID REF ALT CALIDAD INFORMACIÓN FILTRO

Si hay datos de genotipo, se declara una columna FORMATO seguida de nombres de muestra únicos. Todos estos nombres de columnas también deben estar separados por pestañas.

Líneas de datos

Cada línea de datos representa una posición en el genoma. Los datos corresponden a las columnas especificadas en el encabezado y deben estar separados por tabulaciones y terminar con una nueva línea.

A continuación se muestran las columnas y sus valores esperados. En todos los casos, los valores FALTANTES deben representarse con un punto ('.').

    #CHROM - Identificador de cromosomas. Los ejemplos incluyen 7, chr7, X o chrX.

    POS - Posición de referencia. Ordenados numéricamente en orden ascendente por cromosoma.

    ID: identificadores únicos separados por punto y coma. No se permiten espacios en blanco.

    REF - Base de referencia (ACGT). Las inserciones se pueden representar con un punto.

    ALT: Bases alternativas separadas por comas (ACGT). Eliminaciones representadas por un punto.

    QUAL: puntuación de calidad en una escala logarítmica. 100 significa 1 entre 10^10 posibilidades de error.

    FILTER: indica qué filtros han fallado (separados por punto y coma), PASA o FALTA.

    INFO: información a nivel de sitio (no de muestra) en formato de nombre-valor separado por punto y coma.

    FORMATO: declaraciones de nombres de campos a nivel de muestra separadas por punto y coma.

    <SAMPLE DATA>: datos de campo a nivel de muestra separados por punto y coma correspondientes a declaraciones de campo FORMATO.


    Position and Ref/Alt Information

Below are some notes to help understand the first 5 columns about the above file.

    All of the variants occur on Chromosome 20 on the NCBI36 (hg18).

    There are 5 positions identified (14370, 18330, 1110696, 1230237, 1234567).

    Three of the variants have IDs including 2 dbSNP records (rs6054257, rs6040355).

    The first two positions (14370, 17330) are simple single-base pair substitutions.

    The third position has 2 alternate alleles specified (G and T) that replace the ref (A).

    The fourth position represents a deletion of a T since the alt allele is missing (“.”).

    The fifth row has 2 alt alleles, the first is a deletion of TC and second is insertion of a T.

QUAL and FILTER columns

The QUAL column indicates the quality level of the data at that site. The FILTER column designates what filters can be applied. The 2nd row (position 17330), has triggered the q10 filter, which is described in the meta section as “Quality below 10”.

Each bioinformatics pipeline treats these columns differently, so you will need to consult your pipeline’s subject matter experts on how to best interpret this information.

INFO column

The info column includes position-level information for that data row and can be thought as aggregate data that includes all of the sample-level information specified.

FORMAT column

The format column specifies the sample-level fields to expect under each sample. Each row has the same format fields (GT, GQ, DP, and HQ) except for the last row which does not have HQ. 

Each of these fields is described in the Meta section as the following:

    GT (Genotype) indicates which alleles separated by / (unphased) or | (phased).

    GQ is Genotype Quality which is a single integer.

    DP is Read Depth which is a single integer.

    HQ is Haplotype Quality and has 2 integers separated by a comma.
