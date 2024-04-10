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

    CUAL: puntuación de calidad en una escala logarítmica. 100 significa 1 entre 10^10 posibilidades de error.

    FILTRO: indica qué filtros han fallado (separados por punto y coma), PASA o FALTA.

    INFORMACIÓN: información a nivel de sitio (no de muestra) en formato de nombre-valor separado por punto y coma.

    FORMATO: declaraciones de nombres de campos a nivel de muestra separadas por punto y coma.

    <DATOS DE MUESTRA>: datos de campo a nivel de muestra separados por punto y coma correspondientes a declaraciones de campo FORMATO.
