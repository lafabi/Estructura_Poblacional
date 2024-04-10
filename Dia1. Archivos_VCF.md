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


Información de posición y referencia/alt.

A continuación hay algunas notas que le ayudarán a comprender las primeras 5 columnas sobre el archivo anterior.

    Todas las variantes ocurren en el cromosoma 20 del NCBI36 (hg18).

    Hay 5 posiciones identificadas (14370, 18330, 1110696, 1230237, 1234567).

    Tres de las variantes tienen ID que incluyen 2 registros dbSNP (rs6054257, rs6040355).

    Las dos primeras posiciones (14370, 17330) son sustituciones simples de un solo par de bases.

    La tercera posición tiene 2 alelos alternativos especificados (G y T) que reemplazan a la referencia (A).

    La cuarta posición representa una eliminación de una T ya que falta el alelo alt (“.”).

    La quinta fila tiene 2 alelos alt, la primera es una eliminación de TC y la segunda es la inserción de una T.

Columnas CALIDAD y FILTRO

La columna CUAL indica el nivel de calidad de los datos en ese sitio. La columna FILTRO designa qué filtros se pueden aplicar. La segunda fila (posición 17330) ha activado el filtro q10, que se describe en la metasección como "Calidad inferior a 10".

Cada canal bioinformático trata estas columnas de manera diferente, por lo que deberá consultar a los expertos en la materia de su canal sobre cómo interpretar mejor esta información.

columna INFORMACIÓN

La columna de información incluye información a nivel de posición para esa fila de datos y puede considerarse como datos agregados que incluyen toda la información a nivel de muestra especificada.

columna FORMATO

La columna de formato especifica los campos a nivel de muestra que se esperan en cada muestra. Cada fila tiene los mismos campos de formato (GT, GQ, DP y HQ) excepto la última fila que no tiene HQ.

Cada uno de estos campos se describe en la sección Meta de la siguiente manera:

    GT (Genotipo) indica qué alelos separados por / (sin fase) o | (por fases).

    GQ es Calidad del Genotipo, que es un número entero único.

    DP es profundidad de lectura, que es un número entero único.

    HQ es calidad de haplotipo y tiene 2 números enteros separados por una coma.

Ícono de validado por la comunidad
