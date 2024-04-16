# Terminal bash: comandos básicos

**¿Qué es la terminal Bash o Shell?**

La "Shell" es una interfaz de usuario que opera a través de líneas de comando. Estos comandos son, de por sí, programas que actúan como interlocutores o intérpretes entre sistemas complejos y nosotros. Entre sus principales atributos destacan su capacidad de multitarea y multiusuario. La característica multitarea, como su nombre lo indica, permite realizar más de una función o tarea simultáneamente en varias ventanas o pestañas. Por otra parte, el rasgo multiusuario permite realizar distintas tareas a más de una persona en la misma terminal simultánemente.

Como mencionamos con anterioridad, nos podemos comunicar con el sistema operativo a través de un interlocutor; los comandos. Ahora, comenzaremos a aprender un lenguaje que nos permitirá interactuar y dar instrucciones al sistema para la ejecución de una tarea. Usualmente, los comandos son acrónimos de palabras en inglés y se definen según la función que desempeñan.

Primero hagamos un saludo al sistema operativo, haciendo que se imprima en la pantalla a través del comando ``echo`` (es un eco de lo que quieres imprimir en pantalla).

```
echo Hello World
```


## Navegar entre sistemas de archivos

Muchas veces nos resultará útil saber en qué carpeta estamos (desde ahora llamaremos a las carpetas directorios). Para ello utilizaremos el comando ``pwd`` (print working directory).

```
 pwd
```
```/home/usuario```

Generalmente cuando ingresamos por primera vez a la terminal, el sistema nos redirige al home del usuario. Será conveniente cambiarnos de directorio y para ello usaremos el comando ``cd`` (change directory).

```
cd Documentos
```

Ahora verifiquemos que nos cambiamos de directorio con pwd.

```
pwd
```
```/home/usuario/Documentos```

Dado que nos cambiamos a Documentos, ahora creémos un subdirectorio llamado Genobiostoic dentro de Documentos con el comando ``mkdir`` (make directory).

```
mkdir Genobiostoic
```

>**Nota**: En ocasiones, puede suceder que creamos un directorio en un lugar equivocado, o que simplemente queramos deshacernos de un directorio en desuso. El comando para eliminar un directorio dependerá de si este se encuentra o no vacío. Si un directorio **se encuentra vacío**, entonces utilizaremos el comando ``rmdir``, acrónimo de “remove directory”, para eliminarlo. Ahora bien, si el directorio que queremos eliminar **contiene archivos u otros directorios**, debemos utilizar el comando ``rm``, acrónimo de "remove", junto con otra instrucción, ``-r``, que indica al comando rm que la remoción debe ser *recursiva* para todo el contenido del directorio. Así, para remover directorios que no se encuentren vacíos, utilizaremos ``rm -r``. **¡Cuidado! La remoción de directorios mediante rm -r es riesgosa, pues no tiene vuelta atrás. Asegúrate de utilizar este comando con precaución, procurando que los archivos que eliminarás sean los correctos.** 


Verifiquemos que se haya creado el directorio Genobiostoic a través del comando ``ls`` (list).

```
ls
```

Aquí observamos que se ha creado un nuevo directorio llamado Genobiostoic, cuyo directorio parental es Documentos. Entremos al directorio con ``cd`` y verifiquemos el sistema anidado y jerárquico de directorios con ``pwd``.

```
pwd
```
```/home/usuario/Documentos/Genobiostoic```

En la terminal, las rutas relativas y absolutas se utilizan para especificar la ubicación de un archivo o directorio en el sistema de archivos. Una ruta absoluta, como su nombre lo indica, es una ruta completa que comienza desde la raíz o home del sistema de archivos y especifica la ubicación exacta de un archivo o directorio. En Linux toda ruta absoluta comienza con un slash (/). Por el contrario, una ruta relativa especifica la ubicación de un archivo o directorio en relación con la ubicación actual. No comienza con un slash (/) y depende la ruta donde se encuentre el usuario.


## Descargar archivos

Ahora que sabemos cómo crear directorios, movernos entre ellos y listar su contenido, aprendamos a descargar archivos desde repositorios remotos a través de líneas de comando. Para ello usaremos ```wget```. En esta oportunidad descargaremos el libro de literatura inglesa clásica Dr. Jeckyll y Mr.Hyde en el directorio Genobiostoic a través de este comando:


```
wget https://www.gutenberg.org/cache/epub/43/pg43.txt
```
Verifiquemos que se haya descargado el libro con el comando ```ls```.


El comando ```mv``` (move) permite cambiar el nombre de un archivo o directorio. Cambiemos el nombre del libro desde pg43.txt a libro.txt.


```
mv pg43.txt libro.txt
```

En este punto hemos dado una instrucción, mover o cambiar el nombre de un archivo (pg43.txt) por otro nombre de archivo (libro.txt).
El comando ```mv``` también permite mudar archivos entre directorios. Si queremos mudar "libro.txt" desde Genobiostoic a Practica1 escribamos el siguiente comando:

```
mv libro.txt Terminal/Libro

```

Luego, verifiquemos que efectivamente hayamos mudado el archivo libro.txt al directorio indicado.

```
cd Terminal/Libro
ls
```
```libro.txt```


Si queremos imprimir **todo** el libro en pantalla ocupamos el comando ```cat```, este nos mostrará todo el libro sin la posibilidad de editarlo.

A continuación haremos una práctica más extensa con diversos comandos para interactuar con distintos tipos de archivos.


## Interactuar con archivos

Lo primero que haremos en esta sección será crear un archivo con un editor de texto. La terminal shell trae un editor de texto instalado llamado ``nano``. Basta con invocarlo tipeando nano en la terminal.

```
nano
```

>**Nota**: Aunque nano es uno de los programas más comúnmente utilizados, existen una serie de otros editores de texto que pueden ser descargados en tu máquina local, tales como *Atom*, *Sublime Text* y *Visual Studio Code*. Te invitamos a probar, aparte de esta inducción, los distintos editores y evaluar cuál te permite trabajar de manera más cómoda. 

A continuación se abre una ventana en la que podemos pegar, tipear y/o editar el texto que deseemos. Presionando ^x salimos del editor de texto nano y regresamos a la terminal. Ahora que sabemos invocar el editor de texto procedamos a crear un documento de texto plano **.txt** llamado *verde*. Para ello copiemos el poema de Federico García Lorca (Disponible al final de esta sección) y seguidamente invoquemos ```nano``` y creémos al mismo tiempo el documento con su nombre y extensión:

```
nano verde.txt
```

Ahora peguemos el poema de Federico Garcia Lorca dentro del documento en nano, guardemos con ^o y salgamos de nano con ^x. Verifiquemos que se ha creado el documento con el comando ``ls``, debería aparecer el nombre del documento "verde.txt" en el directorio.

```
ls
```

En este punto aprenderemos varios comandos que nos permitirán interactuar con estos archivos. Como vimos, el programa ```nano``` nos confiere el permiso de editar o modificar el documento. Ahora bien, muchas veces queremos ver el documento, es decir, imprimir el contenido en la pantalla sin que sea modificado. Para esto ocupamos el comando ``cat`` (concat):

```
cat verde.txt
```

A continuación el texto contenido en el documento se imprimirá en la pantalla. Muchas veces, cuando trabajamos con datos genómicos muy grandes es poco conveniente imprimir todo su contenido, de hecho la concatenación del documento puede durar varios minutos. Si nos pasa esto podemos detener la impresión en pantalla con ^C. Si deseamos limpiar la pantalla, basta con tipear clear en la terminal:

```
clear
```

Como mencionamos, en ocasiones queremos imprimir el encabezado o cierto número de líneas del documento, para ello utilizamos el comando ``head`` (head):

```
head verde.txt
```

Varios comandos tienen ciertas funciones, por ejemplo ``head`` acompañado de la función ``-n`` (number) permite indicar el número de lineas a imprimir en pantalla. Por lo tanto, el comando ``-n 40`` imprime en pantalla las primeras 40 lineas del documento:

```
head -n 40 verde.txt
```


También, podemos imprimir el final de un documento con el comando ``tail``.

```
tail verde.txt
```

Al igual que ``head``, ``tail`` acepta la función ``-n`` y podemos imprimir el número de lineas que deseemos, contando desde el final hacia el principio del documento.


Cuando se trabajan con datos Genómicos, muchas veces queremos buscar, ubicar o filtrar ciertos loci, cromosomas, individuos, etc, dentro de un archivo. Para ello ocupamos el comando ``grep`` (globally search for regular expression and print out). En este ejemplo buscaremos la palabra verde en el texto verde.txt, y para ello tipearemos:

```
grep verde verde.txt
```

Aquí, vale la pena resaltar que dimos una instrucción distinta, ordenamos buscar (``grep``) una palabra (verde) en un archivo (verde.txt). De esta forma diseñamos la sintaxis de un programa un poco más complejo que los anteriores.

Ahora vayamos más allá y busquemos (``grep``) un palabra (*Compadre*) y contemos cuántas veces se repite con ``wc`` (word count) y la función ``-l``, que permite contar el número de líneas.

```
grep "Compadre" verde.txt | wc -l
```

