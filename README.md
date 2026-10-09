# Simulador del dogma central de la biología molecular

Práctica 1 de **Bioinformática** (Escuela de Ingeniería Informática, ULPGC).

Programa de consola en Python que, a partir de una molécula de ADN, simula paso a paso el flujo de
la información genética:

1. **Replicación** (ADN → ADN)
2. **Transcripción** (ADN → ARN)
3. **Traducción** (ARN → proteína)

Al final, un resumen alinea codón a codón la hebra codificante, la hebra molde, el ARNm, los
anticodones y los aminoácidos, para seguir el flujo completo.

## Ejecución

Necesita Python 3.8 o posterior. El simulador de texto no usa librerías externas; solo la salida
en imagen (`--imagen`) necesita `matplotlib` (`pip install -r requirements.txt`).

```bash
python3 simulador.py                    # menú: gen de ejemplo, secuencia propia, FASTA o aleatoria
python3 simulador.py ATGGCCAAGTTTGGCTAA # una secuencia (hebra codificante, 5'→3')
python3 simulador.py -f gen.fasta       # desde un archivo FASTA
python3 simulador.py -a 60              # secuencia aleatoria de 60 pb
python3 simulador.py --sin-pausa        # todas las pantallas seguidas
python3 simulador.py --imagen salida/   # además guarda PNG: un dibujo por paso de la replicación
                                        # y otro con el flujo ADN → ARNm → proteína
```

En la terminal se avanza de pantalla con **Enter** (**q** para salir). Si la salida se redirige a un
archivo (`python3 simulador.py --sin-pausa > salida.txt`), se escribe sin colores.
En [`ejemplos/salida_ejemplo.txt`](ejemplos/salida_ejemplo.txt) está la salida completa del gen de
ejemplo.

## Qué muestra

### 1. Replicación

Una horquilla de replicación avanza de izquierda a derecha. En cada paso se ve qué parte de la
molécula ya está abierta y cuál sigue unida por puentes de hidrógeno (`╎` = 2 en A=T, `┆` = 3 en
G≡C):

```
posición           1         11        21        31        41        51
molde superior  5' AGGAGGTAAACCATGAAAGCAATTTTCGTACTGAAAGGTTGGTGGCGCACTTCCTGAAGC 3'
                   ╎┆┆╎┆┆╎╎╎╎┆┆╎╎┆╎╎╎┆┆╎╎╎╎~~~~~~~~~~~~┆┆╎╎┆┆╎┆┆┆┆┆╎┆╎╎┆┆╎┆╎╎┆┆
rezagada  ←     3' TCCTCCATTTGGTACTTTCGTuaa            ┆┆╎╎┆┆╎┆┆┆┆┆╎┆╎╎┆┆╎┆╎╎┆┆ 5'
                   ├─────────◀Ok1─────────┤            H┆╎T┆┆╎┆┆┆┆┆╎┆╎╎┆┆╎┆╎╎┆┆
líder     →     5' aggAGGTAAACCATGAAAGCAATTTTCGTACTGAAA┆┆╎╎┆┆╎┆┆┆┆┆╎┆╎╎┆┆╎┆╎╎┆┆ 3'
                   ╎┆┆╎┆┆╎╎╎╎┆┆╎╎┆╎╎╎┆┆╎╎╎╎╎╎┆┆╎╎┆╎┆╎╎╎┆┆╎╎┆┆╎┆┆┆┆┆╎┆╎╎┆┆╎┆╎╎┆┆
molde inferior  3' TCCTCCATTTGGTACTTTCGTTAAAAGCATGACTTTCCAACCACCGCGTGAAGGACTTCG 5'
```

- Las ADN polimerasas solo sintetizan en sentido 5'→3' y necesitan un **cebador** de ARN (en
  minúsculas).
- La **cadena líder** crece en el mismo sentido en que avanza la horquilla: es continua y le basta
  un cebador.
- La **cadena rezagada** crece en sentido contrario. Se sintetiza a trozos, los **fragmentos de
  Okazaki** (`◀Ok1`, `◀Ok2`…), cada uno con su cebador. Mientras tanto, el molde que queda al
  descubierto espera protegido por las proteínas SSB (`~`).
- Al terminar, la ADN polimerasa I sustituye los cebadores por ADN y deja mellas (`v ^`), que la
  ADN ligasa sella.
- El resultado son dos moléculas hijas, cada una con una hebra parental y otra nueva (replicación
  **semiconservativa**). El programa comprueba que ambas son idénticas a la original.

| Enzima o proteína | Función |
|---|---|
| DnaA y DnaC | reconocen el origen de replicación y cargan la helicasa |
| Helicasa (DnaB), `H` | rompe los puentes de hidrógeno y separa las hebras |
| ADN girasa (topoisomerasa), `T` | alivia la tensión por delante de la horquilla |
| Proteínas SSB, `~` | mantienen separadas las hebras sencillas |
| Primasa | sintetiza los cebadores de ARN |
| ADN polimerasa III | alarga las hebras nuevas en sentido 5'→3' |
| ADN polimerasa I | sustituye los cebadores de ARN por ADN |
| ADN ligasa | sella las mellas entre fragmentos |

### 2. Transcripción

La ARN polimerasa recorre la hebra molde (la inferior) en sentido 3'→5' y sintetiza el ARNm en
sentido 5'→3'. Aplica la complementariedad A→U, T→A, G→C y C→G. Se ve a mitad de camino, con la
burbuja de transcripción y el ARN recién hecho unido al molde, y al final. El ARNm es igual que la
hebra codificante, pero con U en lugar de T.

### 3. Traducción

Se muestran las tres fases con el ribosoma y sus sitios E, P y A:

- **Iniciación:** el ribosoma se coloca sobre el primer codón AUG y el ARNt iniciador trae la
  formilmetionina (fMet).
- **Elongación:** se dibuja un ciclo completo. Entra el ARNt cuyo anticodón es complementario y
  antiparalelo al codón, se forma el enlace peptídico y el ribosoma avanza. Después se enumeran
  todos los codones leídos.
- **Terminación:** al llegar un codón de terminación (UAA, UAG o UGA), entra un factor de
  liberación y la proteína se suelta.

El resumen final recoge la lectura completa codón a codón.

El gen de ejemplo ([`ejemplos/gen_ejemplo.fasta`](ejemplos/gen_ejemplo.fasta)) codifica el péptido
líder del operón *trp* de *E. coli*: **MKAIFVLKGWWRTS**.

## Simplificaciones

- Se simula una bacteria (*E. coli*) con una sola horquilla. Como su cromosoma es circular, el
  cebador del extremo de cada hebra se rellena desde el otro extremo y no queda ningún hueco.
- Los tamaños se reducen para que quepan en pantalla. Los cebadores tienen 3 nt (unos 10 en la
  célula) y los fragmentos de Okazaki unas decenas de nucleótidos (1000-2000 en *E. coli*).
- La ADN polimerasa I y la ligasa actúan al final. En la célula lo hacen mientras avanza la
  horquilla.
- Se transcribe la secuencia completa, sin buscar promotor ni terminador.
- Se usa el código genético estándar.

## Estructura

```
simulador.py           programa principal: menú, opciones y paso de pantallas
dogma/adn.py           complementariedad, código genético y lectura de secuencias
dogma/replicacion.py   modelo de la replicación (ADN → ADN)
dogma/transcripcion.py modelo de la transcripción (ADN → ARN)
dogma/traduccion.py    modelo de la traducción (ARN → proteína)
dogma/dibujo.py        dibujo de cada paso en la terminal
dogma/imagen.py        salida opcional en imagen PNG (matplotlib)
ejemplos/              gen de ejemplo y su salida
tests/test_dogma.py    pruebas
```

Cada modelo guarda una lista de estados (una "foto" de cada paso con lo que ha ocurrido) y
`dibujo.py` convierte cada estado en una pantalla.

## Pruebas

```bash
python3 -m unittest -v
```

## Autores

Lucía Hernández Nieves, Darío José Iglesias San José, Jesús Varela Vasquez
