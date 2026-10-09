"""
Reglas de complementariedad, código genético y lectura de secuencias.

Convención: la hebra superior del ADN se escribe 5'→3' de izquierda a derecha y
su complementaria (hebra inferior) se escribe debajo, alineada, así que se lee
3'→5'. Las dos hebras son antiparalelas.
"""

import random

# Complementariedad de bases
COMPLEMENTO_ADN = {"A": "T", "T": "A", "G": "C", "C": "G"}
ADN_A_ARN = {"A": "U", "T": "A", "G": "C", "C": "G"}      # base del molde -> base del ARN
COMPLEMENTO_ARN = {"A": "U", "U": "A", "G": "C", "C": "G"}
PUENTES_H = {"A": 2, "T": 2, "U": 2, "G": 3, "C": 3}       # A=T forma 2 puentes; G≡C, 3


def complementaria(adn):
    """Hebra complementaria, alineada debajo de la original (se lee 3'→5')."""
    return "".join(COMPLEMENTO_ADN[base] for base in adn)


def arn_complementario(molde):
    """ARN que se aparea con un molde de ADN: A→U, T→A, G→C, C→G."""
    return "".join(ADN_A_ARN[base] for base in molde)


def anticodon(codon):
    """Anticodón del ARNt, escrito 3'→5' para que quede alineado con el codón."""
    return "".join(COMPLEMENTO_ARN[base] for base in codon)


# Código genético estándar. Las bases de cada posición van en el orden U, C, A, G.
BASES = "UCAG"
TABLA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODIGO_GENETICO = {
    b1 + b2 + b3: TABLA[16 * i + 4 * j + k]
    for i, b1 in enumerate(BASES)
    for j, b2 in enumerate(BASES)
    for k, b3 in enumerate(BASES)
}
CODONES_STOP = ("UAA", "UAG", "UGA")

AMINOACIDOS = {   # código de una letra -> (abreviatura, nombre)
    "A": ("Ala", "alanina"), "R": ("Arg", "arginina"), "N": ("Asn", "asparagina"),
    "D": ("Asp", "ácido aspártico"), "C": ("Cys", "cisteína"), "Q": ("Gln", "glutamina"),
    "E": ("Glu", "ácido glutámico"), "G": ("Gly", "glicina"), "H": ("His", "histidina"),
    "I": ("Ile", "isoleucina"), "L": ("Leu", "leucina"), "K": ("Lys", "lisina"),
    "M": ("Met", "metionina"), "F": ("Phe", "fenilalanina"), "P": ("Pro", "prolina"),
    "S": ("Ser", "serina"), "T": ("Thr", "treonina"), "W": ("Trp", "triptófano"),
    "Y": ("Tyr", "tirosina"), "V": ("Val", "valina"),
}


# --- Entrada de datos ------------------------------------------------------------

def limpiar(texto):
    """Quita espacios y números, pasa a mayúsculas y comprueba que solo haya A, C, G y T."""
    adn = "".join(c for c in texto.upper() if not c.isspace() and not c.isdigit())
    no_validos = sorted(set(adn) - set("ACGT"))
    if no_validos:
        raise ValueError(f"caracteres no válidos ({', '.join(no_validos)}): el ADN solo puede "
                         f"contener A, C, G y T")
    if len(adn) < 12:
        raise ValueError("la secuencia debe tener al menos 12 pares de bases")
    return adn


def leer_fasta(ruta):
    """Secuencia del primer registro de un archivo FASTA."""
    lineas = []
    with open(ruta, encoding="utf-8") as archivo:
        for linea in archivo:
            linea = linea.strip()
            if linea.startswith(">"):
                if lineas:
                    break          # solo se lee el primer registro
                continue
            lineas.append(linea)
    return "".join(lineas)


def gen_aleatorio(longitud, semilla=None):
    """Secuencia aleatoria con forma de gen: inicio sin ATG + ATG + codones + stop + final."""
    rng = random.Random(semilla)
    sin_stop = [c.replace("U", "T") for c, aa in CODIGO_GENETICO.items() if aa != "*"]
    while True:
        inicio = "".join(rng.choice("ACGT") for _ in range(6))
        if "ATG" not in inicio:
            break
    codones = (longitud - 6 - 3) // 3                       # deja 3 nt o más al final
    gen = inicio + "ATG" + "".join(rng.choice(sin_stop) for _ in range(codones - 2))
    gen += rng.choice(["TAA", "TAG", "TGA"])
    return gen + "".join(rng.choice("ACGT") for _ in range(longitud - len(gen)))
