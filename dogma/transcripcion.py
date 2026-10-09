"""
Transcripción (ADN → ARN).

La ARN polimerasa recorre la hebra molde (la inferior) en sentido 3'→5' y
sintetiza el ARN en sentido 5'→3', emparejando cada base del molde con su
complementaria: A→U, T→A, G→C y C→G. Por eso el ARNm tiene la misma secuencia
que la hebra codificante (la superior), pero con U en lugar de T.

Simplificación: se transcribe la secuencia completa. En la célula, la ARN
polimerasa empieza en un promotor y se detiene en un terminador.
"""

from dataclasses import dataclass
from typing import List

from .adn import arn_complementario, complementaria


@dataclass
class Estado:
    titulo: str
    posicion: int        # nucleótidos transcritos hasta ahora
    eventos: List[str]


@dataclass
class Transcripcion:
    codificante: str     # hebra superior, 5'→3'
    molde: str           # hebra inferior, 3'→5'
    arnm: str            # 5'→3'
    estados: List[Estado]


def transcribir(adn):
    molde = complementaria(adn)
    arnm = arn_complementario(molde)
    n = len(adn)
    estados = [
        Estado("Iniciación y elongación", n // 2, [
            "La ARN polimerasa (con el factor σ) se une al principio del gen y abre la doble "
            "hélice; a diferencia de la ADN polimerasa, no necesita cebador.",
            "Lee la hebra molde (la inferior) en sentido 3'→5' y sintetiza el ARN 5'→3' por "
            "complementariedad: A→U, T→A, G→C y C→G."]),
        Estado("Terminación", n, [
            f"La ARN polimerasa llega al final del gen y suelta el ARNm de {n} nucleótidos: es "
            f"igual que la hebra codificante, pero con U en lugar de T."]),
    ]
    return Transcripcion(adn, molde, arnm, estados)
