"""
Traducción (ARN → proteína).

El ribosoma se coloca sobre el primer codón AUG del ARNm y lo lee de tres en
tres bases (codones) en sentido 5'→3'. Para cada codón entra en el sitio A el
ARNt cuyo anticodón es complementario y antiparalelo (5'-AUG-3' / 3'-UAC-5'),
se forma el enlace peptídico y el ribosoma avanza un codón. Al llegar un codón
de terminación (UAA, UAG o UGA) entra un factor de liberación y la proteína se
suelta. En bacterias el primer aminoácido es la formilmetionina (fMet).
"""

from dataclasses import dataclass
from typing import List, Optional

from .adn import AMINOACIDOS, CODIGO_GENETICO, CODONES_STOP, anticodon


@dataclass
class Estado:
    titulo: str
    fase: str            # "inicio", "elongacion", "terminacion", "sin_stop" o "sin_inicio"
    sitio_p: int         # índice del codón que ocupa el sitio P del ribosoma
    peptido: str         # aminoácidos unidos hasta ahora (código de una letra)
    eventos: List[str]


@dataclass
class Traduccion:
    arnm: str
    inicio: int               # posición del AUG (-1 si no hay)
    codones: List[str]        # del AUG al codón stop (incluido)
    proteina: str             # código de una letra
    stop: Optional[str]
    estados: List[Estado]

    def tres_letras(self):
        nombres = [AMINOACIDOS[aa][0] for aa in self.proteina]
        if nombres:
            nombres[0] = "fMet"
        return nombres


def traducir(arnm):
    inicio = arnm.find("AUG")
    if inicio == -1:
        estado = Estado("Sin codón de inicio", "sin_inicio", 0, "", [
            "El ARNm no tiene ningún codón AUG: el ribosoma no puede empezar y no se sintetiza "
            "ninguna proteína."])
        return Traduccion(arnm, -1, [], "", None, [estado])

    codones, stop = [], None
    for i in range(inicio, len(arnm) - 2, 3):
        codones.append(arnm[i:i + 3])
        if codones[-1] in CODONES_STOP:
            stop = codones[-1]
            break
    proteina = "".join(CODIGO_GENETICO[c] for c in codones if c not in CODONES_STOP)

    estados = [Estado("Iniciación", "inicio", 0, proteina[:1], [
        f"La subunidad pequeña (30S) se une al ARNm y se coloca sobre el primer AUG (posición "
        f"{inicio + 1}). El ARNt iniciador (anticodón 3'-UAC-5') trae la formilmetionina (fMet) "
        f"al sitio P y se une la subunidad grande (50S)."])]

    if len(proteina) >= 2:
        # Se dibuja un ciclo de elongación; con 3 aminoácidos ya se ven ocupados E, P y A
        p = 1 if len(proteina) >= 3 else 0
        codon = codones[p + 1]
        leidos = [f"{c}→{AMINOACIDOS[CODIGO_GENETICO[c]][0]}" for c in codones[1:len(proteina)]]
        if len(leidos) > 16:
            leidos = leidos[:12] + [f"… ({len(leidos)} codones en total)"]
        estados.append(Estado("Elongación", "elongacion", p, proteina[:p + 2], [
            f"Entra en el sitio A el ARNt cuyo anticodón (3'-{anticodon(codon)}-5') es "
            f"complementario del codón (5'-{codon}-3'). Se forma el enlace peptídico y el "
            f"ribosoma avanza un codón; el ARNt vacío sale por el sitio E.",
            "El ciclo se repite con cada codón: " + ", ".join(leidos) + "."]))

    if stop:
        estados.append(Estado("Terminación", "terminacion", len(proteina) - 1, proteina, [
            f"Ningún ARNt reconoce el codón {stop}: entra un factor de liberación, la proteína "
            f"de {len(proteina)} aminoácidos se suelta y el ribosoma se separa."]))
    else:
        estados.append(Estado("El ARNm se acaba sin codón de terminación", "sin_stop",
                              len(proteina) - 1, proteina, [
            "El ribosoma llega al final del ARNm sin encontrar un codón de terminación: la "
            "proteína queda incompleta."]))
    return Traduccion(arnm, inicio, codones, proteina, stop, estados)
