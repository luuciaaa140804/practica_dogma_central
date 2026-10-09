"""
Replicación del ADN (ADN → ADN) en una bacteria (E. coli).

Una horquilla de replicación avanza de izquierda a derecha:

    molde superior  5' ─────────────► 3'   su hebra nueva crece hacia la izquierda (←):
                                           CADENA REZAGADA, a trozos (fragmentos de Okazaki)
    molde inferior  3' ◄───────────── 5'   su hebra nueva crece hacia la derecha (→):
                                           CADENA LÍDER, continua

Las ADN polimerasas solo alargan en sentido 5'→3' y necesitan un cebador de ARN
que les dé un extremo 3'-OH. La hebra nueva que copia el molde inferior crece en
el mismo sentido en que avanza la horquilla, así que basta un cebador. La que
copia el molde superior crece en sentido contrario: cada vez que la helicasa deja
al descubierto un tramo de molde, la primasa pone un cebador junto a la horquilla
y la polimerasa sintetiza un fragmento de Okazaki hacia atrás.

El cromosoma bacteriano es circular, así que el final de la molécula está unido
a su principio: el cebador del extremo de cada hebra se rellena desde el otro
extremo y no queda ningún hueco.
"""

from dataclasses import dataclass, field
from typing import List, Optional

from .adn import ADN_A_ARN, COMPLEMENTO_ADN, PUENTES_H, complementaria

# Tipos de nucleótido de las hebras nuevas
CEBADOR = "cebador"   # ARN puesto por la primasa
LIDER = "lider"       # ADN de la cadena líder (ADN polimerasa III)
OKAZAKI = "okazaki"   # ADN de un fragmento de Okazaki (ADN polimerasa III)
RELLENO = "relleno"   # ADN que sustituye a un cebador (ADN polimerasa I)


@dataclass
class Nucleotido:
    base: str
    tipo: str
    fragmento: int      # número de fragmento de Okazaki (0 en la cadena líder)


@dataclass
class Estado:
    """Foto de un paso de la replicación."""
    titulo: str
    horquilla: int                          # posiciones [0, horquilla) ya abiertas
    rezagada: List[Optional[Nucleotido]]    # hebra nueva sobre el molde superior (3'←5')
    lider: List[Optional[Nucleotido]]       # hebra nueva sobre el molde inferior (5'→3')
    eventos: List[str]
    mellas: List[tuple] = field(default_factory=list)   # (hebra, posición tras la mella)


@dataclass
class Replicacion:
    superior: str
    inferior: str
    fragmentos: List[tuple]     # (inicio, fin) de cada fragmento de Okazaki
    cebador: int                # longitud de los cebadores
    estados: List[Estado]

    def hijas(self):
        """Las dos moléculas hijas: cada una conserva una hebra parental (semiconservativa)."""
        final = self.estados[-1]
        nueva_sup = "".join(nt.base for nt in final.rezagada)
        nueva_inf = "".join(nt.base for nt in final.lider)
        return (self.superior, nueva_sup), (nueva_inf, self.inferior)

    def es_correcta(self):
        """Hebras nuevas completas, sin ARN, sin mellas y complementarias a su molde."""
        final = self.estados[-1]
        for molde, hebra in ((self.superior, final.rezagada), (self.inferior, final.lider)):
            for base, nt in zip(molde, hebra):
                if nt is None or nt.tipo == CEBADOR or COMPLEMENTO_ADN[base] != nt.base:
                    return False
        return not final.mellas


def trocear(n, longitud, cebador):
    """Límites (inicio, fin) de los fragmentos de Okazaki que cubren la molécula."""
    limites = list(range(0, n, longitud)) + [n]
    if len(limites) > 2 and limites[-1] - limites[-2] <= cebador:
        del limites[-2]          # un último trozo más corto que el cebador se une al anterior
    return list(zip(limites, limites[1:]))


def replicar(adn, avance=None, long_okazaki=None, cebador=3):
    """
    Simula la replicación y devuelve todos sus pasos. Los tamaños están reducidos para
    que quepan en pantalla (en E. coli los fragmentos de Okazaki miden 1000-2000 nt).
    """
    n = len(adn)
    superior, inferior = adn, complementaria(adn)
    avance = avance or max(4, round(n / 5))       # pb que abre la helicasa en cada paso
    long_okazaki = long_okazaki or 2 * avance
    cebador = min(cebador, long_okazaki - 1)
    fragmentos = trocear(n, long_okazaki, cebador)
    rezagada = [None] * n
    lider = [None] * n
    estados = []

    def foto(titulo, horquilla, eventos, mellas=()):
        estados.append(Estado(titulo, horquilla, list(rezagada), list(lider), eventos,
                              list(mellas)))

    horquilla = 0
    hechos = 0      # fragmentos de Okazaki ya sintetizados
    while horquilla < n:
        antes = horquilla
        horquilla = min(n, horquilla + avance)
        rotos = sum(PUENTES_H[base] for base in adn[antes:horquilla])
        eventos = []
        if antes == 0:
            eventos.append("La proteína DnaA reconoce el origen de replicación (oriC). La ADN "
                           "girasa alivia la tensión por delante de la horquilla y las proteínas "
                           "SSB (~) evitan que las hebras separadas se vuelvan a unir.")
        eventos.append(f"Helicasa: separa las posiciones {antes + 1}-{horquilla} ({rotos} "
                       f"puentes de hidrógeno rotos).")

        # Cadena líder: un solo cebador al principio y síntesis continua hasta la horquilla
        if antes == 0:
            for i in range(cebador):
                lider[i] = Nucleotido(ADN_A_ARN[inferior[i]], CEBADOR, 0)
        nuevos = [i for i in range(antes, horquilla) if lider[i] is None]
        for i in nuevos:
            lider[i] = Nucleotido(COMPLEMENTO_ADN[inferior[i]], LIDER, 0)
        if antes == 0:
            arn = "".join(nt.base for nt in lider[:cebador])
            eventos.append(f"Líder: la primasa pone un cebador de ARN (5'-{arn}-3') y la ADN "
                           f"polimerasa III lo alarga de forma continua (→) hasta la horquilla.")
        else:
            eventos.append(f"Líder: la ADN polimerasa III sigue de forma continua hasta la "
                           f"horquilla (+{len(nuevos)} nt).")

        # Cadena rezagada: un fragmento de Okazaki cuando hay suficiente molde al descubierto
        antes_hechos = hechos
        while hechos < len(fragmentos) and fragmentos[hechos][1] <= horquilla:
            inicio, fin = fragmentos[hechos]
            hechos += 1
            for i in range(fin - cebador, fin):
                rezagada[i] = Nucleotido(ADN_A_ARN[superior[i]], CEBADOR, hechos)
            for i in range(inicio, fin - cebador):
                rezagada[i] = Nucleotido(COMPLEMENTO_ADN[superior[i]], OKAZAKI, hechos)
            eventos.append(f"Rezagada: la primasa pone un cebador junto a la horquilla y la ADN "
                           f"polimerasa III sintetiza hacia atrás (←) el fragmento de Okazaki "
                           f"{hechos} (posiciones {inicio + 1}-{fin}).")
        if hechos == antes_hechos:
            eventos.append("Rezagada: el molde recién abierto espera, cubierto por SSB, al "
                           "siguiente fragmento de Okazaki.")
        titulo = "Iniciación: se abre la horquilla" if antes == 0 else "La horquilla avanza"
        foto(titulo, horquilla, eventos)

    # Maduración: la ADN polimerasa I cambia cada cebador por ADN y deja una mella
    mellas = []
    for i in range(cebador):
        lider[i] = Nucleotido(COMPLEMENTO_ADN[inferior[i]], RELLENO, 0)
    mellas.append(("lider", cebador))
    for numero, (inicio, fin) in enumerate(fragmentos, start=1):
        for i in range(fin - cebador, fin):
            rezagada[i] = Nucleotido(COMPLEMENTO_ADN[superior[i]], RELLENO, numero)
        mellas.append(("rezagada", fin - cebador))
    foto("La ADN polimerasa I sustituye los cebadores", n, [
        "La ADN polimerasa I cambia cada cebador de ARN por ADN (actividad exonucleasa 5'→3'). "
        "Como el cromosoma bacteriano es circular, el cebador de cada extremo se rellena desde "
        "el otro extremo.",
        f"Quedan {len(mellas)} mellas (v ^): falta el enlace fosfodiéster entre dos trozos de ADN "
        f"vecinos."], mellas)
    foto("La ADN ligasa sella las mellas", n, [
        f"La ADN ligasa sella las {len(mellas)} mellas: cada hebra nueva es ya una sola cadena "
        f"continua."])

    return Replicacion(superior, inferior, fragmentos, cebador, estados)
