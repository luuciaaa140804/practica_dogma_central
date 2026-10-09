"""
Dibujo en la terminal de cada paso de la simulación.

Cada pantalla es una lista de líneas de texto. Las secuencias se dibujan
alineadas columna a columna (una columna = una base) para que se vea la
complementariedad; si no caben a lo ancho, se parten en bloques.
"""

import textwrap

from .adn import (AMINOACIDOS, CODIGO_GENETICO, CODONES_STOP, PUENTES_H, anticodon,
                  complementaria)
from .replicacion import CEBADOR, LIDER, RELLENO

USAR_COLOR = True      # simulador.py lo desactiva si la salida no es una terminal

COLORES = {            # códigos de color ANSI
    "titulo": "1;97;44", "seccion": "1;96", "tenue": "2", "negrita": "1", "ok": "1;92",
    "error": "1;91", "molde": "1", "lider": "92", "okazaki1": "96", "okazaki2": "94",
    "cebador": "1;95", "relleno": "93", "puente": "2", "ssb": "36", "mella": "1;91",
    "helicasa": "1;91", "girasa": "1;35", "arn": "1;93", "polimerasa": "1;93",
    "inicio": "1;30;42", "stop": "1;97;41", "codon1": "93", "codon2": "33",
    "anticodon": "95", "aminoacido": "1;92", "ribosoma": "1;96",
}


def pintar(texto, estilo):
    if not USAR_COLOR or not estilo:
        return texto
    return f"\033[{COLORES[estilo]}m{texto}\033[0m"


def celdas(texto, estilo=None):
    """Una celda por carácter: cada celda ocupa una columna de la pantalla."""
    return [pintar(c, estilo) for c in texto]


def cabecera(texto, ancho):
    return [pintar(f" {texto} ".ljust(ancho), "titulo"), ""]


def lista(eventos, ancho):
    lineas = []
    for texto in eventos:
        trozos = textwrap.wrap(texto, ancho - 4)
        lineas.append(pintar(" • ", "seccion") + trozos[0])
        lineas += ["   " + trozo for trozo in trozos[1:]]
    return lineas


def regla(n):
    """Fila con el número de posición (1, 11, 21...) encima de cada décima base."""
    fila = [" "] * n
    for i in range(0, n, 10):
        numero = str(i + 1)
        if i + len(numero) <= n:
            fila[i:i + len(numero)] = celdas(numero, "tenue")
    return fila


def rejilla(filas, ancho):
    """Filas (nombre, extremo izquierdo, celdas, extremo derecho) alineadas columna a columna."""
    n = max(len(fila[2]) for fila in filas)
    sitio = max(20, ancho - 23)
    lineas = []
    for inicio in range(0, n, sitio):
        fin = min(n, inicio + sitio)
        for nombre, izquierdo, contenido, derecho in filas:
            linea = pintar(f"{nombre:<16}", "tenue")
            linea += pintar(f"{izquierdo if inicio == 0 else '':>2} ", "negrita")
            linea += "".join(contenido[inicio:fin])
            if fin == n and derecho:
                linea += pintar(" " + derecho, "negrita")
            lineas.append(linea)
        lineas.append("")
    return lineas


def puente(base):
    """╎ = 2 puentes de hidrógeno (A=T); ┆ = 3 (G≡C)."""
    return "┆" if PUENTES_H[base] == 3 else "╎"


def marca(correcto):
    return pintar("✔", "ok") if correcto else pintar("✘", "error")


# --- Replicación ------------------------------------------------------------------

def celda_nucleotido(nt):
    if nt is None:
        return " "
    if nt.tipo == CEBADOR:
        return pintar(nt.base.lower(), "cebador")
    if nt.tipo == RELLENO:
        return pintar(nt.base, "relleno")
    if nt.tipo == LIDER:
        return pintar(nt.base, "lider")
    return pintar(nt.base, "okazaki1" if nt.fragmento % 2 else "okazaki2")


def corchete(ancho, etiqueta):
    """├──etiqueta──┤ que ocupa 'ancho' columnas (sin etiqueta si no cabe)."""
    if ancho < 2:
        return "│" * ancho
    if len(etiqueta) > ancho - 2:
        etiqueta = ""
    return "├" + etiqueta.center(ancho - 2, "─") + "┤"


ETAPA_REPLICACION = "1 · REPLICACIÓN (ADN → ADN)"


def pantalla_replicacion(rep, estado, paso, total, ancho):
    sup, inf, n = rep.superior, rep.inferior, len(rep.superior)
    lineas = cabecera(f"{ETAPA_REPLICACION} · Paso {paso} de {total} · {estado.titulo}", ancho)
    enlace_sup, rezagada, medio, lider, enlace_inf = ([" "] * n for _ in range(5))
    for i in range(n):
        if i < estado.horquilla:       # zona abierta: cada molde con su hebra nueva
            nt_sup, nt_inf = estado.rezagada[i], estado.lider[i]
            enlace_sup[i] = pintar(puente(sup[i]), "puente") if nt_sup else pintar("~", "ssb")
            rezagada[i] = celda_nucleotido(nt_sup)
            lider[i] = celda_nucleotido(nt_inf)
            enlace_inf[i] = pintar(puente(sup[i]), "puente") if nt_inf else pintar("~", "ssb")
        else:                          # zona sin abrir: las hebras siguen unidas
            for fila in (enlace_sup, rezagada, medio, lider, enlace_inf):
                fila[i] = pintar(puente(sup[i]), "puente")
    for numero, (a, b) in enumerate(rep.fragmentos, start=1):
        if estado.rezagada[a] is not None:
            estilo = "okazaki1" if numero % 2 else "okazaki2"
            medio[a:b] = celdas(corchete(b - a, f"◀Ok{numero}"), estilo)
    if estado.horquilla < n:
        medio[estado.horquilla] = pintar("H", "helicasa")
        if estado.horquilla + 3 < n:
            medio[estado.horquilla + 3] = pintar("T", "girasa")
    for hebra, posicion in estado.mellas:
        if hebra == "rezagada":
            enlace_sup[posicion] = pintar("v", "mella")
        else:
            enlace_inf[posicion] = pintar("^", "mella")
    lineas += rejilla([("posición", "", regla(n), ""),
                       ("molde superior", "5'", celdas(sup, "molde"), "3'"),
                       ("", "", enlace_sup, ""),
                       ("rezagada  ←", "3'", rezagada, "5'"),
                       ("", "", medio, ""),
                       ("líder     →", "5'", lider, "3'"),
                       ("", "", enlace_inf, ""),
                       ("molde inferior", "3'", celdas(inf, "molde"), "5'")], ancho)
    lineas += lista(estado.eventos, ancho)
    if paso == 1:
        lineas += ["", leyenda_replicacion()]
    return lineas


def leyenda_replicacion():
    return (f"  {pintar('ATGC', 'molde')} molde   {pintar('ATGC', 'lider')} líder   "
            f"{pintar('AT', 'okazaki1')}{pintar('GC', 'okazaki2')} Okazaki   "
            f"{pintar('augc', 'cebador')} cebador (ARN)   {pintar('ATGC', 'relleno')} relleno de "
            f"la pol I\n  {pintar('H', 'helicasa')} helicasa   {pintar('T', 'girasa')} girasa   "
            f"{pintar('~', 'ssb')} SSB   {pintar('v ^', 'mella')} mella   "
            f"╎ ┆ 2 / 3 puentes de hidrógeno")


def pantalla_hijas(rep, ancho):
    """Último paso (la ligasa sella las mellas) y las dos moléculas hijas resultantes."""
    (sup, nueva_sup), (nueva_inf, inf) = rep.hijas()
    n, final, total = len(sup), rep.estados[-1], len(rep.estados)
    lineas = cabecera(f"{ETAPA_REPLICACION} · Paso {total} de {total} · {final.titulo}: dos "
                      f"moléculas hijas", ancho)
    puentes = [pintar(puente(b), "puente") for b in sup]
    lineas += rejilla([("posición", "", regla(n), ""),
                       ("hija 1: molde", "5'", celdas(sup, "molde"), "3'"),
                       ("", "", puentes, ""),
                       ("        nueva", "3'", [celda_nucleotido(nt) for nt in final.rezagada], "5'"),
                       ("", "", [" "] * n, ""),
                       ("hija 2: nueva", "5'", [celda_nucleotido(nt) for nt in final.lider], "3'"),
                       ("", "", puentes, ""),
                       ("        molde", "3'", celdas(inf, "molde"), "5'")], ancho)
    identicas = nueva_sup == rep.inferior and nueva_inf == rep.superior
    largo = rep.fragmentos[0][1] - rep.fragmentos[0][0]
    return lineas + lista(final.eventos, ancho) + [
        f"  {marca(rep.es_correcta())} Hebras nuevas complementarias a su molde (A-T, G-C), sin "
        f"cebadores ni mellas.",
        f"  {marca(identicas)} Semiconservativa: cada hija conserva una hebra parental y ambas "
        f"son idénticas a la original.",
        f"  Fragmentos de Okazaki: {len(rep.fragmentos)} de {largo} nt (en E. coli miden "
        f"1000-2000 nt). Cebadores: {len(rep.fragmentos) + 1}."]


# --- Transcripción ----------------------------------------------------------------

def pantalla_transcripcion(tx, estado, paso, total, ancho):
    n, p = len(tx.codificante), estado.posicion
    lineas = cabecera(f"2 · TRANSCRIPCIÓN (ADN → ARN) · Paso {paso} de {total} · "
                      f"{estado.titulo}", ancho)
    apareadas = [pintar(puente(b), "puente") for b in tx.codificante]
    union, arn, polimerasa = [" "] * n, [" "] * n, [" "] * n
    for i in range(p):
        arn[i] = pintar(tx.arnm[i], "arn")
    if p < n:     # la ARN polimerasa sigue sobre el ADN: burbuja de transcripción
        a, b = max(0, p - 8), min(n, max(p + 3, 10))
        for i in range(a, b):
            apareadas[i] = " "                      # las dos hebras de ADN, separadas
        for i in range(max(0, p - 6), p):
            union[i] = pintar(puente(tx.codificante[i]), "puente")   # ARN unido al molde
        arn[p] = pintar("▶", "polimerasa")
        polimerasa[a:b] = celdas(corchete(b - a, "ARN pol").replace("├", "╰").replace("┤", "╯"),
                                 "polimerasa")
    lineas += rejilla([("posición", "", regla(n), ""),
                       ("codificante", "5'", celdas(tx.codificante, "tenue"), "3'"),
                       ("", "", apareadas, ""),
                       ("molde", "3'", celdas(tx.molde, "molde"), "5'"),
                       ("", "", union, ""),
                       ("ARNm", "5'", arn, "3'"),
                       ("", "", polimerasa, "")], ancho)
    return lineas + lista(estado.eventos, ancho)


# --- Traducción -------------------------------------------------------------------

ANCHO_CODON = 6


def triplete(tr, j):
    """Tres bases del ARNm en el marco de lectura del AUG (j = 0 es el propio AUG)."""
    a = tr.inicio + 3 * j
    return "".join(tr.arnm[k] if 0 <= k < len(tr.arnm) else " " for k in range(a, a + 3))


def pantalla_traduccion(tr, estado, paso, total, ancho):
    lineas = cabecera(f"3 · TRADUCCIÓN (ARN → PROTEÍNA) · Paso {paso} de {total} · "
                      f"{estado.titulo}", ancho)
    if estado.fase == "sin_inicio":
        return lineas + lista(estado.eventos, ancho)

    # El ribosoma queda fijo y el ARNm se desplaza: se ven 8 tripletes y los
    # huecos 2, 3 y 4 son los sitios E, P y A.
    p = estado.sitio_p
    huecos = [j for j in range(p - 3, p + 5) if j <= p + 1 or triplete(tr, j).strip()]
    prefijo = "  ARNm  5'-" + ("…" if tr.inicio + 3 * (p - 3) > 0 else "")
    tripletes = []
    for j in huecos:
        if j == 0:
            estilo = "inicio"
        elif 0 < j < len(tr.codones):
            estilo = "stop" if tr.codones[j] in CODONES_STOP else ("codon1" if j % 2 else "codon2")
        else:
            estilo = "tenue"             # regiones que no se traducen
        texto = triplete(tr, j)
        tripletes.append(pintar(texto.rstrip() if j == huecos[-1] else texto, estilo))
    fila_arnm = prefijo + (" " * (ANCHO_CODON - 3)).join(tripletes)
    fila_arnm += ("…" if tr.inicio + 3 * (p + 5) < len(tr.arnm) else "") + "-3'"

    x0 = len(prefijo) + 2 * ANCHO_CODON - 1          # borde izquierdo del ribosoma
    caja = 2 * ANCHO_CODON + 5
    sitios = ["│"] + [" "] * (caja - 2) + ["│"]
    for k, letra in enumerate("EPA"):
        sitios[2 + k * ANCHO_CODON] = letra
    lineas += [" " * x0 + pintar("┌" + "─" * (caja - 2) + "┐  ribosoma 70S", "ribosoma"),
               " " * x0 + pintar("".join(sitios), "ribosoma"),
               fila_arnm]

    # Qué hay en cada sitio: ARNt (con su anticodón) o factor de liberación
    ocupantes = {"P": "arnt"}
    if estado.fase in ("elongacion", "terminacion"):
        ocupantes["A"] = "arnt" if estado.fase == "elongacion" else "factor"
        if p >= 1:
            ocupantes["E"] = "sale"
    filas = {"enlaces": {}, "anticodon": {}, "tallo": {}, "aa": {}}
    for k, sitio in enumerate("EPA"):
        if sitio not in ocupantes:
            continue
        indice = p - 1 + k
        codon = tr.codones[indice]
        if ocupantes[sitio] == "factor":
            filas["anticodon"][k] = ("RF", "ribosoma")
            filas["aa"][k] = ("fin", "stop")
            continue
        if ocupantes[sitio] == "arnt":
            filas["enlaces"][k] = ("".join(puente(b) for b in codon), "puente")
        filas["anticodon"][k] = (anticodon(codon), "anticodon")
        filas["tallo"][k] = (" │", "anticodon")
        if ocupantes[sitio] == "sale":
            filas["aa"][k] = ("sale", "tenue")
        else:
            nombre = "fMet" if indice == 0 else AMINOACIDOS[CODIGO_GENETICO[codon]][0]
            filas["aa"][k] = (nombre, "aminoacido")
    notas = {"anticodon": "  anticodón del ARNt (3'→5')", "aa": "  aminoácido"}
    for clave, contenido in filas.items():
        fila = [" "] * (x0 + caja)
        for k, (texto, estilo) in contenido.items():
            columna = x0 + 1 + k * ANCHO_CODON
            fila[columna:columna + len(texto)] = celdas(texto, estilo)
        lineas.append("".join(fila) + pintar(notas.get(clave, ""), "tenue"))

    nombres = [AMINOACIDOS[aa][0] for aa in estado.peptido]
    nombres[0] = "fMet"
    trozos = textwrap.wrap(" ".join(nombres), ancho - 26)      # una línea por trozo de cadena
    lineas.append("")
    for k, trozo in enumerate(trozos):
        enlace = "─" if k < len(trozos) - 1 else ""
        inicio = "  Cadena polipeptídica: " if k == 0 else " " * 24
        lineas.append(inicio + pintar(trozo.replace(" ", "─") + enlace, "aminoacido"))
    lineas.append("")
    lineas += lista(estado.eventos, ancho)
    if estado.fase == "inicio":
        lineas += ["", f"  {pintar('AUG', 'inicio')} inicio   {pintar('UAA', 'stop')} "
                       f"fin   {pintar('ACC', 'tenue')} no se traduce   "
                       f"{pintar('E P A', 'ribosoma')} sitios del ribosoma   "
                       f"{pintar('RF', 'ribosoma')} factor de liberación"]
    return lineas


# --- Portada y resumen ---------------------------------------------------------------

def portada(adn, descripcion, ancho):
    """Datos de la secuencia y la molécula de ADN de partida."""
    n = len(adn)
    lineas = [pintar(" SIMULADOR DEL DOGMA CENTRAL DE LA BIOLOGÍA MOLECULAR".ljust(ancho), "titulo"),
              pintar(" ADN → ADN (replicación) · ADN → ARN (transcripción) · ARN → proteína "
                     "(traducción)", "seccion"),
              "",
              f"  {descripcion} · {n} pares de bases · bacteria (E. coli)",
              ""]
    lineas += rejilla([("posición", "", regla(n), ""),
                       ("hebra superior", "5'", celdas(adn, "molde"), "3'"),
                       ("", "", [pintar(puente(b), "puente") for b in adn], ""),
                       ("hebra inferior", "3'", celdas(complementaria(adn), "molde"), "5'")],
                      ancho)
    puentes = sum(PUENTES_H[b] for b in adn)
    return lineas + lista([f"Las dos hebras son antiparalelas y complementarias: A=T (2 puentes de "
                           f"hidrógeno) y G≡C (3), {puentes} puentes en total.",
                           "Los fragmentos de Okazaki y los cebadores se dibujan más cortos que en "
                           "la célula para que quepan en pantalla."], ancho)


def pantalla_resumen(rep, tx, tr, ancho):
    lineas = cabecera("RESUMEN · FLUJO DE LA INFORMACIÓN GENÉTICA: ADN → ADN → ARN → PROTEÍNA",
                      ancho)
    proteina = (f"proteína de {len(tr.proteina)} aminoácidos" if tr.proteina
                else "no hay codón de inicio: no se sintetiza proteína")
    lineas += [f"  ADN → ADN    {marca(rep.es_correcta())} replicación semiconservativa: dos "
               f"moléculas idénticas a la original",
               f"  ADN → ARN    {marca(bool(tx.arnm))} transcripción: ARNm de {len(tx.arnm)} "
               f"nucleótidos",
               f"  ARN → prot.  {marca(bool(tr.proteina))} traducción: {proteina}", ""]
    if not tr.codones:
        return lineas

    # Del ADN a la proteína codón a codón (5 columnas por codón)
    filas = {clave: [] for clave in ("cod", "mol", "arn", "anti", "aa", "letra")}

    def poner(clave, texto, estilo, sangria=0):
        filas[clave] += [" "] * sangria + celdas(texto, estilo) + [" "] * (5 - sangria - len(texto))

    for j, codon in enumerate(tr.codones):
        dna = codon.replace("U", "T")
        poner("cod", dna, "tenue")
        poner("mol", complementaria(dna), "molde")
        if codon in CODONES_STOP:
            poner("arn", codon, "stop")
            poner("anti", "", None)
            poner("aa", "Stop", "stop")
            poner("letra", "*", "stop", sangria=1)
        else:
            aa = CODIGO_GENETICO[codon]
            poner("arn", codon, "inicio" if j == 0 else ("codon1" if j % 2 else "codon2"))
            poner("anti", anticodon(codon), "anticodon")
            poner("aa", "fMet" if j == 0 else AMINOACIDOS[aa][0], "aminoacido")
            poner("letra", aa, "aminoacido", sangria=1)
    lineas.append(pintar(" Del gen a la proteína, codón a codón", "seccion"))
    lineas += rejilla([("ADN codificante", "5'", filas["cod"], "3'"),
                       ("ADN molde", "3'", filas["mol"], "5'"),
                       ("ARNm", "5'", filas["arn"], "3'"),
                       ("anticodón ARNt", "3'", filas["anti"], "5'"),
                       ("aminoácido", "", filas["aa"], ""),
                       ("", "", filas["letra"], "")], ancho)
    for texto in ("-".join(tr.tres_letras()), tr.proteina):
        for k, trozo in enumerate(textwrap.wrap(texto, ancho - 14)):
            inicio = "  Proteína: " if k == 0 and texto != tr.proteina else " " * 12
            lineas.append(inicio + pintar(trozo, "aminoacido"))
    if not tr.stop:
        lineas.append("  (incompleta: el ARNm no tiene codón de terminación)")
    return lineas


def todas(adn, descripcion, rep, tx, tr, ancho):
    """Todas las pantallas de la simulación, en orden."""
    pantallas = [portada(adn, descripcion, ancho)]
    pantallas += [pantalla_replicacion(rep, e, i, len(rep.estados), ancho)
                  for i, e in enumerate(rep.estados[:-1], start=1)]
    pantallas.append(pantalla_hijas(rep, ancho))       # el último paso: la ligasa
    pantallas += [pantalla_transcripcion(tx, e, i, len(tx.estados), ancho)
                  for i, e in enumerate(tx.estados, start=1)]
    pantallas += [pantalla_traduccion(tr, e, i, len(tr.estados), ancho)
                  for i, e in enumerate(tr.estados, start=1)]
    pantallas.append(pantalla_resumen(rep, tx, tr, ancho))
    return pantallas
