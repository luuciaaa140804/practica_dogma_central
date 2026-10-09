"""
Salida en imagen (PNG) de la simulación. Es opcional: necesita matplotlib.

    python3 simulador.py --imagen salida/      guarda una imagen por paso de la replicación
                                               y una con el flujo completo ADN → proteína
"""

import os

from .adn import AMINOACIDOS, CODIGO_GENETICO, CODONES_STOP, anticodon, complementaria
from .replicacion import CEBADOR

COLOR = {"molde": "#9e9e9e", "lider": "#2e7d32", "okazaki": ["#1565c0", "#ef6c00"],
         "cebador": "#8e24aa", "relleno": "#00897b"}


def _matplotlib():
    try:
        import matplotlib
        matplotlib.use("Agg")                 # sin ventana: solo se guardan ficheros
        import matplotlib.pyplot as plt
        return plt
    except ImportError:
        raise SystemExit("Error: para --imagen hace falta matplotlib (pip install matplotlib).")


def _tramos(hebra):
    """Agrupa posiciones consecutivas del mismo tipo y fragmento → [(ini, fin, tipo, frag)]."""
    tramos, actual = [], None
    for i, nt in enumerate(list(hebra) + [None]):
        clave = None if nt is None else (nt.tipo, nt.fragmento)
        if actual and actual[2] == clave:
            actual[1] = i + 1
            continue
        if actual and actual[2] is not None:
            tramos.append((actual[0], actual[1], *actual[2]))
        actual = [i, i + 1, clave]
    return tramos


def dibujar_replicacion(rep, estado, paso, total, ruta):
    plt = _matplotlib()
    n, pos = len(rep.superior), estado.horquilla
    y_sup, y_rez, y_dup, y_lid, y_inf = 3.0, 2.2, 1.5, 0.8, 0.0
    fig, ax = plt.subplots(figsize=(13, 4.4))

    def caja(x, y, texto, color):
        ax.text(x, y, texto, ha="center", va="center", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.25", fc=color, ec="none"))

    # Moldes: separados detrás de la horquilla, unidos por delante
    ax.plot([0, pos], [y_sup, y_sup], color=COLOR["molde"], lw=2)
    ax.plot([0, pos], [y_inf, y_inf], color=COLOR["molde"], lw=2)
    if pos < n:
        ax.plot([pos, pos + 1, n], [y_sup, y_dup + .1, y_dup + .1], color=COLOR["molde"], lw=2)
        ax.plot([pos, pos + 1, n], [y_inf, y_dup - .1, y_dup - .1], color=COLOR["molde"], lw=2)
        caja(pos + .5, y_dup, "Helicasa", "#fff59d")
        if pos < n - 8:
            caja(pos + 8, y_dup + .45, "Girasa", "#ffcc80")

    # Líder: flecha continua hacia la horquilla (el cebador, en morado)
    for a, b, tipo, _ in _tramos(estado.lider):
        if tipo == CEBADOR:
            ax.plot([a, b], [y_lid] * 2, color=COLOR["cebador"], lw=7, solid_capstyle="butt")
        else:
            ax.annotate("", xy=(b, y_lid), xytext=(a, y_lid),
                        arrowprops=dict(arrowstyle="-|>", color=COLOR["lider"], lw=2.5))
    # Rezagada: cada fragmento de Okazaki crece alejándose de la horquilla (punta a la izquierda)
    for a, b, tipo, k in _tramos(estado.rezagada):
        if tipo == CEBADOR:
            ax.plot([a, b], [y_rez] * 2, color=COLOR["cebador"], lw=7, solid_capstyle="butt")
        else:
            color = COLOR["okazaki"][k % 2] if tipo != "relleno" else COLOR["relleno"]
            ax.annotate("", xy=(a, y_rez), xytext=(b, y_rez),
                        arrowprops=dict(arrowstyle="-|>", color=color, lw=2.5))
    for k, (ini, fin) in enumerate(rep.fragmentos, start=1):
        if estado.rezagada[ini] is not None:
            ax.text((ini + fin) / 2, y_rez - .4, f"Okazaki {k}", ha="center", fontsize=7,
                    color=COLOR["okazaki"][k % 2])
    for hebra, x in estado.mellas:
        ax.plot([x, x], [y_rez - .15 if hebra == "rezagada" else y_lid - .15] * 2,
                marker="v", color="#e53935", ms=5)

    if estado.mellas == [] and pos < n:
        caja(max(pos - 2, 1), y_lid - .45, "ADN pol III", "#c8e6c9")
    for y, nombre, ext in ((y_sup, "molde sup.", ("5'", "3'")), (y_rez, "rezagada", ("3'", "5'")),
                           (y_lid, "líder", ("5'", "3'")), (y_inf, "molde inf.", ("3'", "5'"))):
        ax.text(-3, y, nombre, ha="right", va="center", fontsize=9)
        ax.text(-1, y, ext[0], ha="right", va="center", fontsize=8)
        ax.text(n + 1, y, ext[1], ha="left", va="center", fontsize=8)

    ax.set_xlim(-12, n + 4)
    ax.set_ylim(-1, 3.8)
    ax.axis("off")
    ax.set_title(f"Replicación · paso {paso} de {total} · {estado.titulo}", fontsize=11)
    ax.text(0.5, -.02, "morado: cebador de ARN · verde: líder · azul/naranja: fragmentos de "
            "Okazaki · ▼ mella", transform=ax.transAxes, ha="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)


def dibujar_flujo(tx, tr, ruta):
    """ADN → ARNm → anticodón → aminoácido, codón a codón."""
    plt = _matplotlib()
    if not tr.codones:
        return False
    cods = tr.codones
    fig, ax = plt.subplots(figsize=(max(8, 0.75 * len(cods) + 3), 3.6))
    filas = [("ADN codificante 5'→3'", 4, lambda c, j: c.replace("U", "T"), "#616161"),
             ("ADN molde 3'→5'", 3, lambda c, j: complementaria(c.replace("U", "T")), "#9e9e9e"),
             ("ARNm 5'→3'", 2, lambda c, j: c, "#8e24aa"),
             ("anticodón ARNt 3'→5'", 1, lambda c, j: "" if c in CODONES_STOP else anticodon(c),
              "#1565c0"),
             ("aminoácido", 0, lambda c, j: "Stop" if c in CODONES_STOP else
              ("fMet" if j == 0 else AMINOACIDOS[CODIGO_GENETICO[c]][0]), "#2e7d32")]
    for nombre, y, f, color in filas:
        ax.text(-0.7, y, nombre, ha="right", va="center", fontsize=9)
        for j, c in enumerate(cods):
            ax.text(j, y, f(c, j), ha="center", va="center", fontsize=9, color=color,
                    family="monospace", fontweight="bold" if y == 0 else "normal")
    ax.set_xlim(-6.5, len(cods))
    ax.set_ylim(-.8, 4.8)
    ax.axis("off")
    ax.set_title("Flujo de la información: ADN → ARN → proteína", fontsize=11)
    ax.text((len(cods) - 1) / 2, -.75, f"Proteína: {tr.proteina}", ha="center", fontsize=10)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    return True


def guardar_imagenes(carpeta, rep, tx, tr):
    """Guarda una imagen por paso de replicación y la del flujo completo. Devuelve las rutas."""
    os.makedirs(carpeta, exist_ok=True)
    rutas = []
    for i, estado in enumerate(rep.estados, start=1):
        ruta = os.path.join(carpeta, f"replicacion_{i:02d}.png")
        dibujar_replicacion(rep, estado, i, len(rep.estados), ruta)
        rutas.append(ruta)
    ruta = os.path.join(carpeta, "flujo_completo.png")
    if dibujar_flujo(tx, tr, ruta):
        rutas.append(ruta)
    return rutas
