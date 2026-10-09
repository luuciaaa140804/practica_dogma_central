#!/usr/bin/env python3
"""
Simulador del dogma central de la biología molecular.

    python3 simulador.py                  menú
    python3 simulador.py ATGGCCAAGTAA...  con una secuencia (hebra codificante, 5'→3')
    python3 simulador.py -f gen.fasta     desde un archivo FASTA
    python3 simulador.py -a 60            con una secuencia aleatoria de 60 pb
    python3 simulador.py --sin-pausa      todas las pantallas seguidas
    python3 simulador.py --imagen salida/ además guarda imágenes PNG (necesita matplotlib)
"""

import argparse
import os
import shutil
import sys

from dogma import adn, dibujo
from dogma.replicacion import replicar
from dogma.traduccion import traducir
from dogma.transcripcion import transcribir

EJEMPLO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ejemplos",
                       "gen_ejemplo.fasta")


def simular(secuencia, descripcion, ancho, carpeta_imagenes=None):
    """Replicación → transcripción → traducción. Devuelve la lista de pantallas."""
    rep = replicar(secuencia)
    hija1 = rep.hijas()[0]
    tx = transcribir(hija1[0])        # se transcribe una de las moléculas recién replicadas
    tr = traducir(tx.arnm)
    if carpeta_imagenes:
        from dogma import imagen      # matplotlib solo se necesita si se piden imágenes
        rutas = imagen.guardar_imagenes(carpeta_imagenes, rep, tx, tr)
        print(f"Guardadas {len(rutas)} imágenes en {carpeta_imagenes}/")
    return dibujo.todas(secuencia, descripcion, rep, tx, tr, ancho)


def mostrar(pantallas, pausa):
    for i, pantalla in enumerate(pantallas):
        if pausa:
            print("\033[2J\033[H", end="")      # borra la terminal
        print("\n".join(pantalla))
        print()
        if pausa and i < len(pantallas) - 1:
            orden = input(dibujo.pintar(" [Enter] siguiente · [q] salir > ", "tenue"))
            if orden.strip().lower() == "q":
                break


def menu():
    print(dibujo.pintar(" SIMULADOR DEL DOGMA CENTRAL DE LA BIOLOGÍA MOLECULAR ", "titulo"))
    print("\n  1) Gen de ejemplo\n  2) Escribir una secuencia\n  3) Leer un archivo FASTA\n"
          "  4) Secuencia aleatoria\n")
    opcion = input(" Opción [1]: ").strip() or "1"
    if opcion == "2":
        return input(" Hebra codificante (5'→3'): "), "secuencia escrita por el usuario"
    if opcion == "3":
        ruta = input(" Ruta del archivo: ").strip()
        return adn.leer_fasta(ruta), f"archivo {ruta}"
    if opcion == "4":
        longitud = input(" Longitud en pb [60]: ").strip() or "60"
        if not longitud.isdigit():
            raise ValueError("la longitud debe ser un número entero")
        return aleatoria(int(longitud))
    return adn.leer_fasta(EJEMPLO), "gen de ejemplo (péptido líder trp de E. coli)"


def aleatoria(longitud):
    if longitud < 15:
        raise ValueError("la secuencia aleatoria debe tener al menos 15 pb")
    return adn.gen_aleatorio(longitud), f"secuencia aleatoria de {longitud} pb"


def main():
    parser = argparse.ArgumentParser(description="Simulador del dogma central de la biología "
                                                 "molecular: replicación, transcripción y "
                                                 "traducción.")
    parser.add_argument("secuencia", nargs="?", help="hebra codificante del ADN, 5'→3'")
    parser.add_argument("-f", "--fasta", help="lee la secuencia de un archivo FASTA")
    parser.add_argument("-a", "--aleatoria", type=int, metavar="N",
                        help="usa una secuencia aleatoria de N pb")
    parser.add_argument("--sin-pausa", action="store_true",
                        help="muestra todas las pantallas seguidas")
    parser.add_argument("--imagen", metavar="CARPETA",
                        help="guarda además imágenes PNG de la simulación (necesita matplotlib)")
    parser.add_argument("--sin-color", action="store_true", help="desactiva los colores")
    args = parser.parse_args()

    if os.name == "nt":
        os.system("")                  # activa los colores ANSI en la consola de Windows
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except AttributeError:
        pass
    terminal = sys.stdin.isatty() and sys.stdout.isatty()
    dibujo.USAR_COLOR = sys.stdout.isatty() and not args.sin_color
    ancho = max(70, min(shutil.get_terminal_size((100, 30)).columns - 1, 140))

    try:
        if args.secuencia:
            texto, descripcion = args.secuencia, "secuencia de la línea de órdenes"
        elif args.fasta:
            texto, descripcion = adn.leer_fasta(args.fasta), f"archivo {args.fasta}"
        elif args.aleatoria:
            texto, descripcion = aleatoria(args.aleatoria)
        elif terminal:
            texto, descripcion = menu()
        else:
            texto, descripcion = adn.leer_fasta(EJEMPLO), "gen de ejemplo"
        secuencia = adn.limpiar(texto)
    except OSError as error:
        sys.exit(f"Error: no se puede leer el archivo ({error.strerror}).")
    except ValueError as error:
        sys.exit(f"Error: {error}.")
    except (KeyboardInterrupt, EOFError):
        sys.exit()

    try:
        mostrar(simular(secuencia, descripcion, ancho, args.imagen), pausa=terminal and not args.sin_pausa)
    except (KeyboardInterrupt, EOFError):
        print()


if __name__ == "__main__":
    main()
