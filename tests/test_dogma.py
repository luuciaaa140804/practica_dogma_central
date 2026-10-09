"""Pruebas del simulador. Se ejecutan con:  python3 -m unittest -v"""

import random
import subprocess
import sys
import unittest
from pathlib import Path

from dogma.adn import (CODIGO_GENETICO, anticodon, arn_complementario, complementaria,
                       gen_aleatorio, leer_fasta, limpiar)
from dogma.replicacion import CEBADOR, LIDER, OKAZAKI, replicar
from dogma.traduccion import traducir
from dogma.transcripcion import transcribir

RAIZ = Path(__file__).resolve().parent.parent
EJEMPLO = leer_fasta(RAIZ / "ejemplos" / "gen_ejemplo.fasta")


class TestADN(unittest.TestCase):
    def test_complementariedad(self):
        self.assertEqual(complementaria("ATGC"), "TACG")
        self.assertEqual(arn_complementario("TACG"), "AUGC")   # el ARN lleva U, no T
        self.assertEqual(anticodon("AUG"), "UAC")             # 5'-AUG-3' / 3'-UAC-5'

    def test_codigo_genetico(self):
        self.assertEqual(len(CODIGO_GENETICO), 64)
        self.assertEqual([c for c, aa in CODIGO_GENETICO.items() if aa == "*"],
                         ["UAA", "UAG", "UGA"])
        for codon, aa in {"AUG": "M", "UGG": "W", "UUU": "F", "GGC": "G", "AAA": "K"}.items():
            self.assertEqual(CODIGO_GENETICO[codon], aa)

    def test_limpiar(self):
        self.assertEqual(limpiar(" atg gcc 12 aag taa "), "ATGGCCAAGTAA")
        with self.assertRaises(ValueError):
            limpiar("ATGGCCAAGUAA")          # U no es una base del ADN
        with self.assertRaises(ValueError):
            limpiar("ATG")                   # demasiado corta


class TestReplicacion(unittest.TestCase):
    def test_muchas_secuencias(self):
        rng = random.Random(1)
        for _ in range(300):
            adn = "".join(rng.choice("ACGT") for _ in range(rng.randint(12, 150)))
            rep = replicar(adn, avance=rng.randint(1, 12), cebador=rng.randint(1, 4))
            self.assertTrue(rep.es_correcta())
            (sup, nueva_sup), (nueva_inf, inf) = rep.hijas()
            self.assertEqual((sup, inf), (adn, complementaria(adn)))   # semiconservativa
            self.assertEqual((nueva_sup, nueva_inf), (complementaria(adn), adn))

    def test_lider_continua_y_rezagada_en_fragmentos(self):
        rep = replicar(EJEMPLO)
        for estado in rep.estados[:-2]:        # pasos en los que avanza la horquilla
            lider = estado.lider[:estado.horquilla]
            # La cadena líder llega hasta la horquilla sin huecos y solo tiene un cebador
            self.assertNotIn(None, lider)
            self.assertEqual([nt.tipo for nt in lider[:rep.cebador]], [CEBADOR] * rep.cebador)
            self.assertTrue(all(nt.tipo == LIDER for nt in lider[rep.cebador:]))
        antes_de_madurar = rep.estados[-3]
        tipos = [nt.tipo for nt in antes_de_madurar.rezagada]
        self.assertEqual(tipos.count(CEBADOR), len(rep.fragmentos) * rep.cebador)
        self.assertEqual(len(rep.fragmentos), 3)
        self.assertTrue(all(t in (CEBADOR, OKAZAKI) for t in tipos))

    def test_cebadores_complementarios(self):
        rep = replicar(EJEMPLO)
        for nt, base in zip(rep.estados[-3].rezagada, rep.superior):
            if nt.tipo == CEBADOR:
                self.assertEqual(nt.base, arn_complementario(base))

    def test_maduracion(self):
        rep = replicar(EJEMPLO)
        con_mellas, sellada = rep.estados[-2], rep.estados[-1]
        self.assertEqual(len(con_mellas.mellas), len(rep.fragmentos) + 1)
        self.assertEqual(sellada.mellas, [])
        self.assertNotIn(CEBADOR, [nt.tipo for nt in sellada.rezagada + sellada.lider])


class TestTranscripcionTraduccion(unittest.TestCase):
    def test_transcripcion(self):
        tx = transcribir(EJEMPLO)
        self.assertEqual(tx.arnm, EJEMPLO.replace("T", "U"))   # codificante con U
        self.assertEqual(tx.arnm, arn_complementario(tx.molde))
        self.assertEqual(tx.estados[-1].posicion, len(EJEMPLO))

    def test_traduccion_del_ejemplo(self):
        tr = traducir(EJEMPLO.replace("T", "U"))
        self.assertEqual(tr.proteina, "MKAIFVLKGWWRTS")     # péptido líder trp de E. coli
        self.assertEqual(tr.stop, "UGA")
        self.assertEqual(tr.tres_letras()[0], "fMet")
        self.assertEqual([e.fase for e in tr.estados], ["inicio", "elongacion", "terminacion"])

    def test_casos_especiales(self):
        self.assertEqual(traducir("GGCCCUUU").proteina, "")          # sin AUG
        sin_stop = traducir("AUGGCCGCC")
        self.assertEqual((sin_stop.proteina, sin_stop.stop), ("MAA", None))
        self.assertEqual(traducir("CCAUGUAA").proteina, "M")

    def test_secuencia_aleatoria(self):
        adn = gen_aleatorio(90, semilla=3)
        self.assertEqual(len(adn), 90)
        self.assertEqual(traducir(adn.replace("T", "U")).proteina[0], "M")


class TestPrograma(unittest.TestCase):
    def test_simulacion_completa(self):
        salida = subprocess.run([sys.executable, str(RAIZ / "simulador.py"), "--sin-pausa"],
                                capture_output=True, text=True, encoding="utf-8")
        self.assertEqual(salida.returncode, 0)
        self.assertIn("MKAIFVLKGWWRTS", salida.stdout)
        self.assertNotIn("✘", salida.stdout)

    def test_secuencia_no_valida(self):
        salida = subprocess.run([sys.executable, str(RAIZ / "simulador.py"), "ATGXXXGGCCAA"],
                                capture_output=True, text=True, encoding="utf-8")
        self.assertNotEqual(salida.returncode, 0)
        self.assertIn("no válidos", salida.stderr)




class TestImagen(unittest.TestCase):
    def test_guarda_imagenes(self):
        try:
            import matplotlib  # noqa: F401
        except ImportError:
            self.skipTest("matplotlib no está instalado")
        import tempfile
        from dogma import imagen
        rep = replicar(EJEMPLO)
        tx = transcribir(EJEMPLO)
        tr = traducir(tx.arnm)
        with tempfile.TemporaryDirectory() as carpeta:
            rutas = imagen.guardar_imagenes(carpeta, rep, tx, tr)
            self.assertEqual(len(rutas), len(rep.estados) + 1)
            self.assertTrue(all(Path(r).stat().st_size > 0 for r in rutas))


if __name__ == "__main__":
    unittest.main()
