import pandas as pd
import re
import numpy as np
from pyscf import gto, dft, scf

PYSCF_AVAILABLE = True



class Orca:
    def __init__(self, ruta_salida):
        self.ruta = ruta_salida
        try:
            with open(ruta_salida, 'r', encoding='utf-8', errors='ignore') as f:
                self.contenido = f.read()
        except FileNotFoundError:
            raise FileNotFoundError(f"No se encontro el archivo de salida en: {ruta_salida}")

    #proxima mejora, ignorar lineas en blanco y comentarios al parsear xyz
    @staticmethod
    def generar_entrada(contenido_xyz, tipo_calculo, metodo, base, palabras_clave,
                        calc_nmr=False):
        palabras_base = f"! {metodo} {base} {palabras_clave}"

        if tipo_calculo == "Optimizacion de Geometria":
            palabras_calculo = "OPT"
        elif tipo_calculo == "Frecuencias Vibracionales (IR)":
            palabras_calculo = "OPT FREQ"
        else:
            palabras_calculo = ""

        if calc_nmr:
            palabras_calculo += " NMR"

        encabezado = f"{palabras_base} {palabras_calculo}\n"
        lineas = contenido_xyz.strip().split('\n')

        try:
            num_atomos = int(lineas[0].strip())
            if len(lineas[1].strip().split()) > 1 and lineas[1].strip().split()[0].isalpha():
                lineas_coords = lineas[1:1 + num_atomos]
            else:
                lineas_coords = lineas[2:2 + num_atomos]
        except (ValueError, IndexError):
            lineas_coords = lineas[2:]

        coords_str = "\n".join(lineas_coords)
        bloque_xyz = f"* xyz 0 1\n{coords_str}\n*\n"
        return encabezado + bloque_xyz

    def verificar_convergencia(self):
        return "THE OPTIMIZATION HAS CONVERGED" in self.contenido

    def extraer_energia_final(self):
        coincidencias = re.findall(r'FINAL SINGLE POINT ENERGY\s+([-\d.]+)', self.contenido)
        if coincidencias:
            return float(coincidencias[-1])
        return None

    def extraer_geometria_optimizada(self):
        patron = r'CARTESIAN COORDINATES \(ANGSTROEM\)\s*\n\s*-+\s*\n((?:\s*\S+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s*\n)+)'
        coincidencias = list(re.finditer(patron, self.contenido))
        if not coincidencias:
            return None

        bloque_coords = coincidencias[-1].group(1).strip()
        lineas_coords = [linea.strip() for linea in bloque_coords.split('\n') if linea.strip()]
        if not lineas_coords:
            return None

        num_atomos = len(lineas_coords)
        bloque_xyz = f"{num_atomos}\nGeometria Optimizada extraida de {self.ruta}\n"
        for linea in lineas_coords:
            partes = linea.split()
            if len(partes) >= 4:
                bloque_xyz += f"{partes[0]:<2} " + " ".join(f"{float(coord):>12.6f}" for coord in partes[1:4]) + "\n"
        return bloque_xyz

    def extraer_espectro_ir(self, factor_escalamiento=1.0):
        patron = r'IR SPECTRUM\s*\n-+\n(?:.|\n)*?-+\n((?:.|\n)*?)(?=\n\s*\*|\n\s*-{2,}\n[A-Z]|\Z)'
        coincidencia = re.search(patron, self.contenido)

        if not coincidencia:
            return pd.DataFrame()

        datos = []
        bloque_datos = coincidencia.group(1).strip()

        for linea in bloque_datos.split('\n'):
            partes = linea.split()
            if len(partes) > 3 and partes[0].endswith(':'):
                try:
                    freq = float(partes[1])
                    intensidad = float(partes[3])
                    if freq > 10.0:
                        datos.append({"Frequency": freq * factor_escalamiento, "Intensity": intensidad})
                except (ValueError, IndexError):
                    continue

        return pd.DataFrame(datos)

    def extraer_componentes_energia(self):
        patrones = {
            "Repulsion Nuclear": r'Nuclear Repulsion\s+:\s*([-\d.]+)',
            "Energia Electronica": r'Electronic Energy\s+:\s*([-\d.]+)',
            "Energia Un Electron": r'One Electron Energy\s+:\s*([-\d.]+)',
            "Energia Dos Electrones": r'Two Electron Energy\s+:\s*([-\d.]+)',
        }
        energias = {}
        for nombre, patron in patrones.items():
            coincidencias = re.findall(patron, self.contenido)
            if coincidencias:
                energias[nombre] = [float(coincidencias[-1])]

        return pd.DataFrame.from_dict(energias, orient='index', columns=['Energia (Hartree)']) if energias else None

    def extraer_cargas_atomicas(self):
        datos_cargas = {}
        for tipo in ['MULLIKEN', 'LOEWDIN']:
            patron = re.compile(rf'{tipo} ATOMIC CHARGES\s*\n-+\n((?:.|\n)*?)(?=\n\n|\Z)')

            coincidencias = list(re.finditer(patron, self.contenido))
            if coincidencias:
                coincidencia_final = coincidencias[-1]
                cargas = []
                for linea in coincidencia_final.group(1).strip().split('\n'):
                    partes = linea.split()
                    if len(partes) == 4 and partes[2] == ':':
                        cargas.append({"Atomo": f"{partes[0]} {partes[1]}", "Carga": float(partes[3])})
                if cargas:
                    datos_cargas[tipo.capitalize()] = pd.DataFrame(cargas)

        return datos_cargas if datos_cargas else None

    def extraer_energias_orbitales(self):
        patron = re.compile(r'ORBITAL ENERGIES\s*\n-+\n((?:.|\n)*?)(?=\n\n|\Z|\*Only the first)')

        coincidencias = list(re.finditer(patron, self.contenido))
        if not coincidencias:
            return None

        coincidencia_final = coincidencias[-1]

        orbitales = []
        for linea in coincidencia_final.group(1).strip().split('\n')[2:]:
            partes = linea.split()
            if len(partes) == 4:
                orbitales.append({
                    "Numero": int(partes[0]),
                    "Ocupacion": float(partes[1]),
                    "Energia (Eh)": float(partes[2]),
                    "Energia (eV)": float(partes[3])
                })

        return pd.DataFrame(orbitales) if orbitales else None

    def extraer_cargas_orbitales_reducidas(self):
        datos_cargas = {}
        for tipo in ['MULLIKEN', 'LOEWDIN']:
            patron = re.compile(
                rf'{tipo} REDUCED ORBITAL CHARGES\s*\n-+\n((?:.|\n)*?)(?=\n\n|\Z|\s*\*+\n|\s*-{{2,}}\n[A-Z])')

            coincidencias = list(re.finditer(patron, self.contenido))
            if coincidencias:
                coincidencia_final = coincidencias[-1]
                cargas_orbitales = []
                atomo_actual = ""

                bloque_texto = coincidencia_final.group(1).strip()
                for linea in bloque_texto.split('\n'):
                    linea_limpia = linea.strip()
                    if not linea_limpia:
                        continue

                    partes_iniciales = linea.split()
                    if partes_iniciales and partes_iniciales[0].isdigit():
                        if len(partes_iniciales) > 1:
                            atomo_actual = f"{partes_iniciales[0]} {partes_iniciales[1]}"

                    pares_orbital_carga = re.findall(r'([a-zA-Z0-9]+)\s*:\s*([\d.-]+)', linea)

                    for orbital, carga in pares_orbital_carga:
                        if orbital in ['s', 'p', 'd', 'f'] and atomo_actual:
                            try:
                                cargas_orbitales.append({
                                    "Atomo": atomo_actual,
                                    "Orbital": orbital,
                                    "Carga": float(carga)
                                })
                            except ValueError:
                                continue

                if cargas_orbitales:
                    datos_cargas[tipo.capitalize()] = pd.DataFrame(cargas_orbitales)

        return datos_cargas if datos_cargas else None

    def extraer_datos_nmr(self):
        patron_bloque = re.compile(
            r'CHEMICAL SHIELDING SUMMARY \(ppm\)\s*\n-+\n\n((?:.|\n)*?)(?=\n\n\s*NMR shielding tensor|\Z|\n\s*-{2,}\n)')
        coincidencia = re.search(patron_bloque, self.contenido)

        if not coincidencia:
            return None

        bloque_texto = coincidencia.group(1).strip()
        lineas = bloque_texto.split('\n')

        if len(lineas) <= 2:
            return None

        datos_nmr = []
        for linea in lineas[2:]:
            partes = linea.split()
            if len(partes) == 4:
                try:
                    datos_nmr.append({
                        "Nucleo": int(partes[0]),
                        "Elemento": partes[1],
                        "Isotropico (ppm)": float(partes[2]),
                        "Anisotropia (ppm)": float(partes[3])
                    })
                except ValueError:
                    continue

        if datos_nmr:
            return pd.DataFrame(datos_nmr)

        return None


class PySCFCalculator:

    @staticmethod
    def calcular_susceptibilidad(xyz_content, metodo='b3lyp', base='def2svp'):
        try:
            lineas = [l.strip() for l in xyz_content.strip().split('\n') if l.strip()]

            try:
                num_atomos = int(lineas[0])
            except (ValueError, IndexError):
                return {"error": "Formato XYZ invalido: primera linea debe ser numero de atomos"}

            if len(lineas) < num_atomos + 2:
                return {"error": "Formato XYZ invalido: faltan lineas de coordenadas"}

            lineas_coords = lineas[2:2 + num_atomos]

            atom_str = ""
            for idx, linea in enumerate(lineas_coords):
                partes = linea.split()
                if len(partes) < 4:
                    return {"error": f"Linea {idx + 3} invalida: {linea}"}
                try:
                    x, y, z = float(partes[1]), float(partes[2]), float(partes[3])
                    atom_str += f"{partes[0]} {x} {y} {z}; "
                except ValueError:
                    return {"error": f"Coordenadas invalidas en linea {idx + 3}: {linea}"}

            base_map = {
                'def2-svp': 'def2svp',
                'def2-tzvp': 'def2tzvp',
                '6-31+g(d,p)': '6-31+g*',
                '6-311++g(d,p)': '6-311++g**',
                'cc-pvdz': 'ccpvdz'
            }
            base_pyscf = base_map.get(base.lower(), base.lower())

            mol = gto.M(
                atom=atom_str,
                basis=base_pyscf,
                unit='Angstrom'
            )

            mf = dft.RKS(mol)
            mf.xc = metodo.lower()

            energia = mf.kernel()

            if not mf.converged:
                return {"error": "SCF no convergio en PySCF"}


            coords = mol.atom_coords()
            charges = mol.atom_charges()

            chi_dia = 0.0

            total_mass = sum(charges)
            com = np.sum(coords * charges[:, np.newaxis], axis=0) / total_mass if total_mass > 0 else np.zeros(3)
            chi_tensor = np.zeros((3, 3))

            for i, (coord, Z) in enumerate(zip(coords, charges)):
                r = coord - com
                r2 = np.dot(r, r)

                for j in range(3):
                    for k in range(3):
                        if j == k:
                            chi_tensor[j, k] -= Z * (r2 - r[j] ** 2) / 6.0
                        else:
                            chi_tensor[j, k] -= Z * r[j] * r[k] / 6.0

            chi_iso = np.trace(chi_tensor) / 3.0

            chi_cgs = chi_iso * 0.78910

            tipo_magnetismo = "Diamagnetico" if chi_cgs < 0 else "Paramagnetico"

            return {
                "tensor": chi_tensor.tolist(),
                "isotropico_au": float(chi_iso),
                "isotropico_cgs": float(chi_cgs),
                "tipo": tipo_magnetismo,
                "energia_scf": float(energia),
                "converged": True,
                "metodo_calculo": "Aproximacion de Pascal (diamagnetica)",
                "nota": "Calculo aproximado basado en geometria molecular. Para resultados precisos usar ORCA con palabras clave NMR."
            }

        except ImportError as e:
            return {
                "error": f"PySCF no esta correctamente instalado: {str(e)}\nIntenta: pip install --upgrade pyscf"
            }
        except Exception as e:
            import traceback
            return {"error": f"Error en calculo: {str(e)}\n\nDetalle:\n{traceback.format_exc()}"}