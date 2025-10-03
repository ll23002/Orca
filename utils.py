import pandas as pd
import re


class OrcaParser:
    """
    Una clase para generar entradas y analizar archivos de salida de ORCA.

    Esta clase encapsula un conjunto de herramientas para la química computacional con ORCA.
    Lee un archivo de salida una sola vez y permite ejecutar múltiples funciones de
    extracción de datos sobre su contenido, lo cual es mucho más eficiente.

    Atributos:
        ruta (str): La ruta al archivo de salida de ORCA.
        contenido (str): El contenido completo del archivo de salida leído.
    """

    def __init__(self, ruta_salida):
        """
        Inicializa el analizador leyendo el archivo de salida de ORCA.

        Args:
            ruta_salida (str): La ruta al archivo .out de ORCA.

        Raises:
            FileNotFoundError: Si el archivo no se encuentra en la ruta especificada.
        """
        self.ruta = ruta_salida
        try:
            with open(ruta_salida, 'r', encoding='utf-8', errors='ignore') as f:
                self.contenido = f.read()
        except FileNotFoundError:
            raise FileNotFoundError(f"No se encontró el archivo de salida en: {ruta_salida}")

    @staticmethod
    def generar_entrada_orca(contenido_xyz, tipo_calculo, metodo, base, palabras_clave):
        """
        Genera una cadena de texto para un archivo de entrada de ORCA de forma robusta.
        Esta función es estática porque no depende de un archivo de salida.
        """
        palabras_base = f"! {metodo} {base} {palabras_clave}"

        if tipo_calculo == "Optimización de Geometría":
            palabras_calculo = "OPT"
        elif tipo_calculo == "Frecuencias Vibracionales (IR)":
            palabras_calculo = "OPT FREQ"
        else:
            palabras_calculo = ""

        encabezado = f"{palabras_base} {palabras_calculo}\n"
        lineas = contenido_xyz.strip().split('\n')

        # --- MEJORA: Lectura robusta de coordenadas ---
        try:
            num_atomos = int(lineas[0].strip())
            # Detecta si la segunda línea ya es una coordenada o un comentario
            if len(lineas[1].strip().split()) > 1 and lineas[1].strip().split()[0].isalpha():
                lineas_coords = lineas[1:1 + num_atomos]
            else:
                lineas_coords = lineas[2:2 + num_atomos]
        except (ValueError, IndexError):
            # Método de respaldo en caso de formato inesperado
            lineas_coords = lineas[2:]

        coords_str = "\n".join(lineas_coords)
        bloque_xyz = f"* xyz 0 1\n{coords_str}\n*\n"
        return encabezado + bloque_xyz

    def verificar_convergencia_optimizacion(self):
        """Verifica si la optimización de geometría ha convergido."""
        return "THE OPTIMIZATION HAS CONVERGED" in self.contenido

    def extraer_energia_final(self):
        """Extrae la energía final de punto único (FINAL SINGLE POINT ENERGY)."""
        coincidencias = re.findall(r'FINAL SINGLE POINT ENERGY\s+([-\d.]+)', self.contenido)
        if coincidencias:
            return float(coincidencias[-1])
        return None

    def extraer_geometria_optimizada(self):
        """Extrae la última geometría del archivo (la optimizada) en formato XYZ."""
        patron = r'CARTESIAN COORDINATES \(ANGSTROEM\)\s*\n\s*-+\s*\n((?:\s*\S+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s*\n)+)'
        coincidencias = list(re.finditer(patron, self.contenido))
        if not coincidencias:
            return None

        bloque_coords = coincidencias[-1].group(1).strip()
        lineas_coords = [linea.strip() for linea in bloque_coords.split('\n') if linea.strip()]
        if not lineas_coords:
            return None

        num_atomos = len(lineas_coords)
        bloque_xyz = f"{num_atomos}\nGeometría Optimizada extraída de {self.ruta}\n"
        for linea in lineas_coords:
            partes = linea.split()
            if len(partes) >= 4:
                bloque_xyz += f"{partes[0]:<2} {float(partes[1]):>12.6f} {float(partes[2]):>12.6f} {float(partes[3]):>12.6f}\n"
        return bloque_xyz

    def extraer_espectro_ir(self, factor_escalamiento=1.0):
        """Extrae las frecuencias vibracionales y las intensidades del espectro IR."""
        patron_ir_bloque = r'IR SPECTRUM\s*\n-+\n(?:.|\n)*?-+\n((?:.|\n)*?)(?=\n\s*\*|\n\s*-{2,}\n[A-Z]|\Z)'
        coincidencia = re.search(patron_ir_bloque, self.contenido)

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
                    # Filtro para ignorar modos de traslación/rotación
                    if freq > 10.0:
                        datos.append({"Frequency": freq * factor_escalamiento, "Intensity": intensidad})
                except (ValueError, IndexError):
                    continue

        return pd.DataFrame(datos)

    def extraer_componentes_energia(self):
        """Extrae los componentes finales de la energía (nuclear, electrónica, etc.)."""
        patrones = {
            "Repulsión Nuclear": r'Nuclear Repulsion\s+:\s*([-\d.]+)',
            "Energía Electrónica": r'Electronic Energy\s+:\s*([-\d.]+)',
            "Energía Un Electrón": r'One Electron Energy\s+:\s*([-\d.]+)',
            "Energía Dos Electrones": r'Two Electron Energy\s+:\s*([-\d.]+)',
        }
        energias = {}
        for nombre, patron in patrones.items():
            coincidencias = re.findall(patron, self.contenido)
            if coincidencias:
                energias[nombre] = [float(coincidencias[-1])]

        return pd.DataFrame.from_dict(energias, orient='index', columns=['Energía (Hartree)']) if energias else None

    def extraer_cargas_atomicas(self):
        """Extrae las cargas atómicas finales de Mulliken y Loewdin."""
        datos_cargas = {}
        for tipo in ['MULLIKEN', 'LOEWDIN']:
            patron = re.compile(rf'{tipo} ATOMIC CHARGES\s*\n-+\n((?:.|\n)*?)(?=\n\n|\Z)')

            # --- CORRECCIÓN: Usar finditer y tomar la última coincidencia ---
            coincidencias = list(re.finditer(patron, self.contenido))
            if coincidencias:
                coincidencia_final = coincidencias[-1]
                cargas = []
                for linea in coincidencia_final.group(1).strip().split('\n'):
                    partes = linea.split()
                    if len(partes) == 4 and partes[2] == ':':
                        cargas.append({"Átomo": f"{partes[0]} {partes[1]}", "Carga": float(partes[3])})
                if cargas:
                    datos_cargas[tipo.capitalize()] = pd.DataFrame(cargas)

        return datos_cargas if datos_cargas else None

    def extraer_energias_orbitales(self):
        """Extrae las energías de los orbitales moleculares finales."""
        patron = re.compile(r'ORBITAL ENERGIES\s*\n-+\n((?:.|\n)*?)(?=\n\n|\Z|\*Only the first)')

        # --- CORRECCIÓN: Usar finditer y tomar la última coincidencia ---
        coincidencias = list(re.finditer(patron, self.contenido))
        if not coincidencias:
            return None

        coincidencia_final = coincidencias[-1]

        orbitales = []
        # Omitir los encabezados de la tabla (índice [2:])
        for linea in coincidencia_final.group(1).strip().split('\n')[2:]:
            partes = linea.split()
            if len(partes) == 4:
                orbitales.append({
                    "Número": int(partes[0]),
                    "Ocupación": float(partes[1]),
                    "Energía (Eh)": float(partes[2]),
                    "Energía (eV)": float(partes[3])
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
                                    "Átomo": atomo_actual,
                                    "Orbital": orbital,
                                    "Carga": float(carga)
                                })
                            except ValueError:
                                continue

                if cargas_orbitales:
                    datos_cargas[tipo.capitalize()] = pd.DataFrame(cargas_orbitales)

        return datos_cargas if datos_cargas else None
