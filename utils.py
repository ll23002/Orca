import pandas as pd
import re


def generar_entrada_orca(contenido_xyz, tipo_calculo, metodo, base, palabras_clave):

    palabras_base = f"! {metodo} {base} {palabras_clave}"

    if tipo_calculo == "Optimización de Geometría":
        palabras_calculo = "OPT"
    elif tipo_calculo == "Frecuencias Vibracionales (IR)":
        palabras_calculo = "OPT FREQ"
    else:
        palabras_calculo = ""

    encabezado = f"{palabras_base} {palabras_calculo}\n"

    lineas = contenido_xyz.strip().split('\n')
    lineas_coords = lineas[2:]
    coords_str = "\n".join(lineas_coords)

    carga = 0
    multiplicidad = 1

    bloque_xyz = f"* xyz {carga} {multiplicidad}\n{coords_str}\n*\n"

    return encabezado + bloque_xyz


def extraer_geometria_optimizada(ruta_salida):
    try:
        with open(ruta_salida, 'r', encoding='utf-8', errors='ignore') as f:
            contenido = f.read()
    except FileNotFoundError:
        return None

    try:
        patrones = [
            r'CARTESIAN COORDINATES \(ANGSTROEM\)\s*\n\s*-+\s*\n((?:\s*\S+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s*\n)+)',
            r'CARTESIAN COORDINATES \(A\.U\.\)\s*\n\s*-+\s*\n((?:\s*\S+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s*\n)+)',
            r'FINAL OPTIMIZED COORDINATES\s*\n.*?\n((?:\s*\S+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s*\n)+)'
        ]

        bloque_coords = None

        for i, patron in enumerate(patrones):
            coincidencias = list(re.finditer(patron, contenido, re.MULTILINE | re.DOTALL))
            if coincidencias:
                bloque_coords = coincidencias[-1].group(1).strip()
                break

        if not bloque_coords:
            return None

        lineas_coords = [linea.strip() for linea in bloque_coords.split('\n') if linea.strip()]

        if not lineas_coords:
            return None

        num_atomos = len(lineas_coords)
        bloque_xyz = f"{num_atomos}\nGeometría Optimizada\n"

        for linea in lineas_coords:
            partes = linea.split()
            if len(partes) >= 4:
                elemento = partes[0]
                x, y, z = float(partes[1]), float(partes[2]), float(partes[3])
                bloque_xyz += f"{elemento:<2} {x:>12.6f} {y:>12.6f} {z:>12.6f}\n"

        return bloque_xyz

    except Exception:
        return None


def extraer_espectro_ir(ruta_salida, factor_escalamiento=1.0):
    try:
        with open(ruta_salida, 'r', encoding='utf-8', errors='ignore') as f:
            contenido = f.read()
    except FileNotFoundError:
        return pd.DataFrame()

    datos = []

    patron_ir = r'(\d+):\s+([\d.]+)\s+[\d.]+\s+([\d.]+)\s+[\d.]+\s+\('

    coincidencias = re.findall(patron_ir, contenido)

    for coincidencia in coincidencias:
        try:
            modo = int(coincidencia[0])
            freq_calculada = float(coincidencia[1])
            intensidad = float(coincidencia[2])

            if freq_calculada > 10.0 and intensidad > 0.1:
                freq_escalada = freq_calculada * factor_escalamiento
                datos.append({
                    "Modo": modo,
                    "Frequency": freq_escalada,
                    "Intensity": intensidad
                })

        except (ValueError, IndexError):
            continue

    if not datos:
        return pd.DataFrame()

    df = pd.DataFrame(datos)
    df = df.sort_values('Frequency').reset_index(drop=True)
    return df[['Frequency', 'Intensity']]


def extraer_energia_final(ruta_salida):
    try:
        with open(ruta_salida, 'r', encoding='utf-8', errors='ignore') as f:
            contenido = f.read()
    except FileNotFoundError:
        return None

    patrones = [
        r'FINAL SINGLE POINT ENERGY\s+([-\d.]+)',
        r'Total Energy\s*:\s*([-\d.]+)',
        r'E\(.*?\)\s*=\s*([-\d.]+)',
        r'Final energy\s*:\s*([-\d.]+)'
    ]

    for patron in patrones:
        coincidencias = re.findall(patron, contenido)
        if coincidencias:
            try:
                return float(coincidencias[-1])
            except ValueError:
                continue

    return None


def extraer_componentes_energia(ruta_salida):
    try:
        with open(ruta_salida, 'r', encoding='utf-8', errors='ignore') as f:
            contenido = f.read()
    except FileNotFoundError:
        return None

    energias = {}

    patrones_energia = {
        "Energía Electrónica": r'Electronic energy\s*:\s*([-\d.]+)',
        "Repulsión Nuclear": r'Nuclear Repulsion\s*:\s*([-\d.]+)',
        "Energía Un Electrón": r'One electron energy\s*:\s*([-\d.]+)',
        "Energía Dos Electrones": r'Two electron energy\s*:\s*([-\d.]+)',
        "Energía Cinética": r'Kinetic Energy\s*:\s*([-\d.]+)',
        "Energía Potencial": r'Potential Energy\s*:\s*([-\d.]+)'
    }

    for nombre, patron in patrones_energia.items():
        coincidencias = re.findall(patron, contenido, re.IGNORECASE)
        if coincidencias:
            try:
                energias[nombre] = [float(coincidencias[-1])]
            except ValueError:
                continue

    if energias:
        df = pd.DataFrame.from_dict(energias, orient='index', columns=['Energía (Hartree)'])
        return df

    return None


def extraer_cargas_atomicas(ruta_salida):
    try:
        with open(ruta_salida, 'r', encoding='utf-8', errors='ignore') as f:
            contenido = f.read()
    except FileNotFoundError:
        return None

    datos_cargas = {}

    def _extraer_bloque_cargas(tipo_carga, patron):
        coincidencias = re.findall(patron, contenido, re.MULTILINE | re.DOTALL)
        if not coincidencias:
            return pd.DataFrame()

        cargas = []
        bloque_cargas = coincidencias[-1]
        lineas = bloque_cargas.strip().split('\n')

        for linea in lineas:
            if not linea.strip():
                continue
            partes = linea.split()
            try:
                if len(partes) >= 3:
                    if ':' in linea:
                        info_atomo = linea.split(':')[0].strip()
                        valor_carga = float(linea.split(':')[1].strip())
                        cargas.append({"Átomo": info_atomo, "Carga": valor_carga})
                    elif len(partes) == 3 and partes[0].isdigit():
                        idx_atomo, elemento, valor_carga = partes[0], partes[1], float(partes[2])
                        cargas.append({"Átomo": f"{idx_atomo} {elemento}", "Carga": valor_carga})
            except (ValueError, IndexError):
                continue

        return pd.DataFrame(cargas)

    patron_mulliken = r'MULLIKEN ATOMIC CHARGES\s*\n\s*-+\s*\n((?:.*?\n)*?)(?=\n\s*\n|\n\s*-|$)'
    patron_loewdin = r'LOEWDIN ATOMIC CHARGES\s*\n\s*-+\s*\n((?:.*?\n)*?)(?=\n\s*\n|\n\s*-|$)'

    datos_cargas['Mulliken'] = _extraer_bloque_cargas('Mulliken', patron_mulliken)
    datos_cargas['Loewdin'] = _extraer_bloque_cargas('Loewdin', patron_loewdin)

    return datos_cargas if not (datos_cargas['Mulliken'].empty and datos_cargas['Loewdin'].empty) else None


def verificar_convergencia_optimizacion(ruta_salida):
    try:
        with open(ruta_salida, 'r', encoding='utf-8', errors='ignore') as f:
            contenido = f.read()
    except FileNotFoundError:
        return False

    patrones_convergencia = [
        "THE OPTIMIZATION HAS CONVERGED",
        "OPTIMIZATION RUN DONE",
        "GEOMETRY OPTIMIZATION COMPLETED",
        "CONVERGED"
    ]

    for patron in patrones_convergencia:
        if re.search(patron, contenido, re.IGNORECASE):
            return True

    return False