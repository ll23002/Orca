# segundo intento
import pandas as pd
import re

def generar_entrada_orca(contenido_xyz, tipo_calculo, metodo, base, palabras_clave):
    """Genera el contenido del archivo de entrada para ORCA."""
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
    bloque_xyz = f"* xyz 0 1\n{coords_str}\n*\n"
    return encabezado + bloque_xyz

def extraer_geometria_optimizada(ruta_salida):
    """Extrae las coordenadas XYZ optimizadas del archivo de salida de ORCA."""
    try:
        with open(ruta_salida, 'r', encoding='utf-8', errors='ignore') as f:
            contenido = f.read()
    except FileNotFoundError: return None

    patron = r'CARTESIAN COORDINATES \(ANGSTROEM\)\s*\n\s*-+\s*\n((?:\s*\S+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s*\n)+)'
    coincidencias = list(re.finditer(patron, contenido))
    if not coincidencias: return None

    bloque_coords = coincidencias[-1].group(1).strip()
    lineas_coords = [linea.strip() for linea in bloque_coords.split('\n') if linea.strip()]
    if not lineas_coords: return None

    num_atomos = len(lineas_coords)
    bloque_xyz = f"{num_atomos}\nGeometría Optimizada\n"
    for linea in lineas_coords:
        partes = linea.split()
        if len(partes) >= 4:
            bloque_xyz += f"{partes[0]:<2} {float(partes[1]):>12.6f} {float(partes[2]):>12.6f} {float(partes[3]):>12.6f}\n"
    return bloque_xyz

def extraer_espectro_ir(ruta_salida, factor_escalamiento=1.0):
    """
    Extrae los datos del espectro IR y aplica un factor de escalamiento.
    Esta versión es más robusta para manejar diferentes formatos de cabecera en el archivo de salida de ORCA.
    """
    try:
        with open(ruta_salida, 'r', encoding='utf-8', errors='ignore') as f:
            contenido = f.read()
    except FileNotFoundError:
        return pd.DataFrame()

    # Patrón robusto para encontrar el bloque de datos del espectro IR.
    # Busca "IR SPECTRUM", salta cualquier línea de cabecera y el separador '---',
    # y luego captura todas las líneas de datos hasta encontrar una línea que empieza con '*' o una sección nueva.
    patron_ir_bloque = r'IR SPECTRUM\s*\n-+\n(?:.|\n)*?-+\n((?:.|\n)*?)(?=\n\s*\*|\n\s*-{2,}\n[A-Z]|\Z)'
    coincidencia = re.search(patron_ir_bloque, contenido)

    if not coincidencia:
        return pd.DataFrame()

    datos = []
    bloque_datos = coincidencia.group(1).strip()

    for linea in bloque_datos.split('\n'):
        partes = linea.split()
        # Una línea de datos válida tiene un formato como "6: 1602.68 ..."
        if len(partes) > 3 and partes[0].endswith(':'):
            try:
                # La frecuencia es la segunda columna (índice 1)
                freq = float(partes[1])
                # La intensidad es la cuarta columna (índice 3)
                intensidad = float(partes[3])

                # Ignorar modos de traslación/rotación que tienen frecuencia ~0
                if freq > 10.0:
                    datos.append({"Frequency": freq * factor_escalamiento, "Intensity": intensidad})
            except (ValueError, IndexError):
                # Ignorar líneas que no se pueden procesar
                continue

    return pd.DataFrame(datos)

def extraer_energia_final(ruta_salida):
    """Extrae la energía final (Single Point Energy) del archivo de salida."""
    try:
        with open(ruta_salida, 'r', encoding='utf-8', errors='ignore') as f:
            contenido = f.read()
    except FileNotFoundError: return None

    coincidencias = re.findall(r'FINAL SINGLE POINT ENERGY\s+([-\d.]+)', contenido)
    if coincidencias: return float(coincidencias[-1])
    return None

def extraer_componentes_energia(ruta_salida):
    """Extrae el desglose de componentes energéticos."""
    try:
        with open(ruta_salida, 'r', encoding='utf-8', errors='ignore') as f:
            contenido = f.read()
    except FileNotFoundError: return None

    patrones = {
        "Repulsión Nuclear": r'Nuclear Repulsion\s+:\s*([-\d.]+)',
        "Energía Electrónica": r'Electronic Energy\s+:\s*([-\d.]+)',
        "Energía Un Electrón": r'One Electron Energy\s+:\s*([-\d.]+)',
        "Energía Dos Electrones": r'Two Electron Energy\s+:\s*([-\d.]+)',
        "Energía Cinética": r'Kinetic Energy\s+:\s*([-\d.]+)',
        "Energía Potencial": r'Potential Energy\s+:\s*([-\d.]+)'
    }
    energias = {}
    for nombre, patron in patrones.items():
        coincidencias = re.findall(patron, contenido)
        if coincidencias: energias[nombre] = [float(coincidencias[-1])]

    return pd.DataFrame.from_dict(energias, orient='index', columns=['Energía (Hartree)']) if energias else None

def extraer_cargas_atomicas(ruta_salida):
    """Extrae las cargas atómicas de Mulliken y Loewdin."""
    try:
        with open(ruta_salida, 'r', encoding='utf-8', errors='ignore') as f:
            contenido = f.read()
    except FileNotFoundError: return None

    datos_cargas = {}
    for tipo in ['MULLIKEN', 'LOEWDIN']:
        patron = re.compile(rf'{tipo} ATOMIC CHARGES\s*\n-+\n((?:.|\n)*?)(?=\n\n|\Z)')
        coincidencia = re.search(patron, contenido)
        if coincidencia:
            cargas = []
            for linea in coincidencia.group(1).strip().split('\n'):
                partes = linea.split()
                if len(partes) == 4 and partes[0].endswith(':'):
                    cargas.append({"Átomo": f"{partes[0][:-1]} {partes[1]}", "Carga": float(partes[3])})
            datos_cargas[tipo.capitalize()] = pd.DataFrame(cargas)

    return datos_cargas if datos_cargas else None

def verificar_convergencia_optimizacion(ruta_salida):
    """Verifica si la optimización de la geometría ha convergido."""
    try:
        with open(ruta_salida, 'r', encoding='utf-8', errors='ignore') as f:
            contenido = f.read()
    except FileNotFoundError: return False
    return "THE OPTIMIZATION HAS CONVERGED" in contenido

def extraer_energias_orbitales(ruta_salida):
    """Extrae las energías de los orbitales."""
    try:
        with open(ruta_salida, 'r', encoding='utf-8', errors='ignore') as f:
            contenido = f.read()
    except FileNotFoundError: return None

    patron = re.compile(r'ORBITAL ENERGIES\s*\n-+\n((?:.|\n)*?)(?=\n\n|\Z)')
    coincidencia = re.search(patron, contenido)
    if not coincidencia: return None

    orbitales = []
    for linea in coincidencia.group(1).strip().split('\n')[2:]:  # Omitir cabeceras
        partes = linea.split()
        if len(partes) == 4:
            orbitales.append({
                "Número": int(partes[0]),
                "Ocupación": float(partes[1]),
                "Energía (Eh)": float(partes[2]),
                "Energía (eV)": float(partes[3])
            })
    return pd.DataFrame(orbitales)

def extraer_cargas_orbitales_reducidas(ruta_salida):
    """Extrae las cargas orbitales reducidas de Mulliken y Loewdin."""
    try:
        with open(ruta_salida, 'r', encoding='utf-8', errors='ignore') as f:
            contenido = f.read()
    except FileNotFoundError: return None

    datos_cargas = {}
    for tipo in ['MULLIKEN', 'LOEWDIN']:
        patron = re.compile(rf'{tipo} REDUCED ORBITAL CHARGES\s*\n-+\n((?:.|\n)*?)(?=\n\n|\Z)')
        coincidencia = re.search(patron, contenido)
        if coincidencia:
            cargas_orbitales = []
            atomo_actual = ""
            for linea in coincidencia.group(1).strip().split('\n'):
                if ":" in linea and not linea.strip().startswith(('s ', 'p ', 'd ')):
                    atomo_actual = linea.split(':')[0].strip()
                elif any(orb in linea for orb in ["s       :", "p       :", "d       :"]):
                    partes = [p for p in linea.split() if p != ':']
                    if len(partes) >= 3:
                        cargas_orbitales.append({
                            "Átomo": atomo_actual,
                            "Orbital": partes[0],
                            "Carga": float(partes[1])
                        })
            datos_cargas[tipo.capitalize()] = pd.DataFrame(cargas_orbitales)

    return datos_cargas if datos_cargas else None

