import pandas as pd
import re
import os


def generate_orca_input(xyz_content, calc_type, method, basis_set, keywords):
    #Genera el contenido del archivo de entrada para ORCA
    base_keywords = f"! {method} {basis_set} {keywords}"

    if calc_type == "Optimización de Geometría":
        calc_keywords = "OPT"
    elif calc_type == "Frecuencias Vibracionales (IR)":
        calc_keywords = "OPT FREQ"
    else:
        calc_keywords = ""

    header = f"{base_keywords} {calc_keywords}\n"

    lines = xyz_content.strip().split('\n')
    coord_lines = lines[2:]
    coords_str = "\n".join(coord_lines)

    charge = 0
    multiplicity = 1

    xyz_block = f"* xyz {charge} {multiplicity}\n{coords_str}\n*\n"

    return header + xyz_block


def debug_file_sections(output_path):
    # Secciones están presentes en el archivo
    try:
        with open(output_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"Archivo no encontrado: {output_path}")
        return

    # Patrones comunes
    patterns_to_check = [
        "CARTESIAN COORDINATES (ANGSTROEM)",
        "CARTESIAN COORDINATES (A.U.)",
        "FINAL SINGLE POINT ENERGY",
        "THE OPTIMIZATION HAS CONVERGED",
        "OPTIMIZATION RUN DONE",
        "IR SPECTRUM",
        "VIBRATIONAL FREQUENCIES",
        "MULLIKEN ATOMIC CHARGES",
        "LOEWDIN ATOMIC CHARGES",
        "Electronic energy",
        "Nuclear Repulsion"
    ]

    print(f"=== DEBUG: Análisis del archivo {output_path} ===")
    for pattern in patterns_to_check:
        matches = len(re.findall(pattern, content, re.IGNORECASE))
        print(f"'{pattern}': {matches} coincidencias")



def parse_optimized_geometry(output_path):
    #Extrae las coordenadas XYZ optimizadas del archivo de salida de ORCA
    try:
        with open(output_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"DEBUG: Archivo no encontrado: {output_path}")
        return None

    debug_file_sections(output_path)

    try:
        # Patrones alternativos para buscar coordenadas
        patterns = [
            r'CARTESIAN COORDINATES \(ANGSTROEM\)\s*\n\s*-+\s*\n((?:\s*\S+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s*\n)+)',
            r'CARTESIAN COORDINATES \(A\.U\.\)\s*\n\s*-+\s*\n((?:\s*\S+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s*\n)+)',
            r'FINAL OPTIMIZED COORDINATES\s*\n.*?\n((?:\s*\S+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s*\n)+)'
        ]

        coord_block = None
        pattern_used = None

        for i, pattern in enumerate(patterns):
            matches = list(re.finditer(pattern, content, re.MULTILINE | re.DOTALL))
            if matches:
                coord_block = matches[-1].group(1).strip()
                pattern_used = i
                print(f"DEBUG: Usando patrón {i}, encontradas {len(matches)} coincidencias")
                break

        if not coord_block:
            print("DEBUG: No se encontraron coordenadas con ningún patrón")
            return None

        # Procesar el bloque de coordenadas
        coord_lines = [line.strip() for line in coord_block.split('\n') if line.strip()]

        if not coord_lines:
            print("DEBUG: Bloque de coordenadas vacío")
            return None

        num_atoms = len(coord_lines)
        xyz_block = f"{num_atoms}\nOptimized Geometry\n"

        print(f"DEBUG: Procesando {num_atoms} átomos")

        for i, line in enumerate(coord_lines):
            parts = line.split()
            if len(parts) >= 4:
                element = parts[0]
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                xyz_block += f"{element:<2} {x:>12.6f} {y:>12.6f} {z:>12.6f}\n"
                print(f"DEBUG: Átomo {i + 1}: {element} {x:.6f} {y:.6f} {z:.6f}")
            else:
                print(f"DEBUG: Línea mal formateada: {line}")

        return xyz_block

    except Exception as e:
        print(f"DEBUG: Error en parse_optimized_geometry: {e}")
        return None


def parse_ir_spectrum(output_path, scaling_factor=1.0):
    #Extrae los datos del espectro IR y aplica un factor de escalamiento
    try:
        with open(output_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"DEBUG: Archivo no encontrado: {output_path}")
        return pd.DataFrame()

    print("DEBUG IR: Iniciando búsqueda de espectro IR...")

    print("DEBUG: Buscando y mostrando TODAS las secciones IR...")

    ir_sections = re.findall(r'(IR SPECTRUM.*?)(?=\n\n[A-Z]|\nTimings|\n\*\*\*\*ORCA|$)', content,
                             re.MULTILINE | re.DOTALL)
    print(f"DEBUG: Encontradas {len(ir_sections)} secciones 'IR SPECTRUM'")

    for i, section in enumerate(ir_sections):
        print(f"\n=== DEBUG: Sección IR SPECTRUM {i + 1} ===")
        lines = section.split('\n')
        for line_num, line in enumerate(lines[:30]):
            print(f"  {line_num:2d}: {line}")
        if len(lines) > 30:
            print(f"  ... ({len(lines) - 30} líneas más)")

    # Buscar secciones VIBRATIONAL FREQUENCIES
    freq_sections = re.findall(r'(VIBRATIONAL FREQUENCIES.*?)(?=\n\n[A-Z]|\nIR SPECTRUM|\nTimings|\n\*\*\*\*ORCA|$)',
                               content, re.MULTILINE | re.DOTALL)
    print(f"\nDEBUG: Encontradas {len(freq_sections)} secciones 'VIBRATIONAL FREQUENCIES'")

    for i, section in enumerate(freq_sections):
        print(f"\n=== DEBUG: Sección VIBRATIONAL FREQUENCIES {i + 1} ===")
        lines = section.split('\n')
        for line_num, line in enumerate(lines[:30]):
            print(f"  {line_num:2d}: {line}")
        if len(lines) > 30:
            print(f"  ... ({len(lines) - 30} líneas más)")

    data = []

    patterns_to_try = [
        r'Mode\s+freq.*?\n.*?-+.*?\n((?:\s*\d+.*?\n)+)',
        r'(\d+)\s+([\d.]+)\s+.*?([\d.]+)',
        r'Mode:\s*(\d+)\s+Frequency:\s*([\d.]+)\s+.*?Intensity:\s*([\d.]+)',
    ]

    for i, pattern in enumerate(patterns_to_try):
        print(f"\nDEBUG: Probando patrón {i + 1}: {pattern[:50]}...")
        matches = re.findall(pattern, content, re.MULTILINE | re.DOTALL)
        print(f"DEBUG: Encontradas {len(matches)} coincidencias con patrón {i + 1}")

        if matches and isinstance(matches[0], str):
            lines = matches[-1].strip().split('\n')
            print(f"DEBUG: Procesando tabla con {len(lines)} líneas")
            for line_num, line in enumerate(lines[:10]):
                print(f"  Línea {line_num}: {line}")
        elif matches:
            print(f"DEBUG: Procesando {len(matches)} grupos")
            for j, match in enumerate(matches[:5]):
                print(f"  Grupo {j}: {match}")

    return pd.DataFrame()


def parse_final_energy(output_path):
    """Extrae la energía final del archivo de salida."""
    try:
        with open(output_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"DEBUG: Archivo no encontrado: {output_path}")
        return None

    # Múltiples patrones para buscar la energía final
    patterns = [
        r'FINAL SINGLE POINT ENERGY\s+([-\d.]+)',
        r'Total Energy\s*:\s*([-\d.]+)',
        r'E\(.*?\)\s*=\s*([-\d.]+)',
        r'Final energy\s*:\s*([-\d.]+)'
    ]

    for pattern in patterns:
        matches = re.findall(pattern, content)
        if matches:
            try:
                energy = float(matches[-1])
                print(f"DEBUG: Energía final encontrada: {energy}")
                return energy
            except ValueError:
                continue

    print("DEBUG: No se encontró energía final")
    return None


def parse_energy_components(output_path):
    try:
        with open(output_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"DEBUG: Archivo no encontrado: {output_path}")
        return None

    energies = {}

    energy_patterns = {
        "Electronic Energy": r'Electronic energy\s*:\s*([-\d.]+)',
        "Nuclear Repulsion": r'Nuclear Repulsion\s*:\s*([-\d.]+)',
        "One Electron Energy": r'One electron energy\s*:\s*([-\d.]+)',
        "Two Electron Energy": r'Two electron energy\s*:\s*([-\d.]+)',
        "Kinetic Energy": r'Kinetic Energy\s*:\s*([-\d.]+)',
        "Potential Energy": r'Potential Energy\s*:\s*([-\d.]+)'
    }

    for name, pattern in energy_patterns.items():
        matches = re.findall(pattern, content, re.IGNORECASE)
        if matches:
            try:
                energies[name] = [float(matches[-1])]
                print(f"DEBUG: {name}: {matches[-1]}")
            except ValueError:
                continue

    if energies:
        df = pd.DataFrame.from_dict(energies, orient='index', columns=['Energy (Hartree)'])
        return df

    print("DEBUG: No se encontraron componentes energéticos")
    return None


def parse_atomic_charges(output_path):
    #Extrae las cargas atómicas de Mulliken y Loewdin
    try:
        with open(output_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"DEBUG: Archivo no encontrado: {output_path}")
        return None

    charge_data = {}

    def _parse_charge_block(charge_type, pattern):
        matches = re.findall(pattern, content, re.MULTILINE | re.DOTALL)
        if not matches:
            print(f"DEBUG: No se encontraron cargas {charge_type}")
            return pd.DataFrame()

        charges = []
        charge_block = matches[-1]
        lines = charge_block.strip().split('\n')

        for line in lines:
            if not line.strip():
                continue
            parts = line.split()
            try:
                if len(parts) >= 3:
                    if ':' in line:
                        atom_info = line.split(':')[0].strip()
                        charge_str = line.split(':')[1].strip()
                        charge_val = float(charge_str)
                        charges.append({"Atom": atom_info, "Charge": charge_val})
                    elif len(parts) == 3 and parts[0].isdigit():
                        atom_idx, element, charge_val = parts[0], parts[1], float(parts[2])
                        charges.append({"Atom": f"{atom_idx} {element}", "Charge": charge_val})
            except (ValueError, IndexError):
                continue

        print(f"DEBUG: Encontradas {len(charges)} cargas {charge_type}")
        return pd.DataFrame(charges)

    # Patrones para Mulliken y Loewdin
    mulliken_pattern = r'MULLIKEN ATOMIC CHARGES\s*\n\s*-+\s*\n((?:.*?\n)*?)(?=\n\s*\n|\n\s*-|$)'
    loewdin_pattern = r'LOEWDIN ATOMIC CHARGES\s*\n\s*-+\s*\n((?:.*?\n)*?)(?=\n\s*\n|\n\s*-|$)'

    charge_data['Mulliken'] = _parse_charge_block('Mulliken', mulliken_pattern)
    charge_data['Loewdin'] = _parse_charge_block('Loewdin', loewdin_pattern)

    return charge_data if not (charge_data['Mulliken'].empty and charge_data['Loewdin'].empty) else None


def check_optimization_convergence(output_path):
    try:
        with open(output_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"DEBUG: Archivo no encontrado: {output_path}")
        return False

    # Verificar convergencia
    convergence_patterns = [
        "THE OPTIMIZATION HAS CONVERGED",
        "OPTIMIZATION RUN DONE",
        "GEOMETRY OPTIMIZATION COMPLETED",
        "CONVERGED"
    ]

    for pattern in convergence_patterns:
        if re.search(pattern, content, re.IGNORECASE):
            print(f"DEBUG: Convergencia detectada con patrón: {pattern}")
            return True

    print("DEBUG: No se detectó convergencia")
    return False