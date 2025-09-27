import pandas as pd
from pathlib import Path
import shutil
import json

# Importaciones correctas de la biblioteca OPI
from opi.core import Calculator
from opi.input.structures.structure import Structure


def ejecutar_y_procesar_orca(nombre_trabajo, contenido_xyz, tipo_calculo, metodo, base, palabras_clave_extra,
                             dir_calculos, factor_escalamiento=1.0):
    """
    Ejecuta un cálculo de ORCA y procesa la salida utilizando la biblioteca OPI,
    con una verificación de convergencia robusta y un parseo directo del JSON de resultados.
    """
    # 1. Configurar directorios y rutas
    working_dir = Path(dir_calculos) / nombre_trabajo
    if working_dir.exists():
        shutil.rmtree(working_dir)
    working_dir.mkdir(parents=True, exist_ok=True)
    ruta_salida = working_dir / f"{nombre_trabajo}.out"
    ruta_prop_json = working_dir / f"{nombre_trabajo}.property.json"

    # 2. Configurar la calculadora
    calc = Calculator(
        basename=nombre_trabajo,
        working_dir=working_dir,
    )

    # 3. Cargar la estructura
    ruta_xyz_entrada = working_dir / f"{nombre_trabajo}_input.xyz"
    with open(ruta_xyz_entrada, "w") as f:
        f.write(contenido_xyz)
    calc.structure = Structure.from_xyz(ruta_xyz_entrada)

    # 4. Añadir palabras clave
    input_keywords = [metodo, base] + palabras_clave_extra.split()
    if tipo_calculo == "Optimización de Geometría":
        input_keywords.append("OPT")
    elif tipo_calculo == "Frecuencias Vibracionales (IR)":
        input_keywords.append("OPT")
        input_keywords.append("FREQ")
    calc.input.add_simple_keywords(*input_keywords)

    # 5. Escribir el input y ejecutar ORCA
    try:
        calc.write_input()
        calc.run()
    except Exception as e:
        log_content = ""
        if ruta_salida.exists():
            with open(ruta_salida, 'r', encoding='utf-8', errors='ignore') as f:
                log_content = f.read()
        return {"error": f"Se produjo una excepción durante la ejecución de ORCA: {e}",
                "log_completo_orca": log_content}

    # 6. Obtener y procesar los resultados
    output = calc.get_output()
    if not output.terminated_normally():
        log_content = ""
        if ruta_salida.exists():
            with open(ruta_salida, 'r', encoding='utf-8', errors='ignore') as f:
                log_content = f.read()
        return {"error": "ORCA no terminó normalmente. Revisa el log.", "log_completo_orca": log_content}

    # Generar el JSON de propiedades
    output.create_property_json()

    resultados = {
        "opt_convergida": False, "xyz_optimizada": None, "energia_final": None,
        "datos_ir": None, "datos_cargas": None, "datos_orbitales": None,
        "log_completo_orca": "", "error": None
    }

    # Leer el log completo para la pestaña de "Datos Técnicos" y para la verificación
    with open(ruta_salida, 'r', encoding='utf-8', errors='ignore') as f:
        log_completo = f.read()
        resultados["log_completo_orca"] = log_completo

    # --- MÉTODO DE VERIFICACIÓN HÍBRIDO Y ROBUSTO ---

    # 1. Verificar la convergencia directamente desde el log de ORCA
    if "THE OPTIMIZATION HAS CONVERGED" in log_completo:
        resultados["opt_convergida"] = True

    # 2. Leer el JSON de propiedades manualmente para extraer los datos
    if ruta_prop_json.exists():
        with open(ruta_prop_json, 'r') as f:
            props = json.load(f)

        # La mayoría de los datos finales están en el último paso de la optimización
        if props and "Geometries" in props and props["Geometries"]:
            final_step = props["Geometries"][-1]

            # Extraer energía final
            if "Energy" in final_step and final_step["Energy"]:
                resultados["energia_final"] = final_step["Energy"][0]["totalEnergy"][0][0]

            # Extraer geometría optimizada del JSON
            if "Geometry" in final_step and "Coordinates" in final_step["Geometry"]:
                coords = final_step["Geometry"]["Coordinates"]
                if "Cartesians" in coords:
                    cartesians = coords["Cartesians"]
                    n_atoms = len(cartesians)

                    # Construir el formato XYZ
                    xyz_lines = [str(n_atoms), "Optimized geometry from ORCA"]

                    for atom_data in cartesians:
                        element = atom_data[0]
                        x, y, z = atom_data[1], atom_data[2], atom_data[3]
                        # Convertir de bohr a angstrom (1 bohr = 0.529177 angstrom)
                        x_ang = x * 0.529177249
                        y_ang = y * 0.529177249
                        z_ang = z * 0.529177249
                        xyz_lines.append(f"{element:2s} {x_ang:15.10f} {y_ang:15.10f} {z_ang:15.10f}")

                    resultados["xyz_optimizada"] = "\n".join(xyz_lines)

            # Extraer cargas de Mulliken
            if "Mulliken_Population_Analysis" in final_step:
                cargas = {}
                mulliken_data = final_step["Mulliken_Population_Analysis"][0]

                # Obtener los símbolos de los átomos de la geometría
                atom_symbols = []
                if "Geometry" in final_step and "Coordinates" in final_step["Geometry"]:
                    cartesians = final_step["Geometry"]["Coordinates"]["Cartesians"]
                    atom_symbols = [atom[0] for atom in cartesians]

                cargas['Mulliken'] = pd.DataFrame({
                    "Átomo": [f'{symbol}{i + 1}' for i, symbol in enumerate(atom_symbols)] if atom_symbols else [
                        f'Atom {i + 1}' for i in range(len(mulliken_data["AtomicCharges"]))],
                    "Carga": [c[0] for c in mulliken_data["AtomicCharges"]]
                })
                resultados["datos_cargas"] = cargas

            # Extraer cargas de Loewdin
            if "Loewdin_Population_Analysis" in final_step:
                cargas = resultados.get("datos_cargas", {})
                loewdin_data = final_step["Loewdin_Population_Analysis"][0]

                # Obtener los símbolos de los átomos de la geometría
                atom_symbols = []
                if "Geometry" in final_step and "Coordinates" in final_step["Geometry"]:
                    cartesians = final_step["Geometry"]["Coordinates"]["Cartesians"]
                    atom_symbols = [atom[0] for atom in cartesians]

                cargas['Loewdin'] = pd.DataFrame({
                    "Átomo": [f'{symbol}{i + 1}' for i, symbol in enumerate(atom_symbols)] if atom_symbols else [
                        f'Atom {i + 1}' for i in range(len(loewdin_data["AtomicCharges"]))],
                    "Carga": [c[0] for c in loewdin_data["AtomicCharges"]]
                })
                resultados["datos_cargas"] = cargas

    # Si por alguna razón el JSON no funcionó, intentar con OPI como respaldo
    if not resultados["xyz_optimizada"]:
        try:
            if hasattr(output, 'structures') and hasattr(output.structures,
                                                         'optimized') and output.structures.optimized:
                resultados["xyz_optimizada"] = output.structures.optimized.to_xyz()
        except Exception as e:
            print(f"Warning: No se pudo extraer geometría con OPI: {e}")

    # Verificar que se extrajo la energía desde el JSON, si no, intentar con OPI
    if not resultados["energia_final"]:
        try:
            if hasattr(output, 'final_energy') and output.final_energy:
                resultados["energia_final"] = output.final_energy
        except Exception as e:
            print(f"Warning: No se pudo extraer energía con OPI: {e}")

    return resultados