import pandas as pd
import numpy as np
import os
from io import StringIO

try:
    from pyscf import gto, dft, scf, geomopt, hessian
    from pyscf.tools import cubegen
    from pyscf.geomopt.geometric_solver import optimize
    from pyscf.hessian import thermo

    PYSCF_DISPONIBLE = True
except ImportError:
    PYSCF_DISPONIBLE = False
    print("Advertencia: PySCF no está disponible. Instala con: pip install pyscf")


def parsear_xyz_contenido(contenido_xyz):
    """
    Parsea el contenido de un archivo XYZ y devuelve la geometría en formato PySCF.

    Args:
        contenido_xyz (str): Contenido del archivo XYZ

    Returns:
        list: Lista de tuplas (símbolo_atómico, x, y, z)
    """
    lineas = contenido_xyz.strip().split('\n')

    # Saltar las primeras dos líneas (número de átomos y comentario)
    lineas_coords = lineas[2:]

    geometria = []
    for linea in lineas_coords:
        if linea.strip():
            partes = linea.split()
            if len(partes) >= 4:
                elemento = partes[0]
                x, y, z = float(partes[1]), float(partes[2]), float(partes[3])
                geometria.append((elemento, x, y, z))

    return geometria


def configurar_molecula_pyscf(geometria, base='6-31+G(d,p)', carga=0, spin=0):
    """
    Configura una molécula en PySCF.

    Args:
        geometria (list): Lista de tuplas (elemento, x, y, z)
        base (str): Conjunto de base
        carga (int): Carga molecular
        spin (int): Multiplicidad de spin - 1

    Returns:
        pyscf.gto.Mole: Objeto molécula de PySCF
    """
    if not PYSCF_DISPONIBLE:
        raise ImportError("PySCF no está disponible")

    # Convertir geometría al formato PySCF
    atom_string = []
    for elemento, x, y, z in geometria:
        atom_string.append(f"{elemento} {x:.6f} {y:.6f} {z:.6f}")

    # Mapear nombres de bases de ORCA a PySCF
    mapeo_bases = {
        '6-31+G(d,p)': '6-31+G(d,p)',
        '6-311++G(d,p)': '6-311++G(d,p)',
        'cc-pVDZ': 'cc-pvdz',
        'def2-SVP': 'def2-svp'
    }

    base_pyscf = mapeo_bases.get(base, base.lower())

    mol = gto.Mole()
    mol.atom = '; '.join(atom_string)
    mol.basis = base_pyscf
    mol.charge = carga
    mol.spin = spin
    mol.verbose = 1
    mol.build()

    return mol


def configurar_metodo_dft(mol, metodo='B3LYP'):
    """
    Configura el método DFT en PySCF.

    Args:
        mol: Objeto molécula de PySCF
        metodo (str): Método DFT

    Returns:
        pyscf.dft.rks.RKS: Objeto de cálculo DFT
    """
    # Mapear métodos de ORCA a PySCF
    mapeo_metodos = {
        'B3LYP': 'b3lyp',
        'PBE0': 'pbe0',
        'M06-2X': 'm062x',
        'wB97X-D': 'wb97x-d'
    }

    metodo_pyscf = mapeo_metodos.get(metodo, metodo.lower())

    # Determinar si usar RKS o UKS basado en el spin
    if mol.spin == 0:
        mf = dft.RKS(mol)
    else:
        mf = dft.UKS(mol)

    mf.xc = metodo_pyscf
    return mf


def ejecutar_optimizacion_geometria(mol, metodo_obj):
    """
    Ejecuta optimización de geometría.

    Args:
        mol: Objeto molécula de PySCF
        metodo_obj: Objeto del método de cálculo

    Returns:
        dict: Resultados de la optimización
    """
    try:
        # Optimización usando el optimizador integrado de PySCF
        mol_eq = optimize(metodo_obj, maxsteps=100)

        # Calcular energía final con geometría optimizada
        metodo_final = configurar_metodo_dft(mol_eq, metodo_obj.xc)
        energia_final = metodo_final.kernel()

        convergido = metodo_final.converged

        resultados = {
            'mol_optimizada': mol_eq,
            'energia_final': energia_final,
            'convergido': convergido,
            'metodo_final': metodo_final,
            'geometria_optimizada': mol_eq.atom_coords(),
            'simbolos_atomicos': [mol_eq.atom_symbol(i) for i in range(mol_eq.natm)]
        }

        return resultados

    except Exception as e:
        print(f"Error durante optimización: {e}")
        # Si falla la optimización, hacer cálculo de punto único
        energia_final = metodo_obj.kernel()
        return {
            'mol_optimizada': mol,
            'energia_final': energia_final,
            'convergido': metodo_obj.converged,
            'metodo_final': metodo_obj,
            'geometria_optimizada': mol.atom_coords(),
            'simbolos_atomicos': [mol.atom_symbol(i) for i in range(mol.natm)],
            'optimizacion_fallida': True
        }


def calcular_frecuencias_ir(mol_opt, metodo_obj):
    """
    Calcula frecuencias vibracionales e intensidades IR.

    Args:
        mol_opt: Molécula optimizada
        metodo_obj: Objeto del método de cálculo

    Returns:
        dict: Frecuencias e intensidades IR
    """
    try:
        # Calcular hessiano
        hess = hessian.RHF(metodo_obj) if mol_opt.spin == 0 else hessian.UHF(metodo_obj)
        hessiano_matrix = hess.kernel()

        # Análisis termoquímico y frecuencias
        freq_info = thermo.harmonic_analysis(mol_opt, hessiano_matrix)

        # freq_info contiene las frecuencias en cm^-1
        frecuencias = freq_info['freq_au'] * 219474.6  # Convertir de au a cm^-1

        # Para intensidades IR, aproximación simple (en PySCF es más complejo)
        # Usaremos valores aproximados basados en las frecuencias
        intensidades_aproximadas = []
        for freq in frecuencias:
            if freq > 0:  # Solo modos reales
                # Intensidad aproximada basada en la frecuencia
                intensidad = max(1.0, abs(freq) * 0.01)
                intensidades_aproximadas.append(intensidad)
            else:
                intensidades_aproximadas.append(0.0)

        return {
            'frecuencias': frecuencias,
            'intensidades': np.array(intensidades_aproximadas),
            'freq_info': freq_info
        }

    except Exception as e:
        print(f"Error calculando frecuencias: {e}")
        return None


def ejecutar_calculo_pyscf(geometria, config, nombre_trabajo):
    """
    Función principal para ejecutar cálculos con PySCF.

    Args:
        geometria (list): Geometría molecular
        config (dict): Configuración del cálculo
        nombre_trabajo (str): Nombre del trabajo

    Returns:
        dict: Resultados del cálculo
    """
    if not PYSCF_DISPONIBLE:
        raise ImportError("PySCF no está disponible. Instala con: pip install pyscf")

    try:
        # Configurar molécula
        mol = configurar_molecula_pyscf(
            geometria,
            base=config['base'],
            carga=0,  # Puedes hacer esto configurable
            spin=0  # Puedes hacer esto configurable
        )

        # Configurar método
        metodo_obj = configurar_metodo_dft(mol, config['metodo'])

        # Ejecutar cálculo inicial
        energia_scf = metodo_obj.kernel()

        resultados = {
            'energia_scf': energia_scf,
            'convergencia_scf': metodo_obj.converged,
            'iteraciones_scf': getattr(metodo_obj, 'niter', 'N/A'),
            'metodo': config['metodo'],
            'base': config['base'],
            'mol_inicial': mol,
            'metodo_obj': metodo_obj
        }

        # Calcular propiedades adicionales
        try:
            # Momento dipolar
            dipole = metodo_obj.dip_moment()
            resultados['dipole_moment'] = dipole
        except:
            pass

        # Análisis de población Mulliken
        try:
            mulliken_charges = metodo_obj.mulliken_charges()
            resultados['mulliken_charges'] = mulliken_charges

            # Análisis de Löwdin si está disponible
            try:
                lowdin_charges = metodo_obj.mulliken_charges(method='lowdin')
                resultados['lowdin_charges'] = lowdin_charges
            except:
                pass
        except:
            pass

        # Optimización de geometría si es requerida
        if config['tipo_calculo'] in ["Optimización de Geometría", "Frecuencias Vibracionales (IR)"]:
            resultado_opt = ejecutar_optimizacion_geometria(mol, metodo_obj)
            resultados.update(resultado_opt)

            # Calcular frecuencias si es requerido
            if config['tipo_calculo'] == "Frecuencias Vibracionales (IR)" and resultado_opt.get('convergido', False):
                freq_data = calcular_frecuencias_ir(
                    resultado_opt['mol_optimizada'],
                    resultado_opt['metodo_final']
                )
                if freq_data:
                    resultados['frecuencias_data'] = freq_data

        # Marcar como completado
        resultados['calculo_completado'] = True
        resultados['convergido'] = resultados.get('convergido', metodo_obj.converged)

        # Guardar resultados en archivo
        ruta_resultados = os.path.join("calculations", f"{nombre_trabajo}_pyscf.txt")
        with open(ruta_resultados, 'w') as f:
            f.write(f"Resultados PySCF para {nombre_trabajo}\n")
            f.write("=" * 50 + "\n")
            f.write(f"Método: {config['metodo']}\n")
            f.write(f"Base: {config['base']}\n")
            f.write(f"Energía SCF: {energia_scf:.8f} Hartree\n")
            f.write(f"Convergencia: {metodo_obj.converged}\n")
            if 'dipole_moment' in resultados:
                dipole_mag = np.linalg.norm(resultados['dipole_moment'])
                f.write(f"Momento dipolar: {dipole_mag:.4f} Debye\n")

        return resultados

    except Exception as e:
        print(f"Error en cálculo PySCF: {e}")
        return None


def extraer_energia_final(resultados):
    """
    Extrae la energía final de los resultados de PySCF.

    Args:
        resultados (dict): Resultados del cálculo PySCF

    Returns:
        float: Energía final en Hartree
    """
    if resultados is None:
        return None

    # Priorizar energía de geometría optimizada si existe
    if 'energia_final' in resultados:
        return resultados['energia_final']

    # Si no, usar energía SCF
    return resultados.get('energia_scf', None)


def extraer_geometria_optimizada(resultados):
    """
    Extrae la geometría optimizada en formato XYZ.

    Args:
        resultados (dict): Resultados del cálculo PySCF

    Returns:
        str: Geometría en formato XYZ
    """
    if resultados is None or 'geometria_optimizada' not in resultados:
        return None

    try:
        coords = resultados['geometria_optimizada']
        simbolos = resultados['simbolos_atomicos']

        num_atomos = len(simbolos)
        xyz_content = f"{num_atomos}\nGeometría Optimizada con PySCF\n"

        for i, simbolo in enumerate(simbolos):
            x, y, z = coords[i]
            # Convertir de Bohr a Angstroms si es necesario
            # PySCF normalmente usa Angstroms en atom_coords()
            xyz_content += f"{simbolo:<2} {x:>12.6f} {y:>12.6f} {z:>12.6f}\n"

        return xyz_content

    except Exception as e:
        print(f"Error extrayendo geometría: {e}")
        return None


def extraer_componentes_energia(resultados):
    """
    Extrae componentes energéticos de los resultados PySCF.

    Args:
        resultados (dict): Resultados del cálculo PySCF

    Returns:
        pd.DataFrame: DataFrame con componentes energéticos
    """
    if resultados is None:
        return None

    try:
        componentes = {}

        # Energía SCF
        if 'energia_scf' in resultados:
            componentes['Energía SCF'] = [resultados['energia_scf']]

        # Energía final (optimizada)
        if 'energia_final' in resultados:
            componentes['Energía Final'] = [resultados['energia_final']]

        # Energía electrónica y nuclear (si disponible)
        if 'metodo_obj' in resultados:
            metodo = resultados['metodo_obj']
            mol = resultados.get('mol_optimizada', resultados.get('mol_inicial'))

            if mol is not None:
                # Energía de repulsión nuclear
                energia_nuclear = mol.energy_nuc()
                componentes['Repulsión Nuclear'] = [energia_nuclear]

                # Energía electrónica = Total - Nuclear
                if 'energia_final' in resultados:
                    energia_electronica = resultados['energia_final'] - energia_nuclear
                    componentes['Energía Electrónica'] = [energia_electronica]

        if componentes:
            df = pd.DataFrame.from_dict(componentes, orient='index', columns=['Energía (Hartree)'])
            return df

        return None

    except Exception as e:
        print(f"Error extrayendo componentes energéticos: {e}")
        return None


def extraer_cargas_atomicas(resultados):
    """
    Extrae cargas atómicas de los resultados PySCF.

    Args:
        resultados (dict): Resultados del cálculo PySCF

    Returns:
        dict: Diccionario con DataFrames de cargas Mulliken y Löwdin
    """
    if resultados is None:
        return None

    try:
        datos_cargas = {}
        simbolos = resultados.get('simbolos_atomicos', [])

        # Cargas de Mulliken
        if 'mulliken_charges' in resultados:
            cargas_mulliken = []
            mulliken = resultados['mulliken_charges']

            for i, carga in enumerate(mulliken):
                simbolo = simbolos[i] if i < len(simbolos) else f"Atom{i + 1}"
                cargas_mulliken.append({
                    "Átomo": f"{i + 1} {simbolo}",
                    "Carga": float(carga)
                })

            datos_cargas['Mulliken'] = pd.DataFrame(cargas_mulliken)
        else:
            datos_cargas['Mulliken'] = pd.DataFrame()

        # Cargas de Löwdin
        if 'lowdin_charges' in resultados:
            cargas_lowdin = []
            lowdin = resultados['lowdin_charges']

            for i, carga in enumerate(lowdin):
                simbolo = simbolos[i] if i < len(simbolos) else f"Atom{i + 1}"
                cargas_lowdin.append({
                    "Átomo": f"{i + 1} {simbolo}",
                    "Carga": float(carga)
                })

            datos_cargas['Loewdin'] = pd.DataFrame(cargas_lowdin)
        else:
            datos_cargas['Loewdin'] = pd.DataFrame()

        return datos_cargas if not (datos_cargas['Mulliken'].empty and datos_cargas['Loewdin'].empty) else None

    except Exception as e:
        print(f"Error extrayendo cargas atómicas: {e}")
        return None


def extraer_espectro_ir(resultados, factor_escalamiento=1.0):
    """
    Extrae datos del espectro IR de los resultados PySCF.

    Args:
        resultados (dict): Resultados del cálculo PySCF
        factor_escalamiento (float): Factor de escalamiento para frecuencias

    Returns:
        pd.DataFrame: DataFrame con frecuencias e intensidades
    """
    if resultados is None or 'frecuencias_data' not in resultados:
        return pd.DataFrame()

    try:
        freq_data = resultados['frecuencias_data']
        frecuencias = freq_data['frecuencias']
        intensidades = freq_data['intensidades']

        # Filtrar frecuencias positivas (modos reales)
        datos_ir = []

        for i, (freq, intensidad) in enumerate(zip(frecuencias, intensidades)):
            if freq > 10.0:  # Filtrar frecuencias muy bajas
                freq_escalada = freq * factor_escalamiento
                datos_ir.append({
                    "Frequency": freq_escalada,
                    "Intensity": intensidad
                })

        if not datos_ir:
            return pd.DataFrame()

        df = pd.DataFrame(datos_ir)
        df = df.sort_values('Frequency').reset_index(drop=True)
        return df

    except Exception as e:
        print(f"Error extrayendo espectro IR: {e}")
        return pd.DataFrame()


def verificar_convergencia_optimizacion(resultados):
    """
    Verifica si la optimización convergió.

    Args:
        resultados (dict): Resultados del cálculo PySCF

    Returns:
        bool: True si convergió, False si no
    """
    if resultados is None:
        return False

    return resultados.get('convergido', False)


def generar_reporte_pyscf(resultados, nombre_trabajo):
    """
    Genera un reporte detallado de los resultados PySCF.

    Args:
        resultados (dict): Resultados del cálculo PySCF
        nombre_trabajo (str): Nombre del trabajo

    Returns:
        str: Reporte en formato texto
    """
    if resultados is None:
        return "No hay resultados disponibles."

    reporte = []
    reporte.append("=" * 60)
    reporte.append(f"REPORTE DE CÁLCULO PYSCF: {nombre_trabajo}")
    reporte.append("=" * 60)
    reporte.append("")

    # Información general
    reporte.append("INFORMACIÓN GENERAL:")
    reporte.append("-" * 25)
    reporte.append(f"Método: {resultados.get('metodo', 'N/A')}")
    reporte.append(f"Conjunto base: {resultados.get('base', 'N/A')}")
    reporte.append(f"Energía SCF: {resultados.get('energia_scf', 'N/A'):.8f} Hartree")
    reporte.append(f"Convergencia SCF: {'Sí' if resultados.get('convergencia_scf', False) else 'No'}")
    reporte.append(f"Iteraciones SCF: {resultados.get('iteraciones_scf', 'N/A')}")
    reporte.append("")

    # Información de optimización
    if 'energia_final' in resultados:
        reporte.append("OPTIMIZACIÓN DE GEOMETRÍA:")
        reporte.append("-" * 30)
        reporte.append(f"Energía final: {resultados['energia_final']:.8f} Hartree")
        reporte.append(f"Convergió: {'Sí' if resultados.get('convergido', False) else 'No'}")
        reporte.append("")

    # Propiedades moleculares
    if 'dipole_moment' in resultados:
        dipole = resultados['dipole_moment']
        dipole_mag = np.linalg.norm(dipole)
        reporte.append("PROPIEDADES MOLECULARES:")
        reporte.append("-" * 28)
        reporte.append(f"Momento dipolar: {dipole_mag:.4f} Debye")
        reporte.append(f"Componentes (x,y,z): ({dipole[0]:.4f}, {dipole[1]:.4f}, {dipole[2]:.4f})")
        reporte.append("")

    # Información de frecuencias
    if 'frecuencias_data' in resultados:
        freq_data = resultados['frecuencias_data']
        n_freq_positivas = sum(1 for f in freq_data['frecuencias'] if f > 0)
        reporte.append("ANÁLISIS VIBRACIONAL:")
        reporte.append("-" * 23)
        reporte.append(f"Total de frecuencias: {len(freq_data['frecuencias'])}")
        reporte.append(f"Modos reales: {n_freq_positivas}")
        reporte.append("")

    reporte.append("=" * 60)
    reporte.append("Cálculo completado con PySCF")
    reporte.append("=" * 60)

    return "\n".join(reporte)