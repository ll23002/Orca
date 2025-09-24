import pandas as pd
import numpy as np
import os

try:
    from pyscf import gto, dft, scf, geomopt, hessian
    from pyscf.tools import cubegen
    from pyscf.geomopt.geometric_solver import optimize
    from pyscf.hessian import thermo

    PYSCF_DISPONIBLE = True
except ImportError:
    PYSCF_DISPONIBLE = False


def parsear_xyz_contenido(contenido_xyz):
    lineas = contenido_xyz.strip().split('\n')
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
    if not PYSCF_DISPONIBLE:
        raise ImportError("PySCF no está disponible")

    atom_string = [f"{elem} {x:.6f} {y:.6f} {z:.6f}" for elem, x, y, z in geometria]

    mapeo_bases = {
        '6-31+G(d,p)': '6-31+G(d,p)', '6-311++G(d,p)': '6-311++G(d,p)',
        'cc-pVDZ': 'cc-pvdz', 'def2-SVP': 'def2-svp'
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
    mapeo_metodos = {
        'B3LYP': 'b3lyp', 'PBE0': 'pbe0', 'M06-2X': 'm062x', 'wB97X-D': 'wb97x-d'
    }
    metodo_pyscf = mapeo_metodos.get(metodo, metodo.lower())

    mf = dft.RKS(mol) if mol.spin == 0 else dft.UKS(mol)
    mf.xc = metodo_pyscf
    return mf


def ejecutar_optimizacion_geometria(mol, metodo_obj):
    try:
        mol_eq = optimize(metodo_obj, maxsteps=100)
        metodo_final = configurar_metodo_dft(mol_eq, metodo_obj.xc)
        energia_final = metodo_final.kernel()

        return {
            'mol_optimizada': mol_eq, 'energia_final': energia_final,
            'convergido': metodo_final.converged, 'metodo_final': metodo_final,
            'geometria_optimizada': mol_eq.atom_coords(),
            'simbolos_atomicos': [mol_eq.atom_symbol(i) for i in range(mol_eq.natm)]
        }
    except Exception as e:
        print(f"Error durante optimización: {e}")
        energia_final = metodo_obj.kernel()
        return {
            'mol_optimizada': mol, 'energia_final': energia_final,
            'convergido': metodo_obj.converged, 'metodo_final': metodo_obj,
            'geometria_optimizada': mol.atom_coords(),
            'simbolos_atomicos': [mol.atom_symbol(i) for i in range(mol.natm)],
            'optimizacion_fallida': True
        }


def calcular_frecuencias_ir(mol_opt, metodo_obj):
    try:
        hess = hessian.RHF(metodo_obj) if mol_opt.spin == 0 else hessian.UHF(metodo_obj)
        hessiano_matrix = hess.kernel()
        freq_info = thermo.harmonic_analysis(mol_opt, hessiano_matrix)
        frecuencias_au = freq_info['freq_au']
        frecuencias_cm = np.real(frecuencias_au) * 219474.6

        intensidades_aproximadas = []
        for freq in frecuencias_cm:
            intensidad = max(1.0, abs(freq) * 0.01) if freq > 0 else 0.0
            intensidades_aproximadas.append(intensidad)

        return {
            'frecuencias': frecuencias_cm,
            'intensidades': np.array(intensidades_aproximadas),
            'freq_info': freq_info
        }
    except Exception as e:
        print(f"Error calculando frecuencias: {e}")
        return None


def calcular_componentes_energeticos(mol, metodo_obj):
    """
    Calcula todos los componentes energéticos incluyendo One Electron, Two Electron y Kinetic Energy
    """
    try:
        # Obtener la matriz de densidad
        dm = metodo_obj.make_rdm1()

        # Integrales de un electrón
        h1e = metodo_obj.get_hcore()  # Hamiltoniano de core (cinético + potencial nuclear)

        # Integrales cinéticas y de potencial nuclear por separado
        t = mol.intor('int1e_kin')  # Energía cinética
        v = mol.intor('int1e_nuc')  # Potencial núcleo-electrón

        # Energía de un electrón
        energia_un_electron = np.trace(np.dot(dm, h1e))

        # Energía cinética
        energia_cinetica = np.trace(np.dot(dm, t))

        # Energía potencial núcleo-electrón
        energia_potencial_nuclear = np.trace(np.dot(dm, v))

        # Energía de dos electrones (repulsión electrón-electrón)
        vhf = metodo_obj.get_veff()  # Potencial de Hartree-Fock/DFT
        j = metodo_obj.get_j()  # Integral de Coulomb

        # Para DFT, la energía de dos electrones incluye intercambio-correlación
        if hasattr(metodo_obj, 'get_veff'):
            # Energía electrón-electrón total (incluye XC para DFT)
            energia_dos_electrones = 0.5 * np.trace(np.dot(dm, vhf - h1e))
        else:
            # Para métodos HF puros
            energia_dos_electrones = 0.5 * np.trace(np.dot(dm, j))

        # Repulsión nuclear
        energia_repulsion_nuclear = mol.energy_nuc()

        # Energía electrónica total
        energia_electronica = energia_un_electron + energia_dos_electrones

        # Energía total
        energia_total = energia_electronica + energia_repulsion_nuclear

        componentes = {
            'Energía Total': energia_total,
            'Energía Electrónica': energia_electronica,
            'Energía Un Electrón': energia_un_electron,
            'Energía Dos Electrones': energia_dos_electrones,
            'Energía Cinética': energia_cinetica,
            'Potencial Núcleo-Electrón': energia_potencial_nuclear,
            'Repulsión Nuclear': energia_repulsion_nuclear
        }

        # Agregar información sobre orbitales si están disponibles
        if hasattr(metodo_obj, 'mo_energy') and metodo_obj.mo_energy is not None:
            # HOMO y LUMO
            mo_occ = metodo_obj.mo_occ
            mo_energy = metodo_obj.mo_energy

            if np.any(mo_occ > 0):
                homo_energy = mo_energy[mo_occ > 0][-1]
                componentes['Energía HOMO'] = homo_energy

            if np.any(mo_occ == 0):
                lumo_energy = mo_energy[mo_occ == 0][0]
                componentes['Energía LUMO'] = lumo_energy
                componentes['Gap HOMO-LUMO'] = lumo_energy - homo_energy

        return componentes

    except Exception as e:
        print(f"Error calculando componentes energéticos: {e}")
        return None


def ejecutar_calculo_pyscf(geometria, config, nombre_trabajo):
    if not PYSCF_DISPONIBLE:
        raise ImportError("PySCF no está disponible. Instala con: pip install pyscf")

    try:
        mol = configurar_molecula_pyscf(
            geometria, base=config['base'], carga=0, spin=0
        )
        metodo_obj = configurar_metodo_dft(mol, config['metodo'])
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

        # Calcular componentes energéticos detallados
        componentes_energia = calcular_componentes_energeticos(mol, metodo_obj)
        if componentes_energia:
            resultados['componentes_energia'] = componentes_energia

        try:
            resultados['dipole_moment'] = metodo_obj.dip_moment()
            resultados['mulliken_charges'] = metodo_obj.mulliken_charges()
            resultados['lowdin_charges'] = metodo_obj.mulliken_charges(method='lowdin')
        except Exception:
            pass

        if config['tipo_calculo'] in ["Optimización de Geometría", "Frecuencias Vibracionales (IR)"]:
            resultado_opt = ejecutar_optimizacion_geometria(mol, metodo_obj)
            resultados.update(resultado_opt)

            # Recalcular componentes energéticos para la geometría optimizada
            if 'metodo_final' in resultado_opt:
                componentes_opt = calcular_componentes_energeticos(
                    resultado_opt['mol_optimizada'],
                    resultado_opt['metodo_final']
                )
                if componentes_opt:
                    resultados['componentes_energia_optimizada'] = componentes_opt

            if config['tipo_calculo'] == "Frecuencias Vibracionales (IR)" and resultado_opt.get('convergido', False):
                freq_data = calcular_frecuencias_ir(
                    resultado_opt['mol_optimizada'], resultado_opt['metodo_final']
                )
                if freq_data:
                    resultados['frecuencias_data'] = freq_data

        resultados['calculo_completado'] = True
        resultados['convergido'] = resultados.get('convergido', metodo_obj.converged)

        # Guardar reporte detallado
        ruta_resultados = os.path.join("calculations", f"{nombre_trabajo}_pyscf.txt")
        with open(ruta_resultados, 'w') as f:
            f.write(generar_reporte_completo(resultados, nombre_trabajo))

        return resultados

    except Exception as e:
        print(f"Error en cálculo PySCF: {e}")
        return None


def extraer_energia_final(resultados):
    if resultados is None: return None
    return resultados.get('energia_final', resultados.get('energia_scf'))


def extraer_geometria_optimizada(resultados):
    if resultados is None or 'geometria_optimizada' not in resultados:
        return None
    try:
        coords = resultados['geometria_optimizada']
        simbolos = resultados['simbolos_atomicos']
        num_atomos = len(simbolos)
        xyz_content = f"{num_atomos}\nGeometría Optimizada con PySCF\n"
        for i, simbolo in enumerate(simbolos):
            x, y, z = coords[i]
            xyz_content += f"{simbolo:<2} {x:>12.6f} {y:>12.6f} {z:>12.6f}\n"
        return xyz_content
    except Exception as e:
        print(f"Error extrayendo geometría: {e}")
        return None


def extraer_componentes_energia(resultados):
    if resultados is None:
        return None

    try:
        # Priorizar componentes de geometría optimizada si están disponibles
        componentes_data = resultados.get('componentes_energia_optimizada',
                                          resultados.get('componentes_energia'))

        if componentes_data:
            # Convertir a DataFrame con formato apropiado
            componentes_df = {}
            for nombre, valor in componentes_data.items():
                if isinstance(valor, (int, float)):
                    componentes_df[nombre] = [float(valor)]
                else:
                    componentes_df[nombre] = [float(valor)]

            df = pd.DataFrame.from_dict(componentes_df, orient='index', columns=['Energía (Hartree)'])
            return df

        # Fallback a método anterior si no hay componentes detallados
        componentes = {}
        if 'energia_scf' in resultados:
            componentes['Energía SCF'] = [resultados['energia_scf']]
        if 'energia_final' in resultados:
            componentes['Energía Final'] = [resultados['energia_final']]

        if 'metodo_obj' in resultados:
            mol = resultados.get('mol_optimizada', resultados.get('mol_inicial'))
            if mol:
                energia_nuclear = mol.energy_nuc()
                componentes['Repulsión Nuclear'] = [energia_nuclear]
                if 'energia_final' in resultados:
                    energia_electronica = resultados['energia_final'] - energia_nuclear
                    componentes['Energía Electrónica'] = [energia_electronica]

        if componentes:
            return pd.DataFrame.from_dict(componentes, orient='index', columns=['Energía (Hartree)'])
        return None

    except Exception as e:
        print(f"Error extrayendo componentes energéticos: {e}")
        return None


def extraer_cargas_atomicas(resultados):
    if resultados is None: return None
    try:
        datos_cargas = {'Mulliken': pd.DataFrame(), 'Loewdin': pd.DataFrame()}
        simbolos = resultados.get('simbolos_atomicos', [])

        if 'mulliken_charges' in resultados:
            cargas_mulliken = [
                {"Átomo": f"{i + 1} {simbolos[i] if i < len(simbolos) else ''}", "Carga": float(carga)}
                for i, carga in enumerate(resultados['mulliken_charges'])
            ]
            datos_cargas['Mulliken'] = pd.DataFrame(cargas_mulliken)

        if 'lowdin_charges' in resultados:
            cargas_lowdin = [
                {"Átomo": f"{i + 1} {simbolos[i] if i < len(simbolos) else ''}", "Carga": float(carga)}
                for i, carga in enumerate(resultados['lowdin_charges'])
            ]
            datos_cargas['Loewdin'] = pd.DataFrame(cargas_lowdin)

        return datos_cargas if not (datos_cargas['Mulliken'].empty and datos_cargas['Loewdin'].empty) else None
    except Exception as e:
        print(f"Error extrayendo cargas atómicas: {e}")
        return None


def extraer_espectro_ir(resultados, factor_escalamiento=1.0):
    if resultados is None or 'frecuencias_data' not in resultados:
        return pd.DataFrame()
    try:
        freq_data = resultados['frecuencias_data']
        datos_ir = [
            {"Frequency": freq * factor_escalamiento, "Intensity": intensidad}
            for freq, intensidad in zip(freq_data['frecuencias'], freq_data['intensidades'])
            if freq > 10.0
        ]
        if not datos_ir: return pd.DataFrame()
        return pd.DataFrame(datos_ir).sort_values('Frequency').reset_index(drop=True)
    except Exception as e:
        print(f"Error extrayendo espectro IR: {e}")
        return pd.DataFrame()


def verificar_convergencia_optimizacion(resultados):
    if resultados is None: return False
    return resultados.get('convergido', False)


def generar_reporte_completo(resultados, nombre_trabajo):
    if resultados is None:
        return "No hay resultados disponibles."

    reporte = []
    reporte.append("=" * 70)
    reporte.append(f"REPORTE DETALLADO DE CÁLCULO PYSCF: {nombre_trabajo}")
    reporte.append("=" * 70)
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

    # Componentes energéticos detallados
    componentes_data = resultados.get('componentes_energia_optimizada',
                                      resultados.get('componentes_energia'))
    if componentes_data:
        reporte.append("ANÁLISIS ENERGÉTICO DETALLADO:")
        reporte.append("-" * 35)
        for nombre, valor in componentes_data.items():
            if isinstance(valor, (int, float)):
                reporte.append(f"{nombre:<25}: {valor:>15.8f} Hartree")
                if nombre in ['Energía HOMO', 'Energía LUMO']:
                    reporte.append(f"{' ' * 25}  {valor * 27.2114:>15.4f} eV")
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

    reporte.append("=" * 70)
    reporte.append("Cálculo completado con PySCF")
    reporte.append("=" * 70)

    return "\n".join(reporte)