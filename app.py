import streamlit as st
import os
import pandas as pd
import numpy as np
import py3Dmol
from stmol import showmol
import matplotlib.pyplot as plt
from utils import (
    ejecutar_calculo_pyscf,
    extraer_geometria_optimizada,
    extraer_energia_final,
    extraer_componentes_energia,
    extraer_cargas_atomicas,
    extraer_espectro_ir,
    parsear_xyz_contenido
)

st.set_page_config(
    page_title="PySCF Molecular Calculator",
    layout="wide",
    initial_sidebar_state="expanded"
)

col1, col2 = st.columns([3, 1])
with col1:
    st.title("🧬 PySCF Molecular Calculator")
    st.markdown("*Calculadora cuántica para análisis molecular*")

# Inicializar estados de sesión
if "calculo_completado" not in st.session_state:
    st.session_state.calculo_completado = False
if "opt_convergida" not in st.session_state:
    st.session_state.opt_convergida = False
if "ultimo_tipo_calculo" not in st.session_state:
    st.session_state.ultimo_tipo_calculo = None
if "xyz_inicial" not in st.session_state:
    st.session_state.xyz_inicial = None
if "xyz_optimizada" not in st.session_state:
    st.session_state.xyz_optimizada = None
if "energia_final" not in st.session_state:
    st.session_state.energia_final = None
if "datos_energia" not in st.session_state:
    st.session_state.datos_energia = None
if "datos_cargas" not in st.session_state:
    st.session_state.datos_cargas = None
if "datos_ir" not in st.session_state:
    st.session_state.datos_ir = None
if "nombre_trabajo" not in st.session_state:
    st.session_state.nombre_trabajo = ""
if "resultados_pyscf" not in st.session_state:
    st.session_state.resultados_pyscf = None

DIR_CALCULOS = "calculations"
os.makedirs(DIR_CALCULOS, exist_ok=True)

with st.sidebar:
    st.markdown("### ⚛️ Panel de Control")
    st.markdown("---")

    st.markdown("#### 📁 **Molécula de Entrada**")
    archivo_subido = st.file_uploader(
        "Selecciona archivo .xyz",
        type=["xyz"],
        help="Carga la geometría inicial de tu molécula"
    )

    if archivo_subido is not None:
        st.session_state.xyz_inicial = archivo_subido.getvalue().decode("utf-8")
        st.session_state.nombre_trabajo = os.path.splitext(archivo_subido.name)[0]
        st.success(f"✅ {archivo_subido.name}")

    st.markdown("---")

    st.markdown("#### 🧮 **Tipo de Cálculo**")
    tipo_calculo = st.radio(
        "Selecciona el tipo de análisis:",
        ["Optimización de Geometría", "Frecuencias Vibracionales (IR)"],
        label_visibility="collapsed",
        help="Selecciona el tipo de análisis cuántico"
    )

    factor_escalamiento = 1.0
    if tipo_calculo == "Frecuencias Vibracionales (IR)":
        st.markdown("##### 📊 Factor de Escalamiento")
        factor_escalamiento = st.slider(
            "Factor IR",
            min_value=0.80, max_value=1.20, value=0.9679, step=0.001,
            help="Corrección para frecuencias calculadas"
        )

    st.markdown("---")

    st.markdown("#### ⚙️ **Configuración Computacional**")

    col_a, col_b = st.columns(2)
    with col_a:
        metodo = st.selectbox("Método", ["B3LYP", "PBE0", "M06-2X", "wB97X-D"])
    with col_b:
        conjunto_base = st.selectbox("Base", ["6-31+G(d,p)", "6-311++G(d,p)", "cc-pVDZ", "def2-SVP"])

    st.markdown("---")

    st.markdown("#### 🚀 **Ejecutar Cálculo**")
    boton_ejecutar = st.button(
        "🎯 **CALCULAR**",
        use_container_width=True,
        type="primary",
        help="Inicia el cálculo cuántico con PySCF"
    )

    if st.session_state.calculo_completado:
        if st.session_state.opt_convergida:
            st.success("✅ Cálculo completado")
        else:
            st.warning("⚠️ No convergió")

if boton_ejecutar:
    if st.session_state.xyz_inicial is None:
        st.sidebar.error("Por favor, carga un archivo .xyz primero.")
    else:
        # Resetear estados
        st.session_state.calculo_completado = False
        st.session_state.opt_convergida = False
        st.session_state.xyz_optimizada = None
        st.session_state.energia_final = None
        st.session_state.datos_energia = None
        st.session_state.datos_cargas = None
        st.session_state.datos_ir = None
        st.session_state.resultados_pyscf = None
        st.session_state.ultimo_tipo_calculo = tipo_calculo

        nombre_trabajo = st.session_state.nombre_trabajo

        with st.spinner(f"Ejecutando PySCF para '{nombre_trabajo}'... Esto puede tardar varios minutos."):
            try:
                # Parsear la geometría inicial
                geometria_inicial = parsear_xyz_contenido(st.session_state.xyz_inicial)

                # Configurar parámetros del cálculo
                config_calculo = {
                    'metodo': metodo,
                    'base': conjunto_base,
                    'tipo_calculo': tipo_calculo,
                    'factor_escalamiento': factor_escalamiento
                }

                # Ejecutar cálculo con PySCF
                resultados = ejecutar_calculo_pyscf(geometria_inicial, config_calculo, nombre_trabajo)

                if resultados is not None:
                    st.session_state.resultados_pyscf = resultados
                    st.session_state.calculo_completado = True
                    st.session_state.opt_convergida = resultados.get('convergido', False)

                    # Extraer datos para la interfaz
                    st.session_state.energia_final = extraer_energia_final(resultados)
                    st.session_state.xyz_optimizada = extraer_geometria_optimizada(resultados)
                    st.session_state.datos_energia = extraer_componentes_energia(resultados)
                    st.session_state.datos_cargas = extraer_cargas_atomicas(resultados)

                    if tipo_calculo == "Frecuencias Vibracionales (IR)":
                        st.session_state.datos_ir = extraer_espectro_ir(resultados, factor_escalamiento)
                else:
                    st.error("Error durante el cálculo con PySCF")

            except Exception as e:
                st.error(f"Ha ocurrido un error durante el cálculo: {str(e)}")
                st.exception(e)

        if st.session_state.calculo_completado:
            if st.session_state.opt_convergida:
                st.success(f"Cálculo para '{nombre_trabajo}' finalizado y convergido exitosamente.")
            else:
                st.warning(f"Cálculo para '{nombre_trabajo}' finalizado, pero la optimización NO convergió.")

        st.rerun()

# Métricas principales
if st.session_state.energia_final is not None:
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("🔋 Energía Final", f"{st.session_state.energia_final:.6f} Eh")
    with col2:
        estado = "✅ Convergido" if st.session_state.opt_convergida else "❌ No convergió"
        st.metric("📊 Estado", estado)
    with col3:
        if st.session_state.xyz_optimizada:
            num_atomos = len([l for l in st.session_state.xyz_optimizada.split('\n')[2:] if l.strip()])
            st.metric("⚛️ Átomos", f"{num_atomos}")
    with col4:
        st.metric("🧮 Método", f"{metodo}/{conjunto_base}")

# Pestañas principales
tabs = st.tabs(["🔬 **Visualización 3D**", "📈 **Espectroscopía**", "⚡ **Análisis Energético**", "🔧 **Datos Técnicos**"])

with tabs[0]:
    if not st.session_state.calculo_completado and st.session_state.xyz_inicial is None:
        st.info("💡 Carga un archivo .xyz en la barra lateral para empezar.")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("### 🧪 **Geometría Inicial**")
        if st.session_state.xyz_inicial:
            vista = py3Dmol.view(width=400, height=400)
            vista.addModel(st.session_state.xyz_inicial, 'xyz')
            vista.setStyle({'stick': {'radius': 0.15}, 'sphere': {'radius': 0.3}})
            vista.setBackgroundColor('white')
            vista.zoomTo()
            showmol(vista, height=450, width=450)
        else:
            st.info("Carga una molécula para visualizar")

    with col2:
        st.markdown("### 🎯 **Geometría Optimizada**")
        if st.session_state.xyz_optimizada:
            if not st.session_state.opt_convergida:
                st.warning("⚠️ Geometría no completamente optimizada")
            vista_opt = py3Dmol.view(width=400, height=400)
            vista_opt.addModel(st.session_state.xyz_optimizada, 'xyz')
            vista_opt.setStyle({'stick': {'radius': 0.15}, 'sphere': {'radius': 0.3}})
            vista_opt.setBackgroundColor('white')
            vista_opt.zoomTo()
            showmol(vista_opt, height=450, width=450)
        elif st.session_state.calculo_completado:
            st.error("❌ No se pudo extraer la geometría optimizada")
        else:
            st.info("Ejecuta un cálculo para ver la optimización")

with tabs[1]:
    if st.session_state.datos_ir is not None and not st.session_state.datos_ir.empty:
        st.markdown("### 📊 **Espectro Infrarrojo (IR)**")

        fig, ax = plt.subplots(figsize=(12, 6))
        ax.stem(st.session_state.datos_ir["Frequency"], st.session_state.datos_ir["Intensity"],
                basefmt=' ', linefmt='red', markerfmt='ro')
        ax.set_xlabel("Número de onda (cm⁻¹)", fontsize=12)
        ax.set_ylabel("Intensidad IR (km/mol)", fontsize=12)
        ax.set_title("Espectro IR Teórico", fontsize=14, fontweight='bold')
        ax.invert_xaxis()
        ax.grid(True, alpha=0.3)
        ax.set_facecolor('#fafafa')

        st.pyplot(fig)

        st.markdown("#### 📋 **Tabla de Frecuencias**")
        df_display = st.session_state.datos_ir.copy()
        df_display.columns = ['Frecuencia (cm⁻¹)', 'Intensidad (km/mol)']
        st.dataframe(
            df_display.style.format({"Frecuencia (cm⁻¹)": "{:.2f}", "Intensidad (km/mol)": "{:.2f}"}),
            use_container_width=True
        )

    elif st.session_state.ultimo_tipo_calculo == "Frecuencias Vibracionales (IR)":
        st.warning("⚠️ No se encontraron datos del espectro IR. Verifica que la optimización haya convergido.")
    else:
        st.info("💡 Selecciona 'Frecuencias Vibracionales (IR)' para ver el espectro.")

with tabs[2]:
    if not st.session_state.calculo_completado:
        st.info("💡 Ejecuta un cálculo para ver el análisis energético.")
    else:
        if st.session_state.datos_energia is not None:
            st.markdown("### ⚡ **Componentes Energéticos**")

            df_energia = st.session_state.datos_energia.copy()
            df_energia['Valor (eV)'] = df_energia['Energía (Hartree)'] * 27.2114  # Conversión a eV

            col1, col2 = st.columns([2, 1])
            with col1:
                st.dataframe(df_energia, use_container_width=True)
            with col2:
                st.markdown("##### 📊 **Conversiones**")
                st.markdown("- **1 Hartree** = 27.21 eV")
                st.markdown("- **1 Hartree** = 627.5 kcal/mol")

        if st.session_state.datos_cargas:
            st.markdown("---")
            st.markdown("### ⚛️ **Análisis de Cargas Atómicas**")

            col1, col2 = st.columns(2)
            with col1:
                if 'Mulliken' in st.session_state.datos_cargas and not st.session_state.datos_cargas['Mulliken'].empty:
                    st.markdown("#### 🔵 **Cargas de Mulliken**")
                    st.dataframe(st.session_state.datos_cargas['Mulliken'], use_container_width=True)

            with col2:
                if 'Loewdin' in st.session_state.datos_cargas and not st.session_state.datos_cargas['Loewdin'].empty:
                    st.markdown("#### 🟢 **Cargas de Loewdin**")
                    st.dataframe(st.session_state.datos_cargas['Loewdin'], use_container_width=True)

with tabs[3]:
    if not st.session_state.calculo_completado:
        st.info("💡 Ejecuta un cálculo para ver los datos técnicos.")
    else:
        col1, col2 = st.columns([2, 1])

        with col1:
            st.markdown("### 📋 **Información del Cálculo PySCF**")
            if st.session_state.resultados_pyscf:
                info_texto = f"""
**Energía SCF:** {st.session_state.resultados_pyscf.get('energia_scf', 'N/A'):.6f} Hartree
**Convergencia SCF:** {'✅ Sí' if st.session_state.resultados_pyscf.get('convergencia_scf', False) else '❌ No'}
**Iteraciones SCF:** {st.session_state.resultados_pyscf.get('iteraciones_scf', 'N/A')}
"""
                if 'dipole_moment' in st.session_state.resultados_pyscf:
                    dipole = st.session_state.resultados_pyscf['dipole_moment']
                    dipole_mag = np.linalg.norm(dipole)
                    info_texto += f"**Momento Dipolar:** {dipole_mag:.4f} Debye\n"

                st.markdown(info_texto)

                # Botón para descargar datos
                if st.session_state.resultados_pyscf:
                    datos_exportar = str(st.session_state.resultados_pyscf)
                    st.download_button(
                        label="💾 Descargar resultados",
                        data=datos_exportar,
                        file_name=f"{st.session_state.nombre_trabajo}_pyscf.txt",
                        mime="text/plain"
                    )

        with col2:
            st.markdown("### ℹ️ **Información del Cálculo**")
            info_data = {
                "Archivo": st.session_state.nombre_trabajo,
                "Método": metodo,
                "Base": conjunto_base,
                "Motor": "PySCF",
                "Tipo": st.session_state.ultimo_tipo_calculo,
                "Estado": "Convergido ✅" if st.session_state.opt_convergida else "No convergió ❌"
            }

            for key, value in info_data.items():
                st.text(f"{key}: {value}")

# Pie de página
st.markdown("---")
st.markdown("*Desarrollado con Streamlit • Cálculos cuánticos con PySCF*")