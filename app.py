import streamlit as st
import os
import subprocess
import pandas as pd
import numpy as np
import plotly.graph_objs as go
import py3Dmol
from stmol import showmol
import matplotlib.pyplot as plt
from utils import (
    generar_entrada_orca,
    extraer_espectro_ir,
    extraer_geometria_optimizada,
    extraer_energia_final,
    extraer_componentes_energia,
    extraer_cargas_atomicas,
    verificar_convergencia_optimizacion,
    extraer_energias_orbitales,
    extraer_cargas_orbitales_reducidas
)

# --- Configuración de la Página ---
st.set_page_config(
    page_title="ORCA Molecular",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- Título de la Aplicación ---
col1, col2 = st.columns([3, 1])
with col1:
    st.title("🧬 ORCA Molecular Calculator")
    st.markdown("*Calculadora cuántica para análisis molecular*")

# --- Inicialización del Estado de la Sesión ---
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
if "datos_orbitales" not in st.session_state:
    st.session_state.datos_orbitales = None
if "datos_cargas_reducidas" not in st.session_state:
    st.session_state.datos_cargas_reducidas = None
if "resumen_log_orca" not in st.session_state:
    st.session_state.resumen_log_orca = None
if "log_completo_orca" not in st.session_state:
    st.session_state.log_completo_orca = None
if "nombre_trabajo" not in st.session_state:
    st.session_state.nombre_trabajo = ""

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
    with st.expander("Parámetros Avanzados", expanded=False):
        col_a, col_b = st.columns(2)
        with col_a:
            metodo = st.selectbox("Método", ["B3LYP", "PBE0", "M06-2X", "wB97X-D"])
        with col_b:
            conjunto_base = st.selectbox("Base", ["6-31+G(d,p)", "6-311++G(d,p)", "cc-pVDZ", "def2-SVP"])

        palabras_clave = st.text_input("Palabras clave extra", "D3BJ TIGHTSCF")

    st.markdown("---")

    st.markdown("#### 🚀 **Ejecutar Cálculo**")
    boton_ejecutar = st.button(
        "🎯 **CALCULAR**",
        use_container_width=True,
        type="primary",
        help="Inicia el cálculo cuántico con ORCA"
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
        for key in st.session_state.keys():
            if key not in ['xyz_inicial', 'nombre_trabajo', 'frames']:
                st.session_state[key] = None if not isinstance(st.session_state[key], bool) else False

        st.session_state.ultimo_tipo_calculo = tipo_calculo

        nombre_trabajo = st.session_state.nombre_trabajo
        contenido_entrada = generar_entrada_orca(st.session_state.xyz_inicial, tipo_calculo, metodo, conjunto_base,
                                                 palabras_clave)
        ruta_entrada = os.path.join(DIR_CALCULOS, f"{nombre_trabajo}.inp")
        ruta_salida = os.path.join(DIR_CALCULOS, f"{nombre_trabajo}.out")

        with open(ruta_entrada, "w") as f:
            f.write(contenido_entrada)

        with st.spinner(f"Ejecutando ORCA para '{nombre_trabajo}'... Esto puede tardar varios minutos."):
            try:
                comando = f"orca {ruta_entrada} > {ruta_salida}"
                subprocess.run(comando, shell=True, check=True, timeout=600)
                st.session_state.calculo_completado = True
            except Exception as e:
                st.error(f"Error al ejecutar ORCA: {e}")

        if os.path.exists(ruta_salida):
            with open(ruta_salida, 'r', encoding='utf-8', errors='ignore') as f:
                lineas = f.readlines()
                st.session_state.resumen_log_orca = "".join(lineas[-50:])
                st.session_state.log_completo_orca = "".join(lineas)

            st.session_state.opt_convergida = verificar_convergencia_optimizacion(ruta_salida)
            st.session_state.xyz_optimizada = extraer_geometria_optimizada(ruta_salida)
            st.session_state.energia_final = extraer_energia_final(ruta_salida)
            st.session_state.datos_energia = extraer_componentes_energia(ruta_salida)
            st.session_state.datos_cargas = extraer_cargas_atomicas(ruta_salida)
            st.session_state.datos_orbitales = extraer_energias_orbitales(ruta_salida)
            st.session_state.datos_cargas_reducidas = extraer_cargas_orbitales_reducidas(ruta_salida)

            if tipo_calculo == "Frecuencias Vibracionales (IR)":
                st.session_state.datos_ir = extraer_espectro_ir(ruta_salida, factor_escalamiento)

        st.rerun()

if st.session_state.energia_final is not None:
    st.markdown("---")
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
    st.markdown("---")

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
            vista.setBackgroundColor('#F7F7F7')
            vista.zoomTo()
            showmol(vista, height=450, width=450)

    with col2:
        st.markdown("### 🎯 **Geometría Optimizada**")
        if st.session_state.xyz_optimizada:
            if not st.session_state.opt_convergida:
                st.warning("⚠️ Geometría no completamente optimizada.")
            vista_opt = py3Dmol.view(width=400, height=400)
            vista_opt.addModel(st.session_state.xyz_optimizada, 'xyz')
            vista_opt.setStyle({'stick': {'radius': 0.15}, 'sphere': {'radius': 0.3}})
            vista_opt.setBackgroundColor('#F7F7F7')
            vista_opt.zoomTo()
            showmol(vista_opt, height=450, width=450)

with tabs[1]:
    if st.session_state.datos_ir is not None and not st.session_state.datos_ir.empty:
        st.markdown("### 📊 **Espectro Infrarrojo (IR)**")
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.stem(st.session_state.datos_ir["Frequency"], st.session_state.datos_ir["Intensity"], basefmt=' ',
                linefmt='red', markerfmt='ro')
        ax.set_xlabel("Número de onda (cm⁻¹)");
        ax.set_ylabel("Intensidad IR (km/mol)");
        ax.set_title("Espectro IR Teórico")
        ax.invert_xaxis();
        ax.grid(True, alpha=0.3)
        st.pyplot(fig)
        st.dataframe(st.session_state.datos_ir.style.format({"Frequency": "{:.2f}", "Intensity": "{:.2f}"}),
                     use_container_width=True)
    elif st.session_state.ultimo_tipo_calculo == "Frecuencias Vibracionales (IR)":
        st.warning(
            "⚠️ No se encontraron datos IR. Verifica que la optimización haya convergido en la pestaña de Datos Técnicos.")
    else:
        st.info("💡 Selecciona 'Frecuencias Vibracionales (IR)' y ejecuta un cálculo para ver el espectro.")

with tabs[2]:
    if not st.session_state.calculo_completado:
        st.info("💡 Ejecuta un cálculo para ver el análisis detallado.")
    else:
        st.markdown("### ⚡ **Componentes Energéticos**")
        if st.session_state.datos_energia is not None:
            st.dataframe(st.session_state.datos_energia, use_container_width=True)

        st.markdown("### 🔋 **Energías Orbitales**")
        if st.session_state.datos_orbitales is not None:
            st.dataframe(st.session_state.datos_orbitales, use_container_width=True)

        st.markdown("### ⚛️ **Análisis de Cargas**")
        col1, col2 = st.columns(2)
        with col1:
            st.markdown("#### Cargas Atómicas")
            if st.session_state.datos_cargas:
                if 'Mulliken' in st.session_state.datos_cargas and not st.session_state.datos_cargas['Mulliken'].empty:
                    st.write("**Cargas de Mulliken**");
                    st.dataframe(st.session_state.datos_cargas['Mulliken'])
                if 'Loewdin' in st.session_state.datos_cargas and not st.session_state.datos_cargas['Loewdin'].empty:
                    st.write("**Cargas de Loewdin**");
                    st.dataframe(st.session_state.datos_cargas['Loewdin'])
        with col2:
            st.markdown("#### Cargas Orbitales Reducidas")
            if st.session_state.datos_cargas_reducidas:
                if 'Mulliken' in st.session_state.datos_cargas_reducidas and not \
                st.session_state.datos_cargas_reducidas['Mulliken'].empty:
                    st.write("**Mulliken (Reducidas)**");
                    st.dataframe(st.session_state.datos_cargas_reducidas['Mulliken'])
                if 'Loewdin' in st.session_state.datos_cargas_reducidas and not st.session_state.datos_cargas_reducidas[
                    'Loewdin'].empty:
                    st.write("**Loewdin (Reducidas)**");
                    st.dataframe(st.session_state.datos_cargas_reducidas['Loewdin'])

with tabs[3]:
    if not st.session_state.calculo_completado:
        st.info("💡 Ejecuta un cálculo para ver los datos técnicos.")
    else:
        st.markdown("### 📋 **Log de Salida de ORCA**")
        if st.session_state.log_completo_orca:
            st.download_button(label="💾 Descargar Archivo .out Completo", data=st.session_state.log_completo_orca,
                               file_name=f"{st.session_state.nombre_trabajo}.out")
            with st.expander("📄 Ver Log Completo"):
                st.code(st.session_state.log_completo_orca, language='text')

st.markdown("---")
st.markdown("*Desarrollado con Streamlit • Cálculos cuánticos con ORCA*")

