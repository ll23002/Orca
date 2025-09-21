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
    generate_orca_input,
    parse_ir_spectrum,
    parse_optimized_geometry,
    parse_final_energy,
    parse_energy_components,
    parse_atomic_charges,
    check_optimization_convergence,  # <-- Nueva función importada
    debug_file_sections  # <-- Agregada función de debug
)
import time

# --- Configuración de la Página ---
st.set_page_config(page_title="ORCA Visualizer", layout="wide")
st.title("Interfaz Gráfica para ORCA 🧪")

# --- Inicialización del Estado de la Sesión ---
if "calculation_done" not in st.session_state:
    st.session_state.calculation_done = False
if "opt_converged" not in st.session_state:
    st.session_state.opt_converged = False
if "last_calc_type" not in st.session_state:
    st.session_state.last_calc_type = None
if "initial_xyz" not in st.session_state:
    st.session_state.initial_xyz = None
if "optimized_xyz" not in st.session_state:
    st.session_state.optimized_xyz = None
if "final_energy" not in st.session_state:
    st.session_state.final_energy = None
if "energy_data" not in st.session_state:
    st.session_state.energy_data = None
if "charge_data" not in st.session_state:
    st.session_state.charge_data = None
if "ir_data" not in st.session_state:
    st.session_state.ir_data = None
if "orca_log_summary" not in st.session_state:
    st.session_state.orca_log_summary = None
if "full_orca_log" not in st.session_state:
    st.session_state.full_orca_log = None
if "job_name" not in st.session_state:
    st.session_state.job_name = ""

# --- Directorios ---
CALC_DIR = "calculations"
os.makedirs(CALC_DIR, exist_ok=True)

# --- Barra Lateral de Controles ---
with st.sidebar:
    st.header("Configuración del Cálculo")

    uploaded_file = st.file_uploader("Cargar Molécula (.xyz)", type=["xyz"])

    if uploaded_file is not None:
        st.session_state.initial_xyz = uploaded_file.getvalue().decode("utf-8")
        st.session_state.job_name = os.path.splitext(uploaded_file.name)[0]

    calc_type = st.selectbox(
        "Tipo de Cálculo",
        ["Optimización de Geometría", "Frecuencias Vibracionales (IR)"],
        key="calc_type_selector"
    )

    with st.expander("Parámetros Avanzados"):
        method = st.text_input("Método/Funcional", "B3LYP")
        basis_set = st.text_input("Conjunto Base", "6-31+G(d,p)")
        keywords = st.text_input("Palabras Clave Adicionales", "D3BJ TIGHTSCF")

    scaling_factor = 1.0
    if calc_type == "Frecuencias Vibracionales (IR)":
        scaling_factor = st.number_input(
            "Factor de Escalamiento IR",
            min_value=0.8, max_value=1.2, value=0.9679, step=0.001,
            help="Factor para corregir las frecuencias calculadas. Un valor común para B3LYP/6-31G(d) es ~0.96"
        )

    run_button = st.button("▶️ Iniciar Cálculo", use_container_width=True)

# --- Lógica Principal del Cálculo ---
if run_button:
    if st.session_state.initial_xyz is None:
        st.sidebar.error("Por favor, carga un archivo .xyz primero.")
    else:
        # Resetear los resultados de cálculos anteriores
        st.session_state.calculation_done = False
        st.session_state.opt_converged = False
        st.session_state.optimized_xyz = None
        st.session_state.final_energy = None
        st.session_state.energy_data = None
        st.session_state.charge_data = None
        st.session_state.ir_data = None
        st.session_state.orca_log_summary = None
        st.session_state.full_orca_log = None
        st.session_state.last_calc_type = calc_type

        # 1. Preparar archivos de entrada y salida
        job_name = st.session_state.job_name
        input_content = generate_orca_input(st.session_state.initial_xyz, calc_type, method, basis_set, keywords)
        input_path = os.path.join(CALC_DIR, f"{job_name}.inp")
        output_path = os.path.join(CALC_DIR, f"{job_name}.out")

        with open(input_path, "w") as f:
            f.write(input_content)

        # 2. Ejecutar ORCA
        with st.spinner(f"Ejecutando ORCA para '{job_name}'... Esto puede tardar varios minutos."):
            try:
                command = f"orca {input_path} > {output_path}"
                subprocess.run(command, shell=True, check=True, timeout=600)
                st.session_state.calculation_done = True

            except FileNotFoundError:
                st.error(
                    "Error: El ejecutable 'orca' no se encontró. Asegúrate de que ORCA esté instalado y añadido al PATH del sistema.")
            except subprocess.CalledProcessError as e:
                st.error(f"ORCA devolvió un error. Revisa el archivo de salida para más detalles: {output_path}")
            except subprocess.TimeoutExpired:
                st.error("El cálculo de ORCA ha superado el tiempo límite de 10 minutos y ha sido detenido.")
            except Exception as e:
                st.error(f"Ha ocurrido un error inesperado: {e}")

        # 3. Leer y Parsear los resultados
        if os.path.exists(output_path):
            with open(output_path, 'r') as f:
                lines = f.readlines()
                st.session_state.orca_log_summary = "".join(lines[-50:])
                st.session_state.full_orca_log = "".join(lines)

            # ===== AQUÍ SE AGREGA EL DEBUG =====
            st.write("🔍 **Ejecutando análisis de debug del archivo de salida...**")
            debug_file_sections(output_path)  # Esto te mostrará qué hay en tu archivo
            # ===================================

            # Chequear la convergencia y parsear los datos
            st.session_state.opt_converged = check_optimization_convergence(output_path)
            st.session_state.optimized_xyz = parse_optimized_geometry(output_path)
            st.session_state.final_energy = parse_final_energy(output_path)
            st.session_state.energy_data = parse_energy_components(output_path)
            st.session_state.charge_data = parse_atomic_charges(output_path)
            if calc_type == "Frecuencias Vibracionales (IR)":
                st.session_state.ir_data = parse_ir_spectrum(output_path, scaling_factor)
        else:
            st.session_state.orca_log_summary = "El archivo de salida no fue creado."

        if st.session_state.calculation_done:
            if st.session_state.opt_converged:
                st.success(f"Cálculo para '{job_name}' finalizado y convergido exitosamente.")
            else:
                st.warning(f"Cálculo para '{job_name}' finalizado, pero la optimización NO convergió.")

        st.rerun()

# --- Área Principal de Visualización con Pestañas ---
tab1, tab2, tab3 = st.tabs(
    ["Resultados del Cálculo (Real)", "Análisis Detallado (Real)", "Prototipos Visuales"]
)

with tab1:
    st.header("Visualización Molecular")
    if not st.session_state.calculation_done and st.session_state.initial_xyz is None:
        st.info("Carga un archivo .xyz en la barra lateral para empezar.")

    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Molécula Inicial")
        if st.session_state.initial_xyz:
            view = py3Dmol.view(width=400, height=400)
            view.addModel(st.session_state.initial_xyz, 'xyz')
            view.setStyle({'stick': {}})
            view.zoomTo()
            showmol(view, height=500, width=500)
        else:
            st.write("No se ha cargado ninguna molécula.")

    with col2:
        st.subheader("Geometría Final")
        if st.session_state.optimized_xyz:
            if not st.session_state.opt_converged:
                st.warning("Esta es la última geometría antes de que el cálculo fallara.")
            view_opt = py3Dmol.view(width=400, height=400)
            view_opt.addModel(st.session_state.optimized_xyz, 'xyz')
            view_opt.setStyle({'stick': {}})
            view_opt.zoomTo()
            showmol(view_opt, height=500, width=500)
        elif st.session_state.calculation_done:
            st.warning("No se pudo extraer una geometría final del archivo de salida.")
        else:
            st.write("Ejecuta un cálculo para ver el resultado.")

    if st.session_state.final_energy is not None:
        st.metric("Energía Total Final (Hartree)", f"{st.session_state.final_energy:.6f}")

    if st.session_state.orca_log_summary:
        with st.expander("Ver Resumen del Log de ORCA (últimas 50 líneas)"):
            st.code(st.session_state.orca_log_summary, language='text')

with tab2:
    st.header("Análisis Detallado y Depuración")

    if not st.session_state.calculation_done:
        st.info("Ejecuta un cálculo para ver los resultados detallados y las herramientas de depuración.")
    else:
        # Herramientas de Depuración
        st.subheader("Herramientas de Depuración")
        if st.session_state.full_orca_log:
            st.download_button(
                label="📥 Descargar Log Completo (.out)",
                data=st.session_state.full_orca_log,
                file_name=f"{st.session_state.job_name}.out",
                mime="text/plain"
            )
            with st.expander("Ver Log Completo en el Navegador"):
                st.code(st.session_state.full_orca_log, language='text')

        # Advertencia de Convergencia
        if not st.session_state.opt_converged:
            st.error(
                "¡Atención! La optimización de la geometría no convergió. Los resultados numéricos a continuación pueden no ser significativos. Revisa el log completo para encontrar la causa del error.")

        st.markdown("---")

        # Desglose de Energía
        if st.session_state.energy_data is not None:
            st.subheader("Desglose de Componentes de Energía")
            st.dataframe(st.session_state.energy_data)

        # Cargas Atómicas
        if st.session_state.charge_data:
            st.subheader("Cargas Atómicas Parciales")
            if 'Mulliken' in st.session_state.charge_data and not st.session_state.charge_data['Mulliken'].empty:
                st.write("**Cargas de Mulliken**")
                st.dataframe(st.session_state.charge_data['Mulliken'])
            if 'Loewdin' in st.session_state.charge_data and not st.session_state.charge_data['Loewdin'].empty:
                st.write("**Cargas de Loewdin**")
                st.dataframe(st.session_state.charge_data['Loewdin'])

        # Espectro IR
        if st.session_state.ir_data is not None and not st.session_state.ir_data.empty:
            st.subheader("Espectro Infrarrojo (IR)")
            fig, ax = plt.subplots()
            ax.stem(st.session_state.ir_data["Frequency"], st.session_state.ir_data["Intensity"])
            ax.set_xlabel("Número de onda (cm⁻¹)")
            ax.set_ylabel("Intensidad (km/mol)")
            ax.set_title("Espectro IR Teórico (Escalado)")
            ax.invert_xaxis()
            st.pyplot(fig)
            with st.expander("Ver Tabla de Frecuencias"):
                st.dataframe(st.session_state.ir_data.style.format({"Frequency": "{:.2f}", "Intensity": "{:.2f}"}))
        elif st.session_state.last_calc_type == "Frecuencias Vibracionales (IR)":
            st.warning(
                "No se encontraron datos del espectro IR. Esto es normal si la optimización de geometría no convergió.")

with tab3:
    st.header("Prototipos Visuales (Datos Sintéticos)")
    st.info(
        "ℹ️ El contenido de esta pestaña es una **demostración** con datos generados aleatoriamente. No utiliza los resultados de tu cálculo con ORCA.")

    # Simulación Dinámica (Sintética)
    st.subheader("Animación de Simulación Dinámica (Ejemplo)")

    n_steps, n_poly, n_water = 50, 15, 30
    np.random.seed(42)
    if "frames" not in st.session_state:
        frames = []
        for step in range(n_steps):
            poly_chain = np.cumsum(np.random.randn(n_poly, 3) * 0.2, axis=0)
            waters = np.random.rand(n_water, 3) * 6 - 3
            frames.append((poly_chain.copy(), waters.copy()))
        st.session_state.frames = frames

    current_frame = st.slider("Frame de Animación", 0, n_steps - 1, 0, key="frame_slider_demo")

    poly, waters = st.session_state.frames[current_frame]
    trace_poly = go.Scatter3d(x=poly[:, 0], y=poly[:, 1], z=poly[:, 2], mode='markers+lines',
                              marker=dict(size=8, color='blue'), line=dict(color='darkblue', width=4),
                              name='Polímero (Ejemplo)')
    trace_water = go.Scatter3d(x=waters[:, 0], y=waters[:, 1], z=waters[:, 2], mode='markers',
                               marker=dict(size=6, color='red'), name='Agua (Ejemplo)')
    layout = go.Layout(
        title=f'Frame: {current_frame}',
        scene=dict(xaxis=dict(title='X', range=[-4, 4]), yaxis=dict(title='Y', range=[-4, 4]),
                   zaxis=dict(title='Z', range=[-4, 4]), aspectmode='cube'),
        margin=dict(l=0, r=0, b=0, t=40),
        height=500
    )
    fig_3d = go.Figure(data=[trace_poly, trace_water], layout=layout)
    st.plotly_chart(fig_3d, use_container_width=True)