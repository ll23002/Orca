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
    verificar_convergencia_optimizacion
)

st.set_page_config(page_title="Visualizador ORCA", layout="wide")
st.title("Interfaz Gráfica para ORCA 🧪")

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
if "resumen_log_orca" not in st.session_state:
    st.session_state.resumen_log_orca = None
if "log_completo_orca" not in st.session_state:
    st.session_state.log_completo_orca = None
if "nombre_trabajo" not in st.session_state:
    st.session_state.nombre_trabajo = ""

DIR_CALCULOS = "calculations"
os.makedirs(DIR_CALCULOS, exist_ok=True)

with st.sidebar:
    st.header("Configuración del Cálculo")

    archivo_subido = st.file_uploader("Cargar Molécula (.xyz)", type=["xyz"])

    if archivo_subido is not None:
        st.session_state.xyz_inicial = archivo_subido.getvalue().decode("utf-8")
        st.session_state.nombre_trabajo = os.path.splitext(archivo_subido.name)[0]

    tipo_calculo = st.selectbox(
        "Tipo de Cálculo",
        ["Optimización de Geometría", "Frecuencias Vibracionales (IR)"],
        key="selector_tipo_calculo"
    )

    with st.expander("Configuración Avanzada"):
        metodo = st.text_input("Método/Funcional", "B3LYP")
        conjunto_base = st.text_input("Conjunto Base", "6-31+G(d,p)")
        palabras_clave = st.text_input("Palabras Clave Adicionales", "D3BJ TIGHTSCF")

    factor_escalamiento = 1.0
    if tipo_calculo == "Frecuencias Vibracionales (IR)":
        factor_escalamiento = st.number_input(
            "Factor de Escalamiento IR",
            min_value=0.8, max_value=1.2, value=0.9679, step=0.001,
            help="Factor para corregir las frecuencias calculadas. Valor común para B3LYP/6-31G(d): ~0.96"
        )

    boton_ejecutar = st.button("▶️ Run", use_container_width=True)

if boton_ejecutar:
    if st.session_state.xyz_inicial is None:
        st.sidebar.error("Por favor, carga un archivo .xyz primero.")
    else:
        st.session_state.calculo_completado = False
        st.session_state.opt_convergida = False
        st.session_state.xyz_optimizada = None
        st.session_state.energia_final = None
        st.session_state.datos_energia = None
        st.session_state.datos_cargas = None
        st.session_state.datos_ir = None
        st.session_state.resumen_log_orca = None
        st.session_state.log_completo_orca = None
        st.session_state.ultimo_tipo_calculo = tipo_calculo

        nombre_trabajo = st.session_state.nombre_trabajo
        contenido_entrada = generar_entrada_orca(st.session_state.xyz_inicial, tipo_calculo, metodo, conjunto_base,
                                                 palabras_clave)
        ruta_entrada = os.path.join(DIR_CALCULOS, f"{nombre_trabajo}.inp")
        ruta_salida = os.path.join(DIR_CALCULOS, f"{nombre_trabajo}.out")

        with open(ruta_entrada, "w") as f:
            f.write(contenido_entrada)

        # Ejecutar ORCA
        with st.spinner(f"Ejecutando ORCA para '{nombre_trabajo}'... Esto puede tardar varios minutos."):
            try:
                comando = f"orca {ruta_entrada} > {ruta_salida}"
                subprocess.run(comando, shell=True, check=True, timeout=600)
                st.session_state.calculo_completado = True

            except FileNotFoundError:
                st.error(
                    "Error: El ejecutable 'orca' no se encontró. Asegúrate de que ORCA esté instalado y en el PATH del sistema.")
            except subprocess.CalledProcessError:
                st.error(f"ORCA devolvió un error. Revisa el archivo de salida para más detalles: {ruta_salida}")
            except subprocess.TimeoutExpired:
                st.error("El cálculo de ORCA ha superado el tiempo límite de 10 minutos y ha sido detenido.")
            except Exception as e:
                st.error(f"Ha ocurrido un error inesperado: {e}")

        if os.path.exists(ruta_salida):
            with open(ruta_salida, 'r') as f:
                lineas = f.readlines()
                st.session_state.resumen_log_orca = "".join(lineas[-50:])
                st.session_state.log_completo_orca = "".join(lineas)

            st.session_state.opt_convergida = verificar_convergencia_optimizacion(ruta_salida)
            st.session_state.xyz_optimizada = extraer_geometria_optimizada(ruta_salida)
            st.session_state.energia_final = extraer_energia_final(ruta_salida)
            st.session_state.datos_energia = extraer_componentes_energia(ruta_salida)
            st.session_state.datos_cargas = extraer_cargas_atomicas(ruta_salida)
            if tipo_calculo == "Frecuencias Vibracionales (IR)":
                st.session_state.datos_ir = extraer_espectro_ir(ruta_salida, factor_escalamiento)
        else:
            st.session_state.resumen_log_orca = "El archivo de salida no fue creado."

        if st.session_state.calculo_completado:
            if st.session_state.opt_convergida:
                st.success(f"Cálculo para '{nombre_trabajo}' finalizado y convergido exitosamente.")
            else:
                st.warning(f"Cálculo para '{nombre_trabajo}' finalizado, pero la optimización NO convergió.")

        st.rerun()

# Pestañas
tab1, tab2, tab3 = st.tabs(
    ["Resultados del Cálculo", "Análisis Detallado", "Demostración Visual"]
)

with tab1:
    st.header("Visualización Molecular")
    if not st.session_state.calculo_completado and st.session_state.xyz_inicial is None:
        st.info("Carga un archivo .xyz en la barra lateral para empezar.")

    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Molécula Inicial")
        if st.session_state.xyz_inicial:
            vista = py3Dmol.view(width=400, height=400)
            vista.addModel(st.session_state.xyz_inicial, 'xyz')
            vista.setStyle({'stick': {}})
            vista.zoomTo()
            showmol(vista, height=500, width=500)
        else:
            st.write("No se ha cargado ninguna molécula.")

    with col2:
        st.subheader("Geometría Optimizada")
        if st.session_state.xyz_optimizada:
            if not st.session_state.opt_convergida:
                st.warning("Esta es la última geometría antes de que el cálculo fallara.")
            vista_opt = py3Dmol.view(width=400, height=400)
            vista_opt.addModel(st.session_state.xyz_optimizada, 'xyz')
            vista_opt.setStyle({'stick': {}})
            vista_opt.zoomTo()
            showmol(vista_opt, height=500, width=500)
        elif st.session_state.calculo_completado:
            st.warning("No se pudo extraer una geometría optimizada del archivo de salida.")
        else:
            st.write("Ejecuta un cálculo para ver el resultado.")

    if st.session_state.energia_final is not None:
        st.metric("Energía Total Final (Hartree)", f"{st.session_state.energia_final:.6f}")

    if st.session_state.resumen_log_orca:
        with st.expander("Ver Resumen del Log de ORCA"):
            st.code(st.session_state.resumen_log_orca, language='text')

with tab2:
    st.header("Análisis Detallado")

    if not st.session_state.calculo_completado:
        st.info("Ejecuta un cálculo para ver los resultados detallados.")
    else:
        # Herramientas de análisis
        st.subheader("Descarga y Análisis")
        if st.session_state.log_completo_orca:
            st.download_button(
                label="📥 Descargar Log Completo (.out)",
                data=st.session_state.log_completo_orca,
                file_name=f"{st.session_state.nombre_trabajo}.out",
                mime="text/plain"
            )
            with st.expander("Ver Log Completo en el Navegador"):
                st.code(st.session_state.log_completo_orca, language='text')

        if not st.session_state.opt_convergida:
            st.error(
                "¡Atención! La optimización de la geometría no convergió. Los resultados pueden no ser confiables.")

        st.markdown("---")

        # Componentes de Energía
        if st.session_state.datos_energia is not None:
            st.subheader("Componentes de Energía")
            st.dataframe(st.session_state.datos_energia)

        # Cargas Atómicas
        if st.session_state.datos_cargas:
            st.subheader("Cargas Atómicas Parciales")
            if 'Mulliken' in st.session_state.datos_cargas and not st.session_state.datos_cargas['Mulliken'].empty:
                st.write("**Cargas de Mulliken**")
                st.dataframe(st.session_state.datos_cargas['Mulliken'])
            if 'Loewdin' in st.session_state.datos_cargas and not st.session_state.datos_cargas['Loewdin'].empty:
                st.write("**Cargas de Loewdin**")
                st.dataframe(st.session_state.datos_cargas['Loewdin'])

        # Espectro IR
        if st.session_state.datos_ir is not None and not st.session_state.datos_ir.empty:
            st.subheader("Espectro Infrarrojo (IR)")
            fig, ax = plt.subplots()
            ax.stem(st.session_state.datos_ir["Frequency"], st.session_state.datos_ir["Intensity"])
            ax.set_xlabel("Número de onda (cm⁻¹)")
            ax.set_ylabel("Intensidad (km/mol)")
            ax.set_title("Espectro IR Teórico (Escalado)")
            ax.invert_xaxis()
            st.pyplot(fig)
            with st.expander("Ver Tabla de Frecuencias"):
                st.dataframe(st.session_state.datos_ir.style.format({"Frequency": "{:.2f}", "Intensity": "{:.2f}"}))
        elif st.session_state.ultimo_tipo_calculo == "Frecuencias Vibracionales (IR)":
            st.warning("No se encontraron datos del espectro IR. Esto es normal si la optimización no convergió.")

with tab3:
    st.header("Demostración Visual")
    st.info(
        "ℹ️ Esta pestaña contiene una **demostración** con datos sintéticos para mostrar capacidades de visualización.")

    # Simulación Dinámica
    st.subheader("Simulación Molecular Dinámica")

    n_pasos, n_polimero, n_agua = 50, 15, 30
    np.random.seed(42)
    if "cuadros" not in st.session_state:
        cuadros = []
        for paso in range(n_pasos):
            cadena_polimero = np.cumsum(np.random.randn(n_polimero, 3) * 0.2, axis=0)
            aguas = np.random.rand(n_agua, 3) * 6 - 3
            cuadros.append((cadena_polimero.copy(), aguas.copy()))
        st.session_state.cuadros = cuadros

    cuadro_actual = st.slider("Cuadro de Animación", 0, n_pasos - 1, 0, key="slider_cuadro_demo")

    polimero, aguas = st.session_state.cuadros[cuadro_actual]
    traza_polimero = go.Scatter3d(x=polimero[:, 0], y=polimero[:, 1], z=polimero[:, 2], mode='markers+lines',
                                  marker=dict(size=8, color='blue'), line=dict(color='darkblue', width=4),
                                  name='Polímero (Ejemplo)')
    traza_agua = go.Scatter3d(x=aguas[:, 0], y=aguas[:, 1], z=aguas[:, 2], mode='markers',
                              marker=dict(size=6, color='red'), name='Agua (Ejemplo)')
    layout = go.Layout(
        title=f'Cuadro: {cuadro_actual}',
        scene=dict(xaxis=dict(title='X', range=[-4, 4]), yaxis=dict(title='Y', range=[-4, 4]),
                   zaxis=dict(title='Z', range=[-4, 4]), aspectmode='cube'),
        margin=dict(l=0, r=0, b=0, t=40),
        height=500
    )
    fig_3d = go.Figure(data=[traza_polimero, traza_agua], layout=layout)
    st.plotly_chart(fig_3d, use_container_width=True)