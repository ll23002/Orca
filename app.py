import streamlit as st
import os
import pandas as pd
import py3Dmol
from stmol import showmol
import matplotlib.pyplot as plt
# Importar la nueva función única de utils
from utils import ejecutar_y_procesar_orca

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
def inicializar_estado():
    estados_iniciales = {
        "calculo_completado": False,
        "opt_convergida": False,
        "ultimo_tipo_calculo": None,
        "xyz_inicial": None,
        "xyz_optimizada": None,
        "energia_final": None,
        "datos_energia": None,
        "datos_cargas": None,
        "datos_ir": None,
        "datos_orbitales": None,
        "log_completo_orca": None,
        "nombre_trabajo": "",
        "calculo_en_progreso": False
    }

    for key, default_value in estados_iniciales.items():
        if key not in st.session_state:
            st.session_state[key] = default_value


inicializar_estado()

DIR_CALCULOS = "calculations"
os.makedirs(DIR_CALCULOS, exist_ok=True)

# --- Barra Lateral ---
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
    # Deshabilitar botón si hay cálculo en progreso
    boton_ejecutar = st.button(
        "🎯 **CALCULAR**" if not st.session_state.calculo_en_progreso else "⏳ Calculando...",
        use_container_width=True,
        type="primary",
        help="Inicia el cálculo cuántico con ORCA",
        disabled=st.session_state.calculo_en_progreso
    )

    if st.session_state.calculo_completado:
        if st.session_state.opt_convergida:
            st.success("✅ Cálculo completado")
        else:
            st.warning("⚠️ No convergió")

# --- Lógica de Ejecución del Cálculo (CORREGIDA) ---
if boton_ejecutar and not st.session_state.calculo_en_progreso:
    if st.session_state.xyz_inicial is None:
        st.sidebar.error("Por favor, carga un archivo .xyz primero.")
    else:
        # Marcar que el cálculo está en progreso
        st.session_state.calculo_en_progreso = True

        # Limpiar solo los resultados anteriores, NO los datos de entrada
        resultados_a_limpiar = [
            "calculo_completado", "opt_convergida", "xyz_optimizada",
            "energia_final", "datos_energia", "datos_cargas",
            "datos_ir", "datos_orbitales", "log_completo_orca"
        ]

        for key in resultados_a_limpiar:
            if key in st.session_state:
                if isinstance(st.session_state[key], bool):
                    st.session_state[key] = False
                else:
                    st.session_state[key] = None

        st.session_state.ultimo_tipo_calculo = tipo_calculo
        nombre_trabajo = st.session_state.nombre_trabajo

        # Mostrar progreso
        with st.spinner(f"🔄 Ejecutando cálculo {tipo_calculo} para {nombre_trabajo}..."):
            st.info("🔹 Iniciando cálculo ORCA...")

            # Llamar a la función principal
            resultados = ejecutar_y_procesar_orca(
                nombre_trabajo=nombre_trabajo,
                contenido_xyz=st.session_state.xyz_inicial,
                tipo_calculo=tipo_calculo,
                metodo=metodo,
                base=conjunto_base,
                palabras_clave_extra=palabras_clave,
                dir_calculos=DIR_CALCULOS,
                factor_escalamiento=factor_escalamiento
            )

        # Marcar que el cálculo ya no está en progreso
        st.session_state.calculo_en_progreso = False

        # Procesar resultados
        if resultados and "error" in resultados and resultados["error"]:
            st.error(f"❌ Error en el cálculo: {resultados['error']}")
            # Aún así guardamos el log si está disponible
            if "log_completo_orca" in resultados:
                st.session_state.log_completo_orca = resultados["log_completo_orca"]
        else:
            st.success("✅ Cálculo finalizado exitosamente")
            st.session_state.calculo_completado = True

            # Actualizar todos los resultados en el estado de la sesión
            if resultados:
                for key, value in resultados.items():
                    if key != "error":  # No guardar el error si no hay error
                        st.session_state[key] = value

                # Debug: Mostrar qué se guardó
                st.info(
                    f"🔹 Resultados guardados: Energía={st.session_state.energia_final}, Convergió={st.session_state.opt_convergida}")

        # Forzar actualización de la interfaz
        st.rerun()

# --- Visualización de Resultados ---
# Verificar si tenemos resultados para mostrar
mostrar_metricas = (st.session_state.energia_final is not None or
                    st.session_state.calculo_completado)

if mostrar_metricas:
    st.markdown("---")
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        energia_str = f"{st.session_state.energia_final:.6f} Eh" if st.session_state.energia_final else "No disponible"
        st.metric("🔋 Energía Final", energia_str)
    with col2:
        if st.session_state.opt_convergida:
            estado = "✅ Convergido"
        elif st.session_state.calculo_completado:
            estado = "❌ No convergió"
        else:
            estado = "⏳ En progreso"
        st.metric("📊 Estado", estado)
    with col3:
        if st.session_state.xyz_optimizada:
            try:
                num_atomos = len([l for l in st.session_state.xyz_optimizada.split('\n')[2:] if l.strip()])
                st.metric("⚛️ Átomos", f"{num_atomos}")
            except:
                st.metric("⚛️ Átomos", "Error")
        else:
            st.metric("⚛️ Átomos", "N/A")
    with col4:
        st.metric("🧮 Método", f"{metodo}/{conjunto_base}")
    st.markdown("---")

# Crear tabs siempre, pero mostrar contenido apropiado
tabs = st.tabs(["🔬 **Visualización 3D**", "📈 **Espectroscopía**", "⚡ **Análisis Energético**", "🔧 **Datos Técnicos**"])

with tabs[0]:
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("### 🧪 **Geometría Inicial**")
        if st.session_state.xyz_inicial:
            try:
                vista = py3Dmol.view(width=400, height=400)
                vista.addModel(st.session_state.xyz_inicial, 'xyz')
                vista.setStyle({'stick': {'radius': 0.15}, 'sphere': {'radius': 0.3}})
                vista.setBackgroundColor('#F7F7F7')
                vista.zoomTo()
                showmol(vista, height=450, width=450)
            except Exception as e:
                st.error(f"Error mostrando geometría inicial: {e}")
        else:
            st.info("💡 Carga un archivo .xyz en la barra lateral para empezar.")

    with col2:
        st.markdown("### 🎯 **Geometría Optimizada**")
        if st.session_state.xyz_optimizada:
            if not st.session_state.opt_convergida:
                st.warning("⚠️ Geometría no completamente optimizada.")
            try:
                vista_opt = py3Dmol.view(width=400, height=400)
                vista_opt.addModel(st.session_state.xyz_optimizada, 'xyz')
                vista_opt.setStyle({'stick': {'radius': 0.15}, 'sphere': {'radius': 0.3}})
                vista_opt.setBackgroundColor('#F7F7F7')
                vista_opt.zoomTo()
                showmol(vista_opt, height=450, width=450)
            except Exception as e:
                st.error(f"Error mostrando geometría optimizada: {e}")
        elif st.session_state.calculo_completado:
            st.warning("⚠️ No se pudo obtener la geometría optimizada")
        else:
            st.info("💡 Ejecuta un cálculo para ver la geometría optimizada.")

with tabs[1]:
    if st.session_state.datos_ir is not None and not st.session_state.datos_ir.empty:
        st.markdown("### 📊 **Espectro Infrarrojo (IR)**")
        try:
            fig, ax = plt.subplots(figsize=(12, 6))
            ax.stem(st.session_state.datos_ir["Frequency"], st.session_state.datos_ir["Intensity"],
                    basefmt=' ', linefmt='red', markerfmt='ro')
            ax.set_xlabel("Número de onda (cm⁻¹)")
            ax.set_ylabel("Intensidad IR (km/mol)")
            ax.set_title("Espectro IR Teórico")
            ax.invert_xaxis()
            ax.grid(True, alpha=0.3)
            st.pyplot(fig)
            st.dataframe(st.session_state.datos_ir.style.format({"Frequency": "{:.2f}", "Intensity": "{:.2f}"}),
                         use_container_width=True)
        except Exception as e:
            st.error(f"Error mostrando espectro IR: {e}")
    elif st.session_state.ultimo_tipo_calculo == "Frecuencias Vibracionales (IR)" and st.session_state.calculo_completado:
        st.warning(
            "⚠️ No se encontraron datos IR. Verifica que la optimización haya convergido en la pestaña de Datos Técnicos.")
    else:
        st.info("💡 Selecciona 'Frecuencias Vibracionales (IR)' y ejecuta un cálculo para ver el espectro.")

with tabs[2]:
    if st.session_state.calculo_completado:
        st.markdown("### ⚡ **Componentes Energéticos**")
        if st.session_state.datos_energia is not None:
            st.dataframe(st.session_state.datos_energia, use_container_width=True)
        else:
            st.info("No hay datos energéticos adicionales disponibles")

        st.markdown("### 🔋 **Energías Orbitales**")
        if st.session_state.datos_orbitales is not None:
            st.dataframe(st.session_state.datos_orbitales, use_container_width=True)
        else:
            st.info("No hay datos orbitales disponibles")

        st.markdown("### ⚛️ **Análisis de Cargas**")
        if st.session_state.datos_cargas:
            st.markdown("#### Cargas Atómicas")
            if 'Mulliken' in st.session_state.datos_cargas and not st.session_state.datos_cargas['Mulliken'].empty:
                st.write("**Cargas de Mulliken**")
                st.dataframe(st.session_state.datos_cargas['Mulliken'])
            if 'Loewdin' in st.session_state.datos_cargas and not st.session_state.datos_cargas['Loewdin'].empty:
                st.write("**Cargas de Loewdin**")
                st.dataframe(st.session_state.datos_cargas['Loewdin'])
        else:
            st.info("No hay análisis de cargas disponible")
    else:
        st.info("💡 Ejecuta un cálculo para ver el análisis detallado.")

with tabs[3]:
    if st.session_state.log_completo_orca:
        st.markdown("### 📋 **Log de Salida de ORCA**")
        st.download_button(
            label="💾 Descargar Archivo .out Completo",
            data=st.session_state.log_completo_orca,
            file_name=f"{st.session_state.nombre_trabajo}.out",
            mime="text/plain"
        )
        with st.expander("📄 Ver Log Completo"):
            st.code(st.session_state.log_completo_orca, language='text')
    else:
        st.info("💡 Ejecuta un cálculo para ver los datos técnicos.")

# Debug info (remover en producción)
if st.checkbox("🔧 Mostrar estado debug", value=False):
    st.markdown("### Estado de la sesión")
    debug_info = {
        "calculo_completado": st.session_state.calculo_completado,
        "opt_convergida": st.session_state.opt_convergida,
        "energia_final": st.session_state.energia_final,
        "tiene_xyz_optimizada": st.session_state.xyz_optimizada is not None,
        "tiene_log": st.session_state.log_completo_orca is not None,
        "datos_cargas": st.session_state.datos_cargas is not None
    }
    st.json(debug_info)

st.markdown("---")
st.markdown("*Desarrollado con Streamlit • Cálculos cuánticos con ORCA*")