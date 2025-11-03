# app.py
import streamlit as st
import os
import subprocess
import pandas as pd
import py3Dmol
from stmol import showmol
import matplotlib.pyplot as plt

from documento import generar_reporte_completo
from utils import Orca, PySCFCalculator

st.set_page_config(
    page_title="ORCA Molecular",
    layout="wide",
    initial_sidebar_state="expanded"
)

col1, col2 = st.columns([3, 1])
with col1:
    st.title("🧬 ORCA Molecular Calculator")
    st.markdown("*Calculadora cuántica para análisis molecular*")

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
if "datos_nmr" not in st.session_state:
    st.session_state.datos_nmr = None
if "datos_susceptibilidad" not in st.session_state:
    st.session_state.datos_susceptibilidad = None

if "datos_susceptibilidad" not in st.session_state:
    st.session_state.datos_susceptibilidad = None
if "pdf_generado" not in st.session_state:
    st.session_state.pdf_generado = None

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

    st.markdown("##### 🔬 Propiedades Adicionales")
    calc_nmr = st.checkbox(
        "Calcular Apantallamiento (NMR)",
        help="Calcula las propiedades de RMN (Apantallamiento Isotrópico)",
        value=False
    )

    calc_susceptibilidad = st.checkbox(
        "🧲 Calcular Susceptibilidad Magnética (PySCF)",
        help="Calcula magnetismo/diamagnetismo usando PySCF",
        value=False
    )

    st.markdown("---")

    st.markdown("#### ⚙️ **Configuración Computacional**")
    with st.expander("Parámetros Avanzados", expanded=False):
        col_a, col_b = st.columns(2)
        with col_a:
            metodo = st.selectbox("Método", ["B3LYP", "PBE0", "M06-2X", "wB97X-D"])
        with col_b:
            conjunto_base = st.selectbox("Base", ["def2-SVP", "6-31+G(d,p)", "6-311++G(d,p)", "cc-pVDZ"])

        palabras_clave = st.text_input("Palabras clave extra", "D3BJ TIGHTSCF")

    st.markdown("---")

    st.markdown("#### 🚀 **Ejecutar Cálculo**")
    boton_ejecutar = st.button(
        "🎯 **CALCULAR**",
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
        claves_a_preservar = ['xyz_inicial', 'nombre_trabajo']
        for key in st.session_state.keys():
            if key not in claves_a_preservar:
                st.session_state[key] = None if not isinstance(st.session_state[key], bool) else False

        st.session_state.ultimo_tipo_calculo = tipo_calculo
        nombre_trabajo = st.session_state.nombre_trabajo

        contenido_entrada = Orca.generar_entrada(
            st.session_state.xyz_inicial, tipo_calculo, metodo, conjunto_base, palabras_clave,
            calc_nmr=calc_nmr
        )

        ruta_entrada = os.path.join(DIR_CALCULOS, f"{nombre_trabajo}.inp")
        ruta_salida = os.path.join(DIR_CALCULOS, f"{nombre_trabajo}.out")

        with open(ruta_entrada, "w") as f:
            f.write(contenido_entrada)

        with st.spinner(f"Ejecutando ORCA para '{nombre_trabajo}'... Esto puede tardar varios minutos."):
            try:
                comando = f"orca {ruta_entrada}"
                proceso = subprocess.run(
                    comando,
                    shell=True,
                    check=True,
                    capture_output=True,
                    text=True,
                    timeout=600
                )
                with open(ruta_salida, "w") as f_out:
                    f_out.write(proceso.stdout)
                    f_out.write(proceso.stderr)

                st.session_state.calculo_completado = True
            except subprocess.TimeoutExpired:
                st.error("El cálculo de ORCA tardó demasiado (más de 10 minutos) y fue cancelado.")
            except subprocess.CalledProcessError as e:
                st.error("Error al ejecutar ORCA. Revisa los parámetros y el log.")
                st.code(e.stderr)
            except Exception as e:
                st.error(f"Error inesperado: {e}")

        if os.path.exists(ruta_salida):
            try:
                analizador = Orca(ruta_salida)

                st.session_state.log_completo_orca = analizador.contenido
                st.session_state.resumen_log_orca = "".join(analizador.contenido.splitlines(True)[-50:])

                st.session_state.opt_convergida = analizador.verificar_convergencia()
                st.session_state.xyz_optimizada = analizador.extraer_geometria_optimizada()
                st.session_state.energia_final = analizador.extraer_energia_final()
                st.session_state.datos_energia = analizador.extraer_componentes_energia()
                st.session_state.datos_cargas = analizador.extraer_cargas_atomicas()
                st.session_state.datos_orbitales = analizador.extraer_energias_orbitales()
                st.session_state.datos_cargas_reducidas = analizador.extraer_cargas_orbitales_reducidas()
                st.session_state.datos_nmr = analizador.extraer_datos_nmr()

                if tipo_calculo == "Frecuencias Vibracionales (IR)":
                    st.session_state.datos_ir = analizador.extraer_espectro_ir(factor_escalamiento)

            except Exception as e:
                st.error(f"Ocurrió un error al analizar el archivo de salida: {e}")

        if calc_susceptibilidad:
            xyz_para_pyscf = st.session_state.xyz_optimizada if st.session_state.xyz_optimizada else st.session_state.xyz_inicial

            with st.spinner("🧲 Calculando susceptibilidad magnética con PySCF..."):
                resultados = PySCFCalculator.calcular_susceptibilidad(
                    xyz_para_pyscf,
                    metodo=metodo,
                    base=conjunto_base
                )
                st.session_state.datos_susceptibilidad = resultados

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

tabs = st.tabs(["🔬 **Visualización 3D**", "📈 **Espectroscopía**", "🧲 **Magnetismo**", "⚡ **Análisis Energético**",
                "🔧 **Datos Técnicos**"])

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
    ir_disponible = st.session_state.datos_ir is not None and not st.session_state.datos_ir.empty
    nmr_disponible = st.session_state.datos_nmr is not None and not st.session_state.datos_nmr.empty

    if ir_disponible:
        st.markdown("### 📊 **Espectro Infrarrojo (IR)**")
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.stem(st.session_state.datos_ir["Frequency"], st.session_state.datos_ir["Intensity"], basefmt=' ',
                linefmt='red', markerfmt='ro')
        ax.set_xlabel("Número de onda (cm⁻¹)")
        ax.set_ylabel("Intensidad IR (km/mol)")
        ax.set_title("Espectro IR Teórico")
        ax.invert_xaxis()
        ax.grid(True, alpha=0.3)
        st.pyplot(fig)
        st.dataframe(st.session_state.datos_ir.style.format({"Frequency": "{:.2f}", "Intensity": "{:.2f}"}))

    elif not ir_disponible and st.session_state.ultimo_tipo_calculo == "Frecuencias Vibracionales (IR)":
        st.warning("⚠️ No se encontraron datos IR. Verifica que la optimización haya convergido.")

    if ir_disponible and nmr_disponible:
        st.markdown("---")

    if nmr_disponible:
        st.markdown("### 🛡️ **Apantallamiento Nuclear (NMR)**")
        st.info("Valores de apantallamiento isotrópico (ppm). Valores más altos indican mayor apantallamiento.")
        st.dataframe(st.session_state.datos_nmr.style.format({
            "Núcleo": "{}",
            "Elemento": "{}",
            "Isotrópico (ppm)": "{:.3f}",
            "Anisotropía (ppm)": "{:.3f}"
        }), use_container_width=True)

    if not ir_disponible and not nmr_disponible:
        st.info(
            "💡 Selecciona 'Frecuencias Vibracionales (IR)' y/o 'Calcular Apantallamiento (NMR)' en la barra lateral.")

with tabs[2]:
    if st.session_state.datos_susceptibilidad is None:
        st.info("💡 Activa '🧲 Calcular Susceptibilidad Magnética (PySCF)' en la barra lateral y ejecuta un cálculo.")
    elif "error" in st.session_state.datos_susceptibilidad:
        st.error(f"❌ {st.session_state.datos_susceptibilidad['error']}")
    else:
        datos = st.session_state.datos_susceptibilidad

        st.markdown("### 🧲 **Susceptibilidad Magnética Molecular**")

        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric(
                "χ Isotrópica (CGS)",
                f"{datos['isotropico_cgs']:.2f}",
                help="Susceptibilidad magnética en unidades de 10⁻⁶ cm³/mol"
            )
        with col2:
            st.metric(
                "Tipo de Magnetismo",
                datos['tipo'],
                help="Diamagnético (χ < 0) o Paramagnético (χ > 0)"
            )
        with col3:
            st.metric(
                "χ (a.u.)",
                f"{datos['isotropico_au']:.6f}",
                help="Susceptibilidad en unidades atómicas"
            )

        st.markdown("---")

        st.markdown("#### 📊 **Tensor de Susceptibilidad Magnética (a.u.)**")
        tensor_df = pd.DataFrame(
            datos['tensor'],
            columns=['X', 'Y', 'Z'],
            index=['X', 'Y', 'Z']
        )
        st.dataframe(tensor_df.style.format("{:.6f}"), use_container_width=True)

        st.markdown("#### 📈 **Componentes del Tensor**")
        fig, ax = plt.subplots(figsize=(10, 6))
        componentes = ['χ_XX', 'χ_YY', 'χ_ZZ']
        valores = [datos['tensor'][0][0], datos['tensor'][1][1], datos['tensor'][2][2]]
        colores = ['#FF6B6B' if v < 0 else '#4ECDC4' for v in valores]

        ax.bar(componentes, valores, color=colores, alpha=0.7, edgecolor='black')
        ax.axhline(y=0, color='black', linestyle='--', linewidth=0.8)
        ax.set_ylabel('Susceptibilidad (a.u.)')
        ax.set_title('Componentes Diagonales del Tensor χ')
        ax.grid(True, alpha=0.3, axis='y')
        st.pyplot(fig)

        st.markdown("---")
        st.markdown("#### ℹ️ **Interpretación**")
        if datos['tipo'] == "Diamagnético":
            st.info("""
            **Sustancia Diamagnética** (χ < 0):
            - Repelida débilmente por campos magnéticos
            - Todos los electrones están apareados
            - Ejemplo: H₂O, CH₄, NaCl
            """)
        else:
            st.success("""
            **Sustancia Paramagnética** (χ > 0):
            - Atraída por campos magnéticos
            - Presencia de electrones desapareados
            - Ejemplo: O₂, NO, radicales libres
            """)
        if 'nota' in datos:
            st.info(f"ℹ️ **Nota:** {datos['nota']}")

with tabs[3]:
    if not st.session_state.calculo_completado:
        st.info("💡 Ejecuta un cálculo para ver el análisis detallado.")
    else:
        st.markdown("### ⚡ **Componentes Energéticos**")
        if st.session_state.datos_energia is not None and not st.session_state.datos_energia.empty:
            st.dataframe(st.session_state.datos_energia)

        st.markdown("### 🔋 **Energías Orbitales**")
        if st.session_state.datos_orbitales is not None and not st.session_state.datos_orbitales.empty:
            st.dataframe(st.session_state.datos_orbitales)

        st.markdown("### ⚛️ **Análisis de Cargas**")
        col1, col2 = st.columns(2)
        with col1:
            st.markdown("#### Cargas Atómicas")
            if st.session_state.datos_cargas:
                if 'Mulliken' in st.session_state.datos_cargas:
                    st.write("**Cargas de Mulliken**")
                    st.dataframe(st.session_state.datos_cargas['Mulliken'])
                if 'Loewdin' in st.session_state.datos_cargas:
                    st.write("**Cargas de Loewdin**")
                    st.dataframe(st.session_state.datos_cargas['Loewdin'])
        with col2:
            st.markdown("#### Cargas Orbitales Reducidas")
            if st.session_state.datos_cargas_reducidas:
                if 'Mulliken' in st.session_state.datos_cargas_reducidas:
                    st.write("**Mulliken (Reducidas)**")
                    st.dataframe(st.session_state.datos_cargas_reducidas['Mulliken'])
                if 'Loewdin' in st.session_state.datos_cargas_reducidas:
                    st.write("**Loewdin (Reducidas)**")
                    st.dataframe(st.session_state.datos_cargas_reducidas['Loewdin'])

with tabs[4]:
    if not st.session_state.calculo_completado:
        st.info("💡 Ejecuta un cálculo para generar el reporte.")
    else:
        st.markdown("### 📄 **Generación de Reporte PDF**")
        st.markdown("---")

        col1, col2 = st.columns([2, 1])

        with col1:
            st.markdown("#### 📋 **Contenido del Reporte**")

            secciones_incluidas = []

            if st.session_state.energia_final is not None:
                secciones_incluidas.append("✅ Resultados Energéticos")

            if st.session_state.datos_energia is not None:
                secciones_incluidas.append("✅ Componentes de Energía")

            if st.session_state.datos_ir is not None and not st.session_state.datos_ir.empty:
                secciones_incluidas.append("✅ Espectro Infrarrojo (IR)")

            if st.session_state.datos_nmr is not None and not st.session_state.datos_nmr.empty:
                secciones_incluidas.append("✅ Apantallamiento Nuclear (NMR)")

            if st.session_state.datos_susceptibilidad is not None:
                if 'error' not in st.session_state.datos_susceptibilidad:
                    secciones_incluidas.append("✅ Susceptibilidad Magnética")

            if st.session_state.datos_cargas is not None:
                secciones_incluidas.append("✅ Análisis de Cargas Atómicas")

            if st.session_state.datos_orbitales is not None and not st.session_state.datos_orbitales.empty:
                secciones_incluidas.append("✅ Energías Orbitales (HOMO-LUMO)")

            if secciones_incluidas:
                st.markdown("**El reporte incluirá:**")
                for seccion in secciones_incluidas:
                    st.markdown(f"- {seccion}")
            else:
                st.warning("⚠️ No hay datos suficientes para generar el reporte.")

        with col2:
            st.markdown("#### ⚙️ **Configuración**")

            st.info(f"""
            **Molécula:** {st.session_state.nombre_trabajo}

            **Método:** {metodo}

            **Base:** {conjunto_base}

            **Estado:** {'✅ Convergido' if st.session_state.opt_convergida else '⚠️ No convergido'}
            """)

        st.markdown("---")

        st.markdown("#### 💾 **Descargar Reporte**")

        col_btn1, col_btn2, col_btn3 = st.columns([1, 1, 2])

        with col_btn1:
            if st.button("📥 **Generar PDF**", type="primary", use_container_width=True):
                with st.spinner("🔄 Generando reporte PDF..."):
                    try:
                        pdf_buffer = generar_reporte_completo(
                            nombre_trabajo=st.session_state.nombre_trabajo,
                            metodo=metodo,
                            base=conjunto_base,
                            energia_final=st.session_state.energia_final,
                            convergida=st.session_state.opt_convergida,
                            datos_energia=st.session_state.datos_energia,
                            datos_ir=st.session_state.datos_ir,
                            factor_escalamiento=factor_escalamiento if tipo_calculo == "Frecuencias Vibracionales (IR)" else 1.0,
                            datos_nmr=st.session_state.datos_nmr,
                            datos_susceptibilidad=st.session_state.datos_susceptibilidad,
                            datos_cargas=st.session_state.datos_cargas,
                            datos_orbitales=st.session_state.datos_orbitales
                        )

                        st.session_state.pdf_generado = pdf_buffer
                        st.success("✅ Reporte PDF generado exitosamente!")

                    except Exception as e:
                        st.error(f"❌ Error al generar PDF: {str(e)}")
                        import traceback

                        st.code(traceback.format_exc())

        with col_btn2:
            if 'pdf_generado' in st.session_state and st.session_state.pdf_generado is not None:
                st.download_button(
                    label="📄 Descargar PDF",
                    data=st.session_state.pdf_generado,
                    file_name=f"{st.session_state.nombre_trabajo}_reporte.pdf",
                    mime="application/pdf",
                    use_container_width=True
                )

st.markdown("---")
st.markdown("*Desarrollado con Streamlit • Cálculos cuánticos con ORCA y PySCF*")