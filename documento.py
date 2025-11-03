from reportlab.lib import colors
from reportlab.lib.pagesizes import letter, A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                TableStyle, PageBreak, Image, KeepTogether)
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_JUSTIFY
from reportlab.pdfgen import canvas
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('Agg')  # Backend sin GUI
import io
from datetime import datetime
import numpy as np


class GeneradorReportePDF:
    """Genera reportes PDF profesionales con resultados de calculos moleculares"""

    def __init__(self, nombre_trabajo, metodo, base):
        self.nombre_trabajo = nombre_trabajo
        self.metodo = metodo
        self.base = base
        self.elementos = []
        self.estilos = getSampleStyleSheet()
        self._configurar_estilos()

    def _configurar_estilos(self):
        """Configura estilos personalizados para el documento"""
        # Estilo para titulo principal
        self.estilos.add(ParagraphStyle(
            name='TituloPrincipal',
            parent=self.estilos['Heading1'],
            fontSize=24,
            textColor=colors.HexColor('#1f77b4'),
            spaceAfter=30,
            alignment=TA_CENTER,
            fontName='Helvetica-Bold'
        ))

        # Estilo para subtitulos
        self.estilos.add(ParagraphStyle(
            name='Subtitulo',
            parent=self.estilos['Heading2'],
            fontSize=16,
            textColor=colors.HexColor('#2c3e50'),
            spaceAfter=12,
            spaceBefore=12,
            fontName='Helvetica-Bold'
        ))

        # Estilo para texto normal
        self.estilos.add(ParagraphStyle(
            name='TextoNormal',
            parent=self.estilos['Normal'],
            fontSize=10,
            alignment=TA_JUSTIFY,
            spaceAfter=6
        ))

        # Estilo para metadatos
        self.estilos.add(ParagraphStyle(
            name='Metadata',
            parent=self.estilos['Normal'],
            fontSize=9,
            textColor=colors.grey,
            alignment=TA_CENTER
        ))

    def agregar_portada(self):
        """Agrega portada al documento"""
        # Titulo
        titulo = Paragraph(
            f"<b>Reporte de Calculo Molecular</b>",
            self.estilos['TituloPrincipal']
        )
        self.elementos.append(titulo)
        self.elementos.append(Spacer(1, 0.3 * inch))

        # Nombre del trabajo
        nombre = Paragraph(
            f"<b>Molecula:</b> {self.nombre_trabajo}",
            self.estilos['Heading2']
        )
        self.elementos.append(nombre)
        self.elementos.append(Spacer(1, 0.2 * inch))

        # Metadatos
        fecha = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        metadata = [
            f"<b>Metodo:</b> {self.metodo}",
            f"<b>Conjunto de Base:</b> {self.base}",
            f"<b>Fecha:</b> {fecha}",
            f"<b>Software:</b> ORCA + PySCF"
        ]

        for item in metadata:
            self.elementos.append(Paragraph(item, self.estilos['TextoNormal']))
            self.elementos.append(Spacer(1, 0.1 * inch))

        self.elementos.append(Spacer(1, 0.5 * inch))

        # Linea divisoria
        self.elementos.append(Spacer(1, 0.3 * inch))
        self.elementos.append(PageBreak())

    def agregar_seccion_energia(self, energia_final, convergida, datos_energia):
        """Agrega seccion de resultados energeticos"""
        self.elementos.append(Paragraph("1. Resultados Energeticos", self.estilos['Subtitulo']))

        # Energia final
        estado = "Convergido" if convergida else "No Convergido"
        color_estado = "green" if convergida else "red"

        datos = [
            ["<b>Parametro</b>", "<b>Valor</b>"],
            ["Energia Final", f"{energia_final:.6f} Hartree"],
            ["Estado de Convergencia", f'<font color="{color_estado}">{estado}</font>'],
            ["Energia (eV)", f"{energia_final * 27.211:.4f} eV"],
            ["Energia (kcal/mol)", f"{energia_final * 627.5:.2f} kcal/mol"]
        ]

        tabla = Table(datos, colWidths=[3 * inch, 3 * inch])
        tabla.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#3498db')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 11),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ]))

        self.elementos.append(tabla)
        self.elementos.append(Spacer(1, 0.2 * inch))

        # Componentes de energia
        if datos_energia is not None and not datos_energia.empty:
            self.elementos.append(Paragraph("Componentes Energeticos", self.estilos['Heading3']))

            tabla_datos = [["<b>Componente</b>", "<b>Energia (Hartree)</b>"]]
            for idx, row in datos_energia.iterrows():
                tabla_datos.append([idx, f"{row['Energia (Hartree)']:.6f}"])

            tabla_comp = Table(tabla_datos, colWidths=[3 * inch, 2.5 * inch])
            tabla_comp.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#2ecc71')),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('GRID', (0, 0), (-1, -1), 1, colors.black),
                ('BACKGROUND', (0, 1), (-1, -1), colors.lightgrey),
            ]))

            self.elementos.append(tabla_comp)

        self.elementos.append(Spacer(1, 0.3 * inch))

    def agregar_espectro_ir(self, datos_ir, factor_escalamiento):
        """Agrega grafico de espectro IR"""
        if datos_ir is None or datos_ir.empty:
            return

        self.elementos.append(PageBreak())
        self.elementos.append(Paragraph("2. Espectro Infrarrojo (IR)", self.estilos['Subtitulo']))

        # Generar grafico
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.stem(datos_ir["Frequency"], datos_ir["Intensity"],
                basefmt=' ', linefmt='red', markerfmt='ro')
        ax.set_xlabel("Numero de onda (cm⁻¹)", fontsize=10)
        ax.set_ylabel("Intensidad IR (km/mol)", fontsize=10)
        ax.set_title(f"Espectro IR Teorico (Factor: {factor_escalamiento})", fontsize=12)
        ax.invert_xaxis()
        ax.grid(True, alpha=0.3)

        # Guardar en buffer
        buf = io.BytesIO()
        plt.tight_layout()
        plt.savefig(buf, format='png', dpi=150)
        buf.seek(0)
        plt.close()

        # Agregar imagen
        img = Image(buf, width=6.5 * inch, height=3.25 * inch)
        self.elementos.append(img)
        self.elementos.append(Spacer(1, 0.2 * inch))

        # Tabla de frecuencias principales
        self.elementos.append(Paragraph("Frecuencias Principales", self.estilos['Heading3']))

        # Ordenar por intensidad y tomar top 10
        top_freqs = datos_ir.nlargest(10, 'Intensity')

        tabla_datos = [["<b>Frecuencia (cm⁻¹)</b>", "<b>Intensidad (km/mol)</b>"]]
        for _, row in top_freqs.iterrows():
            tabla_datos.append([
                f"{row['Frequency']:.2f}",
                f"{row['Intensity']:.2f}"
            ])

        tabla = Table(tabla_datos, colWidths=[2.5 * inch, 2.5 * inch])
        tabla.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#e74c3c')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('BACKGROUND', (0, 1), (-1, -1), colors.lightgrey),
        ]))

        self.elementos.append(tabla)
        self.elementos.append(Spacer(1, 0.3 * inch))

    def agregar_datos_nmr(self, datos_nmr):
        """Agrega datos de apantallamiento NMR"""
        if datos_nmr is None or datos_nmr.empty:
            return

        self.elementos.append(PageBreak())
        self.elementos.append(Paragraph("3. Apantallamiento Nuclear (NMR)", self.estilos['Subtitulo']))

        # Convertir DataFrame a lista para tabla
        tabla_datos = [["<b>Nucleo</b>", "<b>Elemento</b>",
                        "<b>Isotropico (ppm)</b>", "<b>Anisotropia (ppm)</b>"]]

        for _, row in datos_nmr.iterrows():
            tabla_datos.append([
                str(row['Nucleo']),
                row['Elemento'],
                f"{row['Isotropico (ppm)']:.3f}",
                f"{row['Anisotropia (ppm)']:.3f}"
            ])

        tabla = Table(tabla_datos, colWidths=[1.2 * inch, 1.2 * inch, 2 * inch, 2 * inch])
        tabla.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#9b59b6')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 9),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('BACKGROUND', (0, 1), (-1, -1), colors.lightgrey),
        ]))

        self.elementos.append(tabla)
        self.elementos.append(Spacer(1, 0.3 * inch))

    def agregar_susceptibilidad(self, datos_susc):
        """Agrega analisis de susceptibilidad magnetica"""
        if datos_susc is None or 'error' in datos_susc:
            return

        self.elementos.append(PageBreak())
        self.elementos.append(Paragraph("4. Susceptibilidad Magnetica", self.estilos['Subtitulo']))

        # Tabla de resultados principales
        datos = [
            ["<b>Propiedad</b>", "<b>Valor</b>"],
            ["Chi Isotropica (a.u.)", f"{datos_susc['isotropico_au']:.6f}"],
            ["Chi Isotropica (CGS)", f"{datos_susc['isotropico_cgs']:.2f} × 10⁻⁶ cm³/mol"],
            ["Tipo de Magnetismo", datos_susc['tipo']],
            ["Metodo de Calculo", datos_susc['metodo_calculo']]
        ]

        tabla = Table(datos, colWidths=[3 * inch, 3 * inch])
        tabla.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#34495e')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('BACKGROUND', (0, 1), (-1, -1), colors.lightgrey),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ]))

        self.elementos.append(tabla)
        self.elementos.append(Spacer(1, 0.2 * inch))

        # Tensor de susceptibilidad
        self.elementos.append(Paragraph("Tensor de Susceptibilidad (a.u.)", self.estilos['Heading3']))

        tensor = np.array(datos_susc['tensor'])
        tabla_datos = [["", "<b>X</b>", "<b>Y</b>", "<b>Z</b>"]]

        for i, label in enumerate(['<b>X</b>', '<b>Y</b>', '<b>Z</b>']):
            fila = [label]
            for j in range(3):
                fila.append(f"{tensor[i][j]:.6f}")
            tabla_datos.append(fila)

        tabla_tensor = Table(tabla_datos, colWidths=[1 * inch, 1.5 * inch, 1.5 * inch, 1.5 * inch])
        tabla_tensor.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#16a085')),
            ('BACKGROUND', (0, 1), (0, -1), colors.HexColor('#16a085')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('TEXTCOLOR', (0, 1), (0, -1), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTNAME', (0, 1), (0, -1), 'Helvetica-Bold'),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('BACKGROUND', (1, 1), (-1, -1), colors.lightgrey),
        ]))

        self.elementos.append(tabla_tensor)
        self.elementos.append(Spacer(1, 0.2 * inch))

        # Grafico de componentes
        fig, ax = plt.subplots(figsize=(6, 3.5))
        componentes = ['χ_XX', 'χ_YY', 'χ_ZZ']
        valores = [tensor[0][0], tensor[1][1], tensor[2][2]]
        colores_bar = ['#e74c3c' if v < 0 else '#2ecc71' for v in valores]

        ax.bar(componentes, valores, color=colores_bar, alpha=0.7, edgecolor='black')
        ax.axhline(y=0, color='black', linestyle='--', linewidth=0.8)
        ax.set_ylabel('Susceptibilidad (a.u.)', fontsize=10)
        ax.set_title('Componentes Diagonales del Tensor χ', fontsize=11)
        ax.grid(True, alpha=0.3, axis='y')

        buf = io.BytesIO()
        plt.tight_layout()
        plt.savefig(buf, format='png', dpi=150)
        buf.seek(0)
        plt.close()

        img = Image(buf, width=5 * inch, height=2.9 * inch)
        self.elementos.append(img)

        # Nota sobre aproximacion
        if 'nota' in datos_susc:
            self.elementos.append(Spacer(1, 0.2 * inch))
            nota = Paragraph(
                f"<i><b>Nota:</b> {datos_susc['nota']}</i>",
                self.estilos['TextoNormal']
            )
            self.elementos.append(nota)

        self.elementos.append(Spacer(1, 0.3 * inch))

    def agregar_cargas(self, datos_cargas):
        """Agrega analisis de cargas atomicas"""
        if datos_cargas is None:
            return

        self.elementos.append(PageBreak())
        self.elementos.append(Paragraph("5. Analisis de Cargas Atomicas", self.estilos['Subtitulo']))

        for tipo, df in datos_cargas.items():
            self.elementos.append(Paragraph(f"Cargas de {tipo}", self.estilos['Heading3']))

            # Limitar a primeros 20 atomos para no llenar mucho
            df_display = df.head(20)

            tabla_datos = [["<b>Atomo</b>", "<b>Carga</b>"]]
            for _, row in df_display.iterrows():
                tabla_datos.append([
                    row['Atomo'],
                    f"{row['Carga']:.4f}"
                ])

            tabla = Table(tabla_datos, colWidths=[2 * inch, 2 * inch])
            tabla.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#f39c12')),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('GRID', (0, 0), (-1, -1), 1, colors.black),
                ('BACKGROUND', (0, 1), (-1, -1), colors.lightgrey),
                ('FONTSIZE', (0, 0), (-1, -1), 8),
            ]))

            self.elementos.append(tabla)

            if len(df) > 20:
                nota = Paragraph(f"<i>(Mostrando primeros 20 de {len(df)} atomos)</i>",
                                 self.estilos['TextoNormal'])
                self.elementos.append(nota)

            self.elementos.append(Spacer(1, 0.2 * inch))

    def agregar_orbitales(self, datos_orbitales):
        """Agrega informacion de orbitales moleculares"""
        if datos_orbitales is None or datos_orbitales.empty:
            return

        self.elementos.append(PageBreak())
        self.elementos.append(Paragraph("6. Energias Orbitales", self.estilos['Subtitulo']))

        # Encontrar HOMO y LUMO
        ocupados = datos_orbitales[datos_orbitales['Ocupacion'] > 0]
        vacios = datos_orbitales[datos_orbitales['Ocupacion'] == 0]

        if not ocupados.empty and not vacios.empty:
            homo = ocupados.iloc[-1]
            lumo = vacios.iloc[0]
            gap = lumo['Energia (eV)'] - homo['Energia (eV)']

            # Tabla HOMO-LUMO
            datos = [
                ["<b>Orbital</b>", "<b>Numero</b>", "<b>Energia (eV)</b>"],
                ["HOMO", str(homo['Numero']), f"{homo['Energia (eV)']:.4f}"],
                ["LUMO", str(lumo['Numero']), f"{lumo['Energia (eV)']:.4f}"],
                ["Gap HOMO-LUMO", "-", f"{gap:.4f}"]
            ]

            tabla = Table(datos, colWidths=[2 * inch, 1.5 * inch, 2 * inch])
            tabla.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#8e44ad')),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('GRID', (0, 0), (-1, -1), 1, colors.black),
                ('BACKGROUND', (0, 1), (-1, -1), colors.lightgrey),
            ]))

            self.elementos.append(tabla)
            self.elementos.append(Spacer(1, 0.3 * inch))

    def generar_pdf(self, nombre_archivo):
        """Genera el archivo PDF final"""
        doc = SimpleDocTemplate(
            nombre_archivo,
            pagesize=letter,
            rightMargin=0.75 * inch,
            leftMargin=0.75 * inch,
            topMargin=0.75 * inch,
            bottomMargin=0.75 * inch
        )

        doc.build(self.elementos)
        return nombre_archivo


def generar_reporte_completo(
        nombre_trabajo,
        metodo,
        base,
        energia_final=None,
        convergida=False,
        datos_energia=None,
        datos_ir=None,
        factor_escalamiento=1.0,
        datos_nmr=None,
        datos_susceptibilidad=None,
        datos_cargas=None,
        datos_orbitales=None
):
    """
    Funcion principal para generar reporte PDF completo

    Returns:
        BytesIO: Buffer con el PDF generado
    """

    generador = GeneradorReportePDF(nombre_trabajo, metodo, base)

    # Agregar portada
    generador.agregar_portada()

    # Agregar secciones segun datos disponibles
    if energia_final is not None:
        generador.agregar_seccion_energia(energia_final, convergida, datos_energia)

    if datos_ir is not None and not datos_ir.empty:
        generador.agregar_espectro_ir(datos_ir, factor_escalamiento)

    if datos_nmr is not None and not datos_nmr.empty:
        generador.agregar_datos_nmr(datos_nmr)

    if datos_susceptibilidad is not None:
        generador.agregar_susceptibilidad(datos_susceptibilidad)

    if datos_cargas is not None:
        generador.agregar_cargas(datos_cargas)

    if datos_orbitales is not None and not datos_orbitales.empty:
        generador.agregar_orbitales(datos_orbitales)

    # Generar PDF en buffer
    buffer = io.BytesIO()
    temp_filename = f"{nombre_trabajo}_reporte.pdf"
    generador.generar_pdf(temp_filename)

    # Leer archivo y guardar en buffer
    with open(temp_filename, 'rb') as f:
        buffer.write(f.read())

    buffer.seek(0)

    # Opcional: limpiar archivo temporal
    import os
    if os.path.exists(temp_filename):
        os.remove(temp_filename)

    return buffer