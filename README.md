# 🧬 ORCA Molecular Calculator

Calculadora cuántica avanzada para análisis molecular con interfaz Streamlit. Combina ORCA y PySCF para cálculos de química computacional.

---

## 📋 Tabla de Contenidos

- [Características](#características)
- [Requisitos](#requisitos)
- [Arquitectura del Software](#arquitectura-del-software)
- [Fundamentos Matemáticos](#fundamentos-matemáticos)
- [Formato de Archivos](#formato-de-archivos)
- [Ejemplos](#ejemplos)


---

## ✨ Características

### Cálculos con ORCA
- ✅ **Optimización de Geometría**: Encuentra la estructura molecular de menor energía
- ✅ **Frecuencias Vibracionales (IR)**: Calcula espectro infrarrojo teórico
- ✅ **Apantallamiento Nuclear (NMR)**: Calcula desplazamientos químicos de RMN
- ✅ **Análisis Energético**: Energías orbitales, repulsión nuclear, energía electrónica
- ✅ **Cargas Atómicas**: Análisis de Mulliken y Löwdin

### Cálculos con PySCF
- ✅ **Susceptibilidad Magnética**: Determina si la molécula es diamagnética o paramagnética (aproximación de Pascal)

### Visualización
- 🎨 **Visualización 3D**: Geometría inicial y optimizada con py3Dmol
- 📊 **Gráficos IR**: Espectro infrarrojo interactivo
- 🧲 **Análisis de Magnetismo**: Tensor de susceptibilidad y componentes

---

## 🔧 Requisitos

### Software Requerido
- **Python 3.8+**
- **ORCA 5.0+** (instalado y accesible desde línea de comandos)
- **Anaconda/Miniconda** (recomendado)

### Librerías Python
```bash
streamlit
pandas
numpy
matplotlib
py3Dmol
stmol
pyscf
```

---

## 🏗️ Arquitectura del Software

### División de Tareas: ORCA vs PySCF

| Tarea | Software | Razón |
|-------|----------|-------|
| **Optimización de Geometría** | ORCA | Motor robusto para convergencia |
| **Frecuencias Vibracionales** | ORCA | Cálculo preciso de Hessianas |
| **Apantallamiento NMR** | ORCA | Implementación GIAO completa |
| **Análisis de Cargas** | ORCA | Métodos de población integrados |
| **Energías Orbitales** | ORCA | Parser de salida completo |
| **Susceptibilidad Magnética** | PySCF | Aproximación cuando ORCA no disponible |

---

## 📐 Fundamentos Matemáticos

### 1. Optimización de Geometría (ORCA)

ORCA minimiza la energía total usando métodos DFT:

```
E_total = E_nuclear-repulsion + E_electronic

E_electronic = E_kinetic + E_nuclear-attraction + E_electron-repulsion + E_XC
```

Donde:
- **E_XC**: Energía de intercambio-correlación (funcional B3LYP, PBE0, etc.)
- Convergencia cuando: `|∇E| < umbral` y `|ΔE| < umbral`

### 2. Frecuencias Vibracionales (ORCA)

Calcula la matriz Hessiana (segundas derivadas de energía):

```
H_ij = ∂²E / ∂x_i ∂x_j
```

Las frecuencias vienen de los eigenvalores:

```
ω_i = √(λ_i / μ_i)
```

**Factor de escalamiento**: Corrige errores sistemáticos del método DFT
```
ω_experimental ≈ factor × ω_calculado
```

Ejemplo: B3LYP/def2-SVP → factor ≈ 0.9679

### 3. Apantallamiento Nuclear - NMR (ORCA)

Usa Gauge-Including Atomic Orbitals (GIAO):

```
σ_núcleo = σ_diamagnético + σ_paramagnético

σ_diamagnético ∝ ⟨ψ| r²/r³ |ψ⟩
σ_paramagnético ∝ Σ_occ,virt ⟨occ|L|virt⟩⟨virt|L|occ⟩ / ΔE
```

Donde:
- **σ**: Constante de apantallamiento (ppm)
- **L**: Operador de momento angular
- **ΔE**: Diferencia de energía entre orbitales

### 4. Susceptibilidad Magnética (PySCF - Aproximación)

**Método**: Aproximación de Pascal (contribución diamagnética)

```
χ_tensor = - (e² / 6mc²) Σ_i Z_i [(r_i²)δ_jk - r_ij r_ik]
```

Donde:
- **Z_i**: Carga nuclear del átomo i
- **r_i**: Vector posición desde el centro de masa
- **δ_jk**: Delta de Kronecker

**Valor isotrópico**:
```
χ_iso = (χ_XX + χ_YY + χ_ZZ) / 3
```

**Conversión a CGS**:
```
χ_CGS (10⁻⁶ cm³/mol) = χ_a.u. × 0.78910
```

**Clasificación**:
- **χ < 0**: Diamagnético (repelido por campo magnético)
- **χ > 0**: Paramagnético (atraído por campo magnético)

⚠️ **Nota**: Esta es una aproximación clásica. Para resultados precisos usar ORCA con palabras clave NMR completas.

### 5. Cargas Atómicas

**Análisis de Mulliken**:
```
q_A = Z_A - Σ_μ∈A [P_μμ + Σ_ν∈B≠A P_μν S_μν]
```

**Análisis de Löwdin**:
```
q_A = Z_A - Σ_μ∈A (PS^(1/2))_μμ
```

Donde:
- **P**: Matriz de densidad
- **S**: Matriz de superposición
- **Z_A**: Carga nuclear

---

## 📄 Formato de Archivos

### Archivo XYZ - Especificaciones

El formato XYZ sigue el estándar de química computacional:

```
<número_de_átomos>
<línea_de_comentario>
<símbolo_1> <x_1> <y_1> <z_1>
<símbolo_2> <x_2> <y_2> <z_2>
...
```

#### ⚠️ **Reglas Importantes**

1. **Línea 1**: DEBE ser un número entero (cantidad de átomos)
2. **Línea 2**: Comentario (puede estar vacía, pero DEBE existir)
3. **Líneas 3+**: Coordenadas atómicas (símbolo y 3 coordenadas)
4. **Separador**: Espacios o tabulaciones
5. **Coordenadas**: En Ångströms (Å)
6. **Símbolos**: Usar nomenclatura estándar (H, C, N, O, etc.)

#### ✅ **Formato Correcto**

```xyz
3
Molecula de agua
O    0.000000    0.000000    0.119262
H    0.000000    0.763239   -0.477047
H    0.000000   -0.763239   -0.477047
```

#### ❌ **Formatos Incorrectos**

**Error 1: Sin línea de comentario**
```xyz
3
O    0.000000    0.000000    0.119262
H    0.000000    0.763239   -0.477047
```
❌ Falta la línea 2 de comentario

**Error 2: Número de átomos incorrecto**
```xyz
2
Agua
O    0.000000    0.000000    0.119262
H    0.000000    0.763239   -0.477047
H    0.000000   -0.763239   -0.477047
```
❌ Dice 2 átomos pero hay 3

**Error 3: Coordenadas faltantes**
```xyz
3
Agua
O    0.000000    0.000000
H    0.000000    0.763239   -0.477047
H    0.000000   -0.763239   -0.477047
```
❌ Primer átomo solo tiene 2 coordenadas

---

## 📚 Ejemplos

### Ejemplo 1: Molécula de Agua (H₂O)

**Archivo: `water.xyz`**
```xyz
3
Water molecule - equilibrium geometry
O     0.000000     0.000000     0.119262
H     0.000000     0.763239    -0.477047
H     0.000000    -0.763239    -0.477047
```

**Configuración recomendada**:
- Método: B3LYP
- Base: def2-SVP
- Cálculo: Frecuencias Vibracionales (IR)
- ☑️ Calcular NMR
- ☑️ Calcular Susceptibilidad

**Resultados esperados**:
- Energía: ~ -76.4 Hartree
- Tipo: Diamagnética (χ < 0)
- 3 frecuencias IR activas

---

### Ejemplo 2: Metano (CH₄)

**Archivo: `methane.xyz`**
```xyz
5
Methane - tetrahedral structure
C     0.000000     0.000000     0.000000
H     0.631391     0.631391     0.631391
H    -0.631391    -0.631391     0.631391
H     0.631391    -0.631391    -0.631391
H    -0.631391     0.631391    -0.631391
```

**Configuración recomendada**:
- Método: PBE0
- Base: 6-31+G(d,p)
- Cálculo: Optimización de Geometría

---

### Ejemplo 3: Benceno (C₆H₆)

**Archivo: `benzene.xyz`**
```xyz
12
Benzene ring
C     0.000000     1.396000     0.000000
C     1.209000     0.698000     0.000000
C     1.209000    -0.698000     0.000000
C     0.000000    -1.396000     0.000000
C    -1.209000    -0.698000     0.000000
C    -1.209000     0.698000     0.000000
H     0.000000     2.479000     0.000000
H     2.147000     1.240000     0.000000
H     2.147000    -1.240000     0.000000
H     0.000000    -2.479000     0.000000
H    -2.147000    -1.240000     0.000000
H    -2.147000     1.240000     0.000000
```

**Configuración recomendada**:
- Método: B3LYP
- Base: def2-TZVP
- Cálculo: Frecuencias Vibracionales
- Factor IR: 0.9679

---

### Ejemplo 4: Radical OH (Paramagnético)

**Archivo: `oh_radical.xyz`**
```xyz
2
Hydroxyl radical - open shell
O     0.000000     0.000000     0.108400
H     0.000000     0.000000    -0.867200
```

⚠️ **Nota**: Para radicales, modificar en el código:
```python
bloque_xyz = f"* xyz 0 2\n{coords_str}\n*\n"  # Multiplicidad 2 (doblete)
```

---

## 📊 Interpretación de Resultados

### Energías
- **Negativas**: Normal en química cuántica (más negativo = más estable)
- **Unidades**: Hartree (1 Hartree = 27.211 eV = 627.5 kcal/mol)

### Frecuencias IR
- **Positivas**: Vibraciones reales
- **Negativas**: Estados de transición o geometría incorrecta
- **Intensidad alta**: Banda fuerte en espectro experimental

### Susceptibilidad Magnética
- **χ < -10**: Fuertemente diamagnética
- **-10 < χ < 0**: Diamagnética
- **χ > 0**: Paramagnética (electrones desapareados)

### Cargas Atómicas
- **Positiva**: Átomo deficiente en electrones
- **Negativa**: Átomo rico en electrones
- **Mulliken vs Löwdin**: Löwdin generalmente más estable

---



## 👨‍💻 Autor

Desarrollado por mi • Cálculos cuánticos con ORCA y PySCF

## 📄 Licencia

Este proyecto es de software propietario, se prohibe su distribución sin licencia.

---
