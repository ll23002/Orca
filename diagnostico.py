import pyscf

print("=== Diagnóstico de PySCF ===")
print(f"Versión: {pyscf.__version__}")
print("\nMódulos en pyscf.prop:")
try:
    import pyscf.prop
    print(dir(pyscf.prop))
except Exception as e:
    print(f"Error: {e}")

print("\n¿Existe pyscf.prop.nmr?")
try:
    from pyscf.prop import nmr
    print("✓ SÍ existe")
    print(f"Funciones en nmr: {dir(nmr)}")
except Exception as e:
    print(f"✗ NO existe: {e}")

print("\n¿Existe pyscf.prop.nmr.rhf?")
try:
    from pyscf.prop.nmr import rhf
    print("✓ SÍ existe")
    print(f"Funciones en rhf: {[x for x in dir(rhf) if not x.startswith('_')]}")
except Exception as e:
    print(f"✗ NO existe: {e}")