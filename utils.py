import pandas as pd
import re
import numpy as np
from pyscf import gto, dft, scf

PYSCF_AVAILABLE = True

class Orca:
    def __init__(self, output_path):
        self.output_path = output_path
        try:
            with open(output_path, 'r', encoding='utf-8', errors='ignore') as f:
                self.content = f.read()
        except FileNotFoundError:
            raise FileNotFoundError(f"Output file not found at: {output_path}")

    @staticmethod
    def generate_input(xyz_content, calculation_type, method, basis_set, extra_keywords, calc_nmr=False):
        base_keywords = f"! {method} {basis_set} {extra_keywords}"

        if "zora" in basis_set.lower():
            base_keywords += " ZORA"

        if calculation_type == "Optimización de Geometría":
            calc_keywords = "OPT"
        elif calculation_type == "Frecuencias Vibracionales (IR)":
            calc_keywords = "OPT FREQ"
        else:
            calc_keywords = ""

        if calc_nmr:
            calc_keywords += " NMR"

        header = f"{base_keywords} {calc_keywords}\n"
        
        lines = xyz_content.strip().split('\n')
        num_atoms = 0
        start_idx = 0
        
        for i, line in enumerate(lines):
            try:
                num_atoms = int(line.strip())
                start_idx = i
                break
            except ValueError:
                continue

        coord_lines = []
        if num_atoms > 0:
            for line in lines[start_idx + 1:]:
                line_str = line.strip()
                if not line_str:
                    continue
                parts = line_str.split()
                if len(parts) >= 4:
                    try:
                        float(parts[1])
                        float(parts[2])
                        float(parts[3])
                        coord_lines.append(line_str)
                        if len(coord_lines) == num_atoms:
                            break
                    except ValueError:
                        pass

        coords_str = "\n".join(coord_lines)
        xyz_block = f"* xyz 0 1\n{coords_str}\n*\n"
        return header + xyz_block

    def check_convergence(self):
        return "THE OPTIMIZATION HAS CONVERGED" in self.content

    def extract_final_energy(self):
        matches = re.findall(r'FINAL SINGLE POINT ENERGY\s+([-\d.]+)', self.content)
        if matches:
            return float(matches[-1])
        return None

    def extract_optimized_geometry(self):
        pattern = r'CARTESIAN COORDINATES \(ANGSTROEM\)\s*\n\s*-+\s*\n((?:\s*\S+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s*\n)+)'
        matches = list(re.finditer(pattern, self.content))
        if not matches:
            return None

        coords_block = matches[-1].group(1).strip()
        coord_lines = [line.strip() for line in coords_block.split('\n') if line.strip()]
        if not coord_lines:
            return None

        num_atoms = len(coord_lines)
        xyz_block = f"{num_atoms}\nOptimized geometry extracted from {self.output_path}\n"
        for line in coord_lines:
            parts = line.split()
            if len(parts) >= 4:
                xyz_block += f"{parts[0]:<2} " + " ".join(f"{float(coord):>12.6f}" for coord in parts[1:4]) + "\n"
        return xyz_block

    def extract_ir_spectrum(self, scaling_factor=1.0):
        pattern = r'IR SPECTRUM\s*\n-+\n(?:.|\n)*?-+\n((?:.|\n)*?)(?=\n\s*\*|\n\s*-{2,}\n[A-Z]|\Z)'
        match = re.search(pattern, self.content)

        if not match:
            return pd.DataFrame()

        data = []
        data_block = match.group(1).strip()

        for line in data_block.split('\n'):
            parts = line.split()
            if len(parts) > 3 and parts[0].endswith(':'):
                try:
                    freq = float(parts[1])
                    intensity = float(parts[3])
                    if freq > 10.0:
                        data.append({"Frequency": freq * scaling_factor, "Intensity": intensity})
                except (ValueError, IndexError):
                    continue

        return pd.DataFrame(data)

    def extract_energy_components(self):
        patterns = {
            "Nuclear Repulsion": r'Nuclear Repulsion\s+:\s*([-\d.]+)',
            "Electronic Energy": r'Electronic Energy\s+:\s*([-\d.]+)',
            "One Electron Energy": r'One Electron Energy\s+:\s*([-\d.]+)',
            "Two Electron Energy": r'Two Electron Energy\s+:\s*([-\d.]+)',
        }
        energies = {}
        for name, pattern in patterns.items():
            matches = re.findall(pattern, self.content)
            if matches:
                energies[name] = [float(matches[-1])]

        return pd.DataFrame.from_dict(energies, orient='index', columns=['Energy (Hartree)']) if energies else None

    def extract_atomic_charges(self):
        charge_data = {}
        for charge_type in ['MULLIKEN', 'LOEWDIN']:
            pattern = re.compile(rf'{charge_type} ATOMIC CHARGES\s*\n-+\n((?:.|\n)*?)(?=\n\n|\Z)')

            matches = list(re.finditer(pattern, self.content))
            if matches:
                final_match = matches[-1]
                charges = []
                for line in final_match.group(1).strip().split('\n'):
                    parts = line.split()
                    if len(parts) == 4 and parts[2] == ':':
                        charges.append({"Atom": f"{parts[0]} {parts[1]}", "Charge": float(parts[3])})
                if charges:
                    charge_data[charge_type.capitalize()] = pd.DataFrame(charges)

        return charge_data if charge_data else None

    def extract_orbital_energies(self):
        pattern = re.compile(r'ORBITAL ENERGIES\s*\n-+\n((?:.|\n)*?)(?=\n\n|\Z|\*Only the first)')

        matches = list(re.finditer(pattern, self.content))
        if not matches:
            return None

        final_match = matches[-1]

        orbitals = []
        for line in final_match.group(1).strip().split('\n')[2:]:
            parts = line.split()
            if len(parts) == 4:
                orbitals.append({
                    "Number": int(parts[0]),
                    "Occupancy": float(parts[1]),
                    "Energy (Eh)": float(parts[2]),
                    "Energy (eV)": float(parts[3])
                })

        return pd.DataFrame(orbitals) if orbitals else None

    def extract_reduced_orbital_charges(self):
        charge_data = {}
        for charge_type in ['MULLIKEN', 'LOEWDIN']:
            pattern = re.compile(
                rf'{charge_type} REDUCED ORBITAL CHARGES\s*\n-+\n((?:.|\n)*?)(?=\n\n|\Z|\s*\*+\n|\s*-{{2,}}\n[A-Z])')

            matches = list(re.finditer(pattern, self.content))
            if matches:
                final_match = matches[-1]
                orbital_charges = []
                current_atom = ""

                text_block = final_match.group(1).strip()
                for line in text_block.split('\n'):
                    clean_line = line.strip()
                    if not clean_line:
                        continue

                    initial_parts = line.split()
                    if initial_parts and initial_parts[0].isdigit():
                        if len(initial_parts) > 1:
                            current_atom = f"{initial_parts[0]} {initial_parts[1]}"

                    orbital_charge_pairs = re.findall(r'([a-zA-Z0-9]+)\s*:\s*([\d.-]+)', line)

                    for orbital, charge in orbital_charge_pairs:
                        if orbital in ['s', 'p', 'd', 'f'] and current_atom:
                            try:
                                orbital_charges.append({
                                    "Atom": current_atom,
                                    "Orbital": orbital,
                                    "Charge": float(charge)
                                })
                            except ValueError:
                                continue

                if orbital_charges:
                    charge_data[charge_type.capitalize()] = pd.DataFrame(orbital_charges)

        return charge_data if charge_data else None

    def extract_nmr_data(self):
        block_pattern = re.compile(
            r'CHEMICAL SHIELDING SUMMARY \(ppm\)\s*\n-+\n\n((?:.|\n)*?)(?=\n\n\s*NMR shielding tensor|\Z|\n\s*-{2,}\n)')
        match = re.search(block_pattern, self.content)

        if not match:
            return None

        text_block = match.group(1).strip()
        lines = text_block.split('\n')

        if len(lines) <= 2:
            return None

        nmr_data = []
        for line in lines[2:]:
            parts = line.split()
            if len(parts) == 4:
                try:
                    nmr_data.append({
                        "Nucleus": int(parts[0]),
                        "Element": parts[1],
                        "Isotropic (ppm)": float(parts[2]),
                        "Anisotropy (ppm)": float(parts[3])
                    })
                except ValueError:
                    continue

        if nmr_data:
            return pd.DataFrame(nmr_data)

        return None


class PySCFCalculator:

    @staticmethod
    def calculate_susceptibility(xyz_content, method='b3lyp', basis='def2svp'):
        try:
            lines = [l.strip() for l in xyz_content.strip().split('\n') if l.strip()]

            num_atoms = 0
            start_idx = 0
            for i, line in enumerate(lines):
                try:
                    num_atoms = int(line.strip())
                    start_idx = i
                    break
                except ValueError:
                    continue

            if num_atoms == 0:
                return {"error": "Invalid XYZ format: could not find number of atoms"}

            coord_lines = []
            for line in lines[start_idx + 1:]:
                line_str = line.strip()
                if not line_str:
                    continue
                parts = line_str.split()
                if len(parts) >= 4:
                    try:
                        float(parts[1])
                        float(parts[2])
                        float(parts[3])
                        coord_lines.append(line_str)
                        if len(coord_lines) == num_atoms:
                            break
                    except ValueError:
                        pass

            if len(coord_lines) != num_atoms:
                return {"error": "Invalid XYZ format: could not parse coordinates"}

            atom_str = ""
            for idx, line in enumerate(coord_lines):
                parts = line.split()
                if len(parts) < 4:
                    return {"error": f"Invalid line {start_idx + 3 + idx}: {line}"}
                try:
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                    atom_str += f"{parts[0]} {x} {y} {z}; "
                except ValueError:
                    return {"error": f"Invalid coordinates in line {start_idx + 3 + idx}: {line}"}

            basis_map = {
                'def2-svp': 'def2svp',
                'def2-tzvp': 'def2tzvp',
                '6-31+g(d,p)': '6-31+g*',
                '6-311++g(d,p)': '6-311++g**',
                'cc-pvdz': 'ccpvdz'
            }
            pyscf_basis = basis_map.get(basis.lower(), basis.lower())

            mol = gto.M(
                atom=atom_str,
                basis=pyscf_basis,
                unit='Angstrom'
            )

            mf = dft.RKS(mol)
            mf.xc = method.lower()

            energy = mf.kernel()

            if not mf.converged:
                return {"error": "SCF did not converge in PySCF"}

            coords = mol.atom_coords()
            charges = mol.atom_charges()

            total_mass = sum(charges)
            com = np.sum(coords * charges[:, np.newaxis], axis=0) / total_mass if total_mass > 0 else np.zeros(3)
            chi_tensor = np.zeros((3, 3))

            for i, (coord, Z) in enumerate(zip(coords, charges)):
                r = coord - com
                r2 = np.dot(r, r)

                for j in range(3):
                    for k in range(3):
                        if j == k:
                            chi_tensor[j, k] -= Z * (r2 - r[j] ** 2) / 6.0
                        else:
                            chi_tensor[j, k] -= Z * r[j] * r[k] / 6.0

            chi_iso = np.trace(chi_tensor) / 3.0
            chi_cgs = chi_iso * 0.78910

            magnetism_type = "Diamagnetic" if chi_cgs < 0 else "Paramagnetic"

            return {
                "tensor": chi_tensor.tolist(),
                "isotropic_au": float(chi_iso),
                "isotropic_cgs": float(chi_cgs),
                "type": magnetism_type,
                "scf_energy": float(energy),
                "converged": True,
                "calculation_method": "Pascal's approximation (diamagnetic)",
                "note": "Approximate calculation based on molecular geometry. For precise results use ORCA with NMR keywords."
            }

        except ImportError as e:
            return {
                "error": f"PySCF is not correctly installed: {str(e)}\nTry: pip install --upgrade pyscf"
            }
        except Exception as e:
            import traceback
            return {"error": f"Error in calculation: {str(e)}\n\nDetails:\n{traceback.format_exc()}"}