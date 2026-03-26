# PTSimSetup

A Python toolkit for calculating and analysing **proton transfer (PT) collective variables (CVs)** along a protein channel, designed for use in QM/MM enhanced-sampling simulations with [PLUMED](https://www.plumed.org/).

The toolkit calculates three CVs per frame:
- **CEC Projection** (`CV`): a geometric projection of the modified Centre of Excess Charge (mCEC) onto the vector between donor and acceptor sites.
- **Principal Curve Connectivity** (`PCC`): a connectivity measure along a principal curve fitted through the protonable heavy atom point cloud.
- **CEC-weighted PCC** (`CECPCC`): a variant of PCC where occupancies are screened by proximity to the CEC.

These CVs are described in [König et al., *J. Phys. Chem. A* 2006, 110(2), 548–563](https://pubs.acs.org/doi/10.1021/jp052328q) and [Li & Voth, *PNAS* 2021, 118(49), e2113141118](https://www.pnas.org/doi/full/10.1073/pnas.2113141118).

---

## Files

| File | Description |
|---|---|
| `pt_setup.py` | Main entry point. Iterates over trajectory frames, computes CVs, writes outputs. |
| `calcCEC.py` | Calculates the Centre of Excess Charge and its projection onto the donor–acceptor axis. |
| `calcPCC.py` | Calculates PCC and CEC-weighted PCC; fits a polynomial principal curve through the heavy atom positions. |
| `plumedWriter.py` | Writes a PLUMED input file (`.dat`) implementing the CEC CV for biasing in a simulation. |

---

## Dependencies

- [MDAnalysis](https://www.mdanalysis.org/)
- [NumPy](https://numpy.org/)
- [SciPy](https://scipy.org/)
- [tqdm](https://tqdm.github.io/)

Install with:
```bash
pip install MDAnalysis numpy scipy tqdm
```

---

## Usage

```bash
python pt_setup.py <reference.pdb> --initial <resid_donor> --target <resid_acceptor> [options]
```

### Required arguments

| Argument | Description |
|---|---|
| `reference_structure` | Reference PDB file |
| `--initial` | Residue ID of the proton donor site |
| `--target` | Residue ID of the proton acceptor site |

### Optional arguments

| Argument | Default | Description |
|---|---|---|
| `--trajectory` | `None` (uses reference PDB) | Trajectory file (`.xtc` or `.pdb`) |
| `--ndx` | `None` | GROMACS index file to prepend to written `.ndx` outputs |
| `--write` | `None` | Write PDB, NDX and PLUMED files for frames where `CV > x` |
| `--write_top` | `None` | Write PDB, NDX and PLUMED files for the top N frames by PCC |
| `--dynamic_nonwaters` | `False` | Dynamically include bridging residues in the QM region |
| `--ligand` | `None` | Residue name of a ligand to include with `--dynamic_nonwaters`, otherwise, no ligand is included. |
| `--interpolate_pc` | `True` | Fit a polynomial principal curve (cubic by default); if `False`, uses linear interpolation |
| `--cylinder_radius_factor` | `2` | Radius of the cylinder used to select water oxygens = donor–acceptor distance / this factor |
| `--verbose` | `False` | Print per-frame CV, time and atom count |
| `--debug` | `False` | Print detailed intermediate quantities and show a 3D debug plot |

### Example

```bash
python pt_setup.py system.pdb --initial 85 --target 212 --trajectory traj.xtc --write 0.4 --verbose
```

---

## Outputs

| File | Description |
|---|---|
| `cv_values.xvg` | Tab-separated: frame index, `CV`, `PCC`, `CECPCC`, dynamic residue list |
| `frame-N.pdb` | Full system PDB for frame N (if `--write` or `--write_top` used) |
| `frame-N.ndx` | GROMACS index file with `[QM-CEC-atoms]` and `[QM-atoms]` groups |
| `frame-N.plumed-cec.dat` | PLUMED input file implementing the CEC CV for frame N |

---

## Supported residue types

The CEC definition uses protonable sidechain heavy atoms and polar hydrogens. The following residue names are currently recognised:

**Donors/acceptors:** `ASP`, `GLU`, `HIS`, `LYS`, `TYR`  
**Protonated variants:** `ASPH`/`ASH`, `GLUH`/`GLH`, `HSD`/`HSE`/`HSP`/`HIP`

---

## Method notes

### Centre of Excess Charge (CEC)
Main reference: [König et al., *J. Phys. Chem. A* 2006, 110(2), 548–563](https://pubs.acs.org/doi/10.1021/jp052328q)
The CEC is computed as:

```
CEC = Σ r_Hi  −  Σ w_j * r_Xj  −  Σ_ij f_sw(d_ij) * (r_Hi − r_Xj)
```

where `r_Hi` are hydrogen positions, `r_Xj` are heavy atom positions with weight factors `w_j`, and `f_sw` is a sigmoidal switching function.

The CEC is then projected onto the donor–acceptor axis using a cosine projection formula involving the weighted distances `RD`, `RE`, and `RDE`, for a normalised progress coordinate.

### Principal Curve Connectivity (PCC)
Main reference: [Li & Voth, *PNAS* 2021, 118(49), e2113141118](https://www.pnas.org/doi/full/10.1073/pnas.2113141118)
A cubic polynomial is fitted to the protonable heavy atom positions as a principal curve between the donor and acceptor centres of mass. Each bead of the curve is assigned an occupancy derived from its water coordination number (via the same switching function). The PCC is the **geometric mean** of pairwise bead-to-bead connectivities and approaches 1 when a continuous hydrated chain exists. Note that the PCC does not consider dipolar alignment, so not all continuous hydrated chains will be ideal water wires.

In the CEC-weighted variant, occupancy contributions are raised to an exponent that decays with distance from the CEC, so beads near the CEC contribute more strongly. This is described in [doi:10.1073/pnas.2113141118](https://doi.org/10.1073/pnas.2113141118): 

---

## Notes and known limitations

- Currently assumes `chainID A` for protein atom selections; modify the selection strings in `pt_setup.py` if your system differs.
- The ligand CEC heavy atom weights in `pt_setup.py` are currently hardcoded (placeholder indices); parametrise these for your specific ligand.
- PLUMED output files use Ångström units (`UNITS LENGTH=A`); ensure your simulation also uses Å or adjust accordingly.
- The `--write_top` flag ranks frames by PCC (not `CV`); this is currently hardcoded and can be changed by editing `writeout_parameter` in `pt_setup.py`.