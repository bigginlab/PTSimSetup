import warnings
warnings.simplefilter("ignore")
import argparse
from tqdm import tqdm

import MDAnalysis as mda
import numpy as np

# import the plumed writer function
from plumedWriter import write_plumed_file
from calcCEC import calculate_cec_projection
from calcPCC import calculate_pcc
from calcPCC import calculate_cec_weighted_pcc
from calcPCC import fit_polynomial_to_point_cloud

rsw=1.25    # switching function parameters for heavy atoms (currently in A)
dsw=0.04    # switching function parameters for heavy atoms (currently in A)
sigma=0.025 # parameter for switching functionfor CV (currently in A)


# Argument parsing using argparse
def parse_args():
    parser = argparse.ArgumentParser(
        description='PTSimSetup: Calculate proton transfer CVs between initial and target residues.'
    )
    parser.add_argument('reference_structure', type=str, help='Reference structure (PDB file)')
    parser.add_argument('--initial', type=str, required=True, help='Initial site residue ID')
    parser.add_argument('--target', type=str, required=True, help='Target site residue ID')
    parser.add_argument('--trajectory', type=str, default=None, help='Trajectory file (.xtc or .pdb, optional)')
    parser.add_argument('--ndx', type=str, default=None, help='Input index file to modify when using write (optional)')
    parser.add_argument('--write_top', type=int, default=None, help='Write pdb and ndx files for the top N frames in terms of CV value (optional)')
    parser.add_argument('--dynamic_nonwaters', action='store_true', help='Include residues connecting the initial and target sites in the QM region (optional)')
    parser.add_argument('--verbose', action='store_true', help='Print out frame number, CV, time and number of atoms (optional)')
    parser.add_argument('--write', type=float, default=None, help='Write every frame where cv > x (optional)')
    parser.add_argument('--debug', action='store_true', help='Print out debug information (optional)')
    parser.add_argument('--cylinder_radius_factor', type=float, default=2, help='Cylinder radius for selecting waters is the distance between initial and target sites divided by this factor (optional, default: 2)')  # This is a placeholder for future use
    parser.add_argument('--interpolate_pc', type=bool, default=True, help='Polynomial fit of principal curve to CEC atoms rather than linear interpolation between initial and target COMs (optional, default: True)')
    parser.add_argument('--ligand', type=str, default=None, help='Residue name of ligand to include in QM region with dynamic nonwaters (optional, default: None)')
    #parser.add_argument('--ligand_cec_heavy_weights', type=dict, default=None, help='include ligand heavy atoms in the CEC definition, using these weights. Dict matches ligand heavy atom indices to weights (optional, default: None)')   
    #parser.add_argument('--ligand_cec_hydrogens', type=dict, default=None, help='include these ligand hydrogen atoms (optional, default: None)')     
    return parser.parse_args()

args = parse_args()

# raise error if reference structure is not a pdb file or a gro file
if args.reference_structure[-4:] != '.pdb':
    print('reference structure must be a pdb file, i havent tested this with gro files')
    exit()

# get arguments (reference structure, initial site, target site)
reference_structure = args.reference_structure
initial_site_resid = args.initial
target_site_resid = args.target
cylinder_radius_factor = args.cylinder_radius_factor

trajectory = args.trajectory if args.trajectory is not None else reference_structure
ndx = args.ndx
dynamic_nonwaters = args.dynamic_nonwaters
verbose = args.verbose
write_top = args.write_top is not None
write = args.write is not None
debug = args.debug
if write_top:
    N = args.write_top
ligand_resname = args.ligand

# objects store info
frame_cv_values = {}
#frame_water_densities = {}
frame_pcc_values = {}
frame_cecpcc_values = {}
frame_cec_atoms = {}
frame_noncec_qm_atoms = {}
frame_dynamic_resids = {}
frame_cec_heavyatoms = {}
frame_cec_hydrogens = {}
frame_weightfactors = {}

# central MDA universe object
u = mda.Universe(reference_structure, trajectory)

# dictinary of residue names and associated hydrogen selection tokens
residue_polar_hydrogens = {'ASP': ['HD*'], 
                           'GLU': ['HE*'], 
                           'HIS': ['HD1', 'HE2'], 
                           'LYS': ['HZ1', 'HZ2', 'HZ3'], 
                           'TYR': ['HH'], 
                           'ASPH': ['HD*'],
                           'ASH': ['HD*'],  
                           'GLUH': ['HE*'], 
                           'GLH': ['HE*'], 
                           'HSD': ['HD1'], 
                           'HSE': ['HE1'], 
                           'HSP': ['HD1', 'HE2'],
                           'HIP': ['HD1', 'HE2']}
                           #'LIG': ['HO', 'HN']} # H1-6

# if ligand being considered, append the ligand heavy atom names to the dicts down here

formatted_strings = []

for key, values in residue_polar_hydrogens.items():
    values_str = " or name ".join(values)
    formatted_str = f"(resname {key} and ( name {values_str} ))"
    formatted_strings.append(formatted_str)
    
## get relevant hydrogen selection tokens from this list of formatted strings
hydrogen_names = " or ".join(formatted_strings)

resudie_heavy_atoms = {'ASP': ['OD1', 'OD2'],
                       'GLU': ['OE1', 'OE2'],
                       'HIS': ['ND1', 'NE2'],
                       'LYS': ['NZ'],
                       'TYR': ['OH'],
                       'ASPH': ['OD1', 'OD2'],
                       'ASH': ['OD1', 'OD2'],
                       'GLUH': ['OE1', 'OE2'],
                       'GLH': ['OE1', 'OE2'],
                       'HSD': ['ND1'],
                       'HSE': ['NE2'],
                       'HSP': ['ND1', 'NE2'],
                       'HIP': ['ND1', 'NE2']}
                       #'LIG': ['O*', 'N*']} # O1, O2, OXT1, OXT2, N1, N2

# if ligand being considered, append the ligand heavy atom names to the dicts down here

formatted_strings = []

for key, values in resudie_heavy_atoms.items():
    values_str = " or name ".join(values)
    formatted_str = f"(resname {key} and ( name {values_str} ))"
    formatted_strings.append(formatted_str)

## get relevant heavy atom selection tokens from this list of formatted strings
heavy_atom_names = " or ".join(formatted_strings)

# temporary ligand defs
ligand_cec_heavy_weights = { # dict of indices mapped to weights
    61249: 0.0,
    61250: 0.0,
    61265: 0.0,
    61266: 0.0,
    61251: 3.0,
    61260: 3.0,
}

ligand_cec_hydrogens = [ # list of atom indices
61252, 61253, 61254, 61261, 61262, 61263
]
 
# Main Loop 
for ts in tqdm(u.trajectory):

    ### Selection logic goes here:

    ## select relevant (protonable) heavy atoms in initial and target sites
    initial_cv_atoms = u.select_atoms('chainID A and protein and not backbone and ((resid %s) and (%s))' % (initial_site_resid, heavy_atom_names))
    target_cv_atoms = u.select_atoms('chainID A and protein and not backbone and ((resid %s) and (%s))' % (target_site_resid, heavy_atom_names))

    # get distance between initial and target sites
    initial_to_target_distance = np.linalg.norm(np.mean(target_cv_atoms.positions, axis=0) - np.mean(initial_cv_atoms.positions, axis=0))
    initial_to_target_vector = np.mean(target_cv_atoms.positions, axis=0) - np.mean(initial_cv_atoms.positions, axis=0)
    
    # define the parameters for the cylinder
    zmax = initial_to_target_distance/2 # cylidner extends +/- zmax from the midpoint
    zmin = zmax
    r = zmax / cylinder_radius_factor # need a more clever way of determining this
    volume = np.pi * r**2 * initial_to_target_distance

    # Rotate the system such that initial to target vector is aligned along z axis (this way the cylindrical selection will connect the two sites)
    # initial to target vector post rotation should be (0, 0, initial_to_target_distance)

    # get the angle and axis of rotation
    angle = np.arccos(np.dot(initial_to_target_vector, np.array([0, 0, 1])) / (np.linalg.norm(initial_to_target_vector) * np.linalg.norm(np.array([0, 0, 1]))))
    axis = np.cross(initial_to_target_vector, np.array([0, 0, 1]))
    angle = angle * 180 / np.pi

    # rotate the system
    u.atoms.rotateby(angle, axis)
    initial_to_target_vector_post_rot = np.mean(target_cv_atoms.positions, axis=0) - np.mean(initial_cv_atoms.positions, axis=0)

    # dynamically select the CEC atoms
    cec_heavy_atoms = u.select_atoms('chainID A and protein and not backbone and ( (resid %s or resid %s) and (%s))' % (initial_site_resid, target_site_resid, heavy_atom_names)) \
    + u.select_atoms('(name O* and (resname TIP3 or resname SOL) and cyzone %s %s -%s ((resid %s or resid %s) and protein))' % (r, zmax, zmin, initial_site_resid, target_site_resid)) 
    cec_hydrogen_atoms = u.select_atoms('chainID A and protein and not backbone and ((resid %s or resid %s) and (%s))' % (initial_site_resid, target_site_resid, hydrogen_names)) \
    + u.select_atoms('(byres (name O* and (resname TIP3 or resname SOL) and cyzone %s %s -%s ((resid %s or resid %s) and protein))) and name H*' % (r, zmax, zmin, initial_site_resid, target_site_resid)) 

    # routine to add additional residues that connect the initial and target sites
    # will add any ASP, GLU, LYS or HIS residues in the cylinder to the selection
    if dynamic_nonwaters == True:
        QMdynamic = u.select_atoms('( (resname ASP or resname GLU or resname LYS or resname HIS or resname %s) and cyzone %s %s -%s ((resid %s or resid %s) and protein) ) and not (resid %s or resid %s)' % (ligand_resname, r, zmax, zmin, initial_site_resid, target_site_resid, initial_site_resid, target_site_resid))
        if len(QMdynamic) > 0:
            # get unique residues
            QMdynamic_resids_list = np.unique([atom.resid for atom in QMdynamic])
            # add the corresponding selection to cec_heavy_atoms
            cec_heavy_atoms += u.select_atoms('chainID A and protein and not backbone and ( (resid %s) and (%s))' % (' or resid '.join([str(resid) for resid in QMdynamic_resids_list]), heavy_atom_names))
            cec_hydrogen_atoms += u.select_atoms('chainID A and protein and not backbone and ( (resid %s) and (%s))' % (' or resid '.join([str(resid) for resid in QMdynamic_resids_list]), hydrogen_names))
        else:
            QMdynamic_resids_list = []

    # undo the rotation
    u.atoms.rotateby(-angle, axis)

    # calc density of cec_heavy_atoms in the cylinder (waters, basically)
    heavy_atom_density = len(cec_heavy_atoms) / volume

    # cec_atoms
    cec_atoms = cec_heavy_atoms + cec_hydrogen_atoms

    # print selections if debug
    if debug:
        # print in format resname+resid
        print('heavy atoms:', [atom.resname + str(atom.resid) + ' ' + str(atom.name) for atom in cec_heavy_atoms])
        print('hydrogens:', [atom.resname + str(atom.resid) + ' ' + str(atom.name) for atom in cec_hydrogen_atoms], '\n')
        print('initial site:', [atom.resname + str(atom.resid) + ' ' + str(atom.name) for atom in initial_cv_atoms])
        print('target site:', [atom.resname + str(atom.resid) + ' ' + str(atom.name) for atom in target_cv_atoms], '\n')
        print('initial to target distance:', initial_to_target_distance)
        print('initial to target vector:', initial_to_target_vector)
        print('initial to target vector post rotation:', initial_to_target_vector_post_rot, '\n')
        print('zmax:', zmax, 'zmin:', zmin, '\n')
        if dynamic_nonwaters == True:
            print('dynamic non-waters included:', [atom.resname + str(atom.resid) + ' ' + str(atom.name) for atom in QMdynamic], '\n')

    # raise error if one of the sites is not recognised (projection CV should work as long as there are 2 heavy atoms in each selection)
    if initial_cv_atoms.residues[0].resname not in ('ASP', 'ASPH', 'GLU', 'GLUH', 'ASH', 'HSP', 'TYR', 'HIP', 'HIS'):
        print('initial site must be an ASP or GLU or HIS residue (currently only ASP and GLU are supported). Detected: ', initial_cv_atoms.residues[0].resname)
        exit()
    if target_cv_atoms.residues[0].resname not in ('ASP', 'ASPH', 'GLU', 'GLUH', 'ASH', 'HSP', 'TYR', 'HIP', 'HIS'):
        print('target site must be an ASP or GLU or HIS residue (currently only ASP and GLU are supported). Detected: ', target_cv_atoms.residues[0].resname)
        exit()

    ### Selection logic ends here

    # logic for the ligand weight factors must go here
    # N1, N2 = 3.0
    # O* = 0

    # heavy atom weight factor is equal to the coordination number of hydrogens in the deprotonated state - maybe the weight should be one for hte protonated oxygen in asp
    weight_factors = {}
    for atom in cec_heavy_atoms:
        if atom.name == 'OW':                           # water
            weight_factors[atom.index] = 2.0
        elif atom.name == 'OD2' or atom.name == 'OE2':  # aspartate oxygens
            weight_factors[atom.index] = 0.0
        elif atom.name == 'OD1' or atom.name == 'OE1':  # glutamate oxygens
            weight_factors[atom.index] = 0.0
        elif atom.name == 'NZ':                         # lysine
            weight_factors[atom.index] = 3.0               
        elif atom.name == 'OH':                         # tyrosine
            weight_factors[atom.index] = 1.0
        elif atom.name == 'OH2':                        # water (alternate)
            weight_factors[atom.index] = 2.0
        elif atom.name == 'NE2':                        # histidine epsilon
            weight_factors[atom.index] = 0.5
        elif atom.name == 'ND1':                        # histidine delta
            weight_factors[atom.index] = 0.5

        else:
            print('atom name not recognized')
            print(atom.name, atom.resname, atom.resid)
            exit()
    
    # routine the select the ligand if present
    ligand_selection = u.select_atoms("resname %s" % ligand_resname)
    ligand_resid = ligand_selection.residues.resids
    if ligand_resid in QMdynamic_resids_list:
        ligand_in_qm_region = True
        # append the indices from the ligand_cec_heavy_weights dict and hydrogens list
        ligand_cec_heavyatoms = list(ligand_cec_heavy_weights.keys())
        cec_heavy_atoms += u.select_atoms(' or '.join(f'index {i}' for i in ligand_cec_heavyatoms))
        cec_hydrogen_atoms += u.select_atoms(' or '.join(f'index {i}' for i in ligand_cec_hydrogens))
        #print('ligand weight factors:', ligand_cec_heavy_weights, '\n')

        # merge the dicts
        weight_factors.update(ligand_cec_heavy_weights)
        # redeclare cec_atoms
        cec_atoms = cec_heavy_atoms + cec_hydrogen_atoms

    # if debug print out the weight factors
    if debug:
        print('weight factors:', weight_factors, '\n')

    CV, CEC = calculate_cec_projection(cec_heavy_atoms=cec_heavy_atoms, cec_hydrogen_atoms=cec_hydrogen_atoms, 
                                  initial_cv_atoms=initial_cv_atoms, target_cv_atoms=target_cv_atoms, 
                                  weight_factors=weight_factors, rsw=rsw, dsw=dsw, sigma=sigma)

    # define the principal curve for the PCC CV by fitting a polynomial to the point cloud of CEC atoms
    n_beads = round(initial_to_target_distance) # n_beads defined such that initial_to target distance is divided into 1A segments
    
    if args.interpolate_pc:
        principal_curve = fit_polynomial_to_point_cloud(cec_heavy_atoms.positions, (initial_cv_atoms.center_of_mass(), target_cv_atoms.center_of_mass()), n_beads, degree=3)
    else:
        # if not interpolating, just use linear interpolation between initial and target COMs
        principal_curve = np.linspace(initial_cv_atoms.center_of_mass(), target_cv_atoms.center_of_mass(), n_beads)

    # principal curve connectivity CVs from https://doi.org/10.1073/pnas.2113141118
    PCC = calculate_pcc(cec_heavy_atoms=cec_heavy_atoms, principal_curve=principal_curve)
    CECPCC = calculate_cec_weighted_pcc(cec_heavy_atoms=cec_heavy_atoms, principal_curve=principal_curve, CEC=CEC)

    # printout frame by frame information if trajectory is passed
    if trajectory[-4:] == '.xtc':
        if verbose == True:
            print('frame:', ts.frame, '\tCV: %.6f' % CV, '\ttime:', ts.time, '\tnumber of atoms:', len(cec_heavy_atoms) + len(cec_hydrogen_atoms), end='\r')
    else:
        print('\n\tCV: %.6f' % CV, '\tnumber of atoms:', len(cec_heavy_atoms) + len(cec_hydrogen_atoms), end='\n\r') 

    # Double check that these selection tokens work for your specific case - they should work for amino acid sidechains
    initial_qm_atoms_noncec = u.select_atoms('chainID A and protein and not (backbone or name HN or name HA or name H) and ((resid %s))' % (initial_site_resid))
    target_qm_atoms_noncec = u.select_atoms('chainID A and protein and not (backbone or name HN or name HA or name H) and ((resid %s))' % (target_site_resid))
    
    if dynamic_nonwaters == True:
        if len(QMdynamic_resids_list) > 0 :
            resid_clause = 'resid ' + ' or resid '.join([str(resid) for resid in QMdynamic_resids_list])
            #selection_str = f'chainID A and protein and not (backbone or name HN or name HA or name H) and ( ({resid_clause}) )'
            selection_str = (
                f"(resname {ligand_resname}) or "
                f"(chainID A and protein and not (backbone or name HN or name HA or name H) and ({resid_clause}))"
            )
      
        else:
            # This ensures the selection is always syntactically valid but matches nothing
            selection_str = 'resname DUM'

        other_qm_atoms_noncec = u.select_atoms(selection_str)
    else:
        QMdynamic_resids_list = []
        # return an empty MDA selectiom (for concatenation later if the flag isnt used)
        other_qm_atoms_noncec = u.select_atoms('resname DUM')

    # if flag -w is passed, write pdb and ndx file for frames where CV > argument
    if write:
        if CV > args.write:

           # write pdb file for the frame in question
            with mda.Writer('frame-%s.pdb' % ts.frame, n_atoms=len(u.atoms)) as W:
                W.write(u.atoms)

            qm_region_atoms = cec_atoms + initial_qm_atoms_noncec + target_qm_atoms_noncec + other_qm_atoms_noncec # the QM region itself has to be defined here (we have only been dealing with CEC atoms)

            # write a gromacs index file with all the QM CEC atoms
            with open('frame-%s.ndx' % ts.frame, 'w') as f:

                # first write out the contents of the input index file, if specified
                if ndx is not None:
                    with open(ndx, 'r') as g:
                        for line in g:
                            f.write(line)

                f.write('[QM-CEC-atoms]\n')
                # must sort so they are ascendiung
                QMlist = np.sort(np.array(cec_atoms.atoms.indices))
                for atom in QMlist:
                    f.write('%s ' % int(atom + 1))
                f.write('\n\n')

                f.write('[QM-atoms]\n')
                # must sort so they are ascendiung and unique (double counting of initial and target CEC atoms)
                QMlist = np.sort(np.unique(np.array(qm_region_atoms.atoms.indices)))
                for atom in QMlist:
                    f.write('%s ' % int(atom + 1))
                f.write('\n\n')

            # write the plumed file
            write_plumed_file(cec_heavy_atoms=cec_heavy_atoms, cec_hydrogen_atoms=cec_hydrogen_atoms, 
                              initial_cv_atoms=initial_cv_atoms, target_cv_atoms=target_cv_atoms, 
                              weight_factors=weight_factors, 
                              outputfile='frame-%s.plumed-cec.dat' % ts.frame)

    # append the ts and proton transfer score to the dict
    frame_cv_values[ts.frame] = CV
    #frame_water_densities[ts.frame] = heavy_atom_density # written to cv file
    frame_pcc_values[ts.frame] = PCC # written to cv file
    frame_cecpcc_values[ts.frame] = CECPCC # written to cv file
    # required for after the fact frame writeout in write_top
    frame_cec_atoms[ts.frame] = cec_atoms
    frame_noncec_qm_atoms[ts.frame] = other_qm_atoms_noncec
    frame_dynamic_resids[ts.frame] = QMdynamic_resids_list
    # for the plumed file
    frame_cec_heavyatoms[ts.frame] = cec_heavy_atoms
    frame_cec_hydrogens[ts.frame] = cec_hydrogen_atoms
    frame_weightfactors[ts.frame] = weight_factors

############################################################################################################################

# routine to write out top N frames
if write_top:
    # Sort frame indices by descending PCC value
    writeout_parameter = frame_pcc_values
    N = min(N, len(writeout_parameter))

    sorted_keys = sorted(writeout_parameter, key=writeout_parameter.get, reverse=True)[:N]

    for frame_index in tqdm(sorted_keys):
        u.trajectory[frame_index]

        # Write PDB file
        with mda.Writer(f'frame-{frame_index}.pdb', n_atoms=len(u.atoms)) as W:
            W.write(u.atoms)

        # Select QM atoms
        initial_qm_atoms_noncec = u.select_atoms(
            f'chainID A and protein and not (backbone or name HN or name HA or name H) and (resid {initial_site_resid})'
        )
        target_qm_atoms_noncec = u.select_atoms(
            f'chainID A and protein and not (backbone or name HN or name HA or name H) and (resid {target_site_resid})'
        )
        other_qm_atoms_noncec = frame_noncec_qm_atoms[frame_index] if dynamic_nonwaters else u.select_atoms('')

        cec_atoms = frame_cec_atoms[frame_index]
        qm_region_atoms = cec_atoms + initial_qm_atoms_noncec + target_qm_atoms_noncec + other_qm_atoms_noncec

        # Write GROMACS index file
        with open(f'frame-{frame_index}.ndx', 'w') as f:
            if ndx is not None:
                with open(ndx, 'r') as g:
                    f.writelines(g.readlines())

            f.write('[QM-CEC-atoms]\n')
            QMlist = np.sort(cec_atoms.atoms.indices)
            f.write(' '.join(str(int(i + 1)) for i in QMlist) + '\n\n')

            f.write('[QM-atoms]\n')
            QMlist = np.sort(np.unique(qm_region_atoms.atoms.indices))
            f.write(' '.join(str(int(i + 1)) for i in QMlist) + '\n\n')

        # Write PLUMED input file
        write_plumed_file(
            cec_heavy_atoms=frame_cec_heavyatoms[frame_index],
            cec_hydrogen_atoms=frame_cec_hydrogens[frame_index],
            initial_cv_atoms=initial_cv_atoms,
            target_cv_atoms=target_cv_atoms,
            weight_factors=frame_weightfactors[frame_index],
            outputfile=f'frame-{frame_index}.plumed-cec.dat'
        )

# write out the cv_values file for all frames (faster than writing out each frame)
with open('cv_values.xvg', 'w') as f:
    for key in frame_cv_values.keys():
        f.write('%s %.6f %.6f %.6f %s\n' % (key, frame_cv_values[key], frame_pcc_values[key], frame_cecpcc_values[key], frame_dynamic_resids[key]))
