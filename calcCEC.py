import numpy as np
def calculate_cec_projection(cec_hydrogen_atoms, cec_heavy_atoms, initial_cv_atoms, target_cv_atoms, weight_factors, rsw, dsw, sigma, debug=False):    

    # get all hydrogen positions
    QMhydrogen_positions = cec_hydrogen_atoms.positions
    QMheavy_positions = cec_heavy_atoms.positions

    # do the above using nested loops
    QMhydrogen_to_QMheavy_vectors = np.zeros((len(cec_hydrogen_atoms), len(cec_heavy_atoms), 3))
    QMhydrogen_to_QMheavy_distances = np.zeros((len(cec_hydrogen_atoms), len(cec_heavy_atoms)))
    QMhydrogen_to_QMheavy_vectors_scaled = np.zeros((len(cec_hydrogen_atoms), len(cec_heavy_atoms), 3))
    for i in range(len(cec_hydrogen_atoms)):
        for j in range(len(cec_heavy_atoms)):
            QMhydrogen_to_QMheavy_vectors[i,j,:] = QMhydrogen_positions[i,:] - QMheavy_positions[j,:]
            QMhydrogen_to_QMheavy_distances[i,j] = np.linalg.norm(QMhydrogen_to_QMheavy_vectors[i,j,:])
            QMhydrogen_to_QMheavy_vectors_scaled[i,j,:] = QMhydrogen_to_QMheavy_vectors[i,j,:] * (1.0 / (1.0 + np.exp((QMhydrogen_to_QMheavy_distances[i,j] - rsw) / dsw)))

    # sum of all hydrogen vectors (A)
    QMhydrogen_sum_vector = np.sum(cec_hydrogen_atoms.positions, axis=0)
    # do the above using a loop
    QMhydrogen_sum_vector = np.zeros(3)
    for atom in cec_hydrogen_atoms:
        QMhydrogen_sum_vector += atom.position
    
    # weighted sum of all heavy atom vectors (B)
    QMheavy_weighted_sum_vector = np.zeros(3)
    for atom in cec_heavy_atoms:
        QMheavy_weighted_sum_vector += weight_factors[atom.index] * atom.position

    # sum of all scaled pairwise vectors (C)
    QMhydrogen_to_QMheavy_vectors_scaled_sum = np.sum(QMhydrogen_to_QMheavy_vectors_scaled, axis=(0,1))

    # center of excess charge = A - B - C
    CEC = (QMhydrogen_sum_vector - QMheavy_weighted_sum_vector - QMhydrogen_to_QMheavy_vectors_scaled_sum) #+ QMheavy_COG

    # project the CEC onto the vector between the initial and target sites
    # vector between initial and target sites
    initial_to_target_vector = np.mean(target_cv_atoms.positions, axis=0) - np.mean(initial_cv_atoms.positions, axis=0)
    # normalise this vector
    initial_to_target_vector = initial_to_target_vector / np.linalg.norm(initial_to_target_vector)

    # distances between all atoms in initial site and the CEC
    initial_site_to_CEC_distances = np.linalg.norm(initial_cv_atoms.positions - CEC, axis=1)

    # distances between all atoms in target site and the CEC
    target_site_to_CEC_distances = np.linalg.norm(target_cv_atoms.positions - CEC, axis=1)

    # sort both these arrays from smallest to largest
    initial_site_to_CEC_distances_sorted = np.sort(initial_site_to_CEC_distances)
    target_site_to_CEC_distances_sorted = np.sort(target_site_to_CEC_distances)

    # update the identities of the atoms in the initial and target sites to be the closest two to the CEC
    initial_cv_atoms = initial_cv_atoms[np.argsort(initial_site_to_CEC_distances)]
    target_cv_atoms = target_cv_atoms[np.argsort(target_site_to_CEC_distances)]
    
    # weight for each distance is 1 / (exp((RD1 - RD2)/sigma) + 1))
    w1 = 1.0 / (1.0 + np.exp((initial_site_to_CEC_distances_sorted[0] - initial_site_to_CEC_distances_sorted[1]) / sigma))
    w2 = 1.0 / (1.0 + np.exp((initial_site_to_CEC_distances_sorted[1] - initial_site_to_CEC_distances_sorted[0]) / sigma))
    w3 = 1.0 / (1.0 + np.exp((target_site_to_CEC_distances_sorted[0] - target_site_to_CEC_distances_sorted[1]) / sigma))
    w4 = 1.0 / (1.0 + np.exp((target_site_to_CEC_distances_sorted[1] - target_site_to_CEC_distances_sorted[0]) / sigma))

    # first two colvar terms 
    RD = w1*initial_site_to_CEC_distances_sorted[0] + w2*initial_site_to_CEC_distances_sorted[1]
    RE = w3*target_site_to_CEC_distances_sorted[0] + w4*target_site_to_CEC_distances_sorted[1]

    RD1E1 = np.linalg.norm(initial_cv_atoms[0].position - target_cv_atoms[0].position)
    RD1E2 = np.linalg.norm(initial_cv_atoms[0].position - target_cv_atoms[1].position)
    RD2E1 = np.linalg.norm(initial_cv_atoms[1].position - target_cv_atoms[0].position)
    RD2E2 = np.linalg.norm(initial_cv_atoms[1].position - target_cv_atoms[1].position)

    # cross term
    RDE = (w1*w3*RD1E1) + (w1*w4*RD1E2) + (w2*w3*RD2E1) + (w2*w4*RD2E2)

    # final CV (projection)
    CV = ((RD * RD) + (RDE * RDE) - (RE * RE)) / (2.0 * (RDE * RDE))

    # if debug flag is set to true, print out everything below
    if debug == True :

        # check atoms
        cec_atoms = cec_heavy_atoms + cec_hydrogen_atoms
        print('\n')
        for atom in cec_atoms:
            print('\tIndex:', atom.index, 'ID:', atom.id, ' resID:', atom.resid, 'resname', atom.resname)
        print('\n')

        # DEBUG PRINTOUTS
        #set print preferences to 3 dp
        np.set_printoptions(precision=3)
        print('\n\tinitial site COG:', np.mean(initial_cv_atoms.positions, axis=0))
        print('\ttarget site COG:', np.mean(target_cv_atoms.positions, axis=0))
        print('\n\tQM hydrogens:', len(cec_hydrogen_atoms), 'QM heavy:', len(cec_heavy_atoms))    
        print('\tmatrix shape:', QMhydrogen_to_QMheavy_vectors.shape)
        print('\tN cross terms:', len(cec_hydrogen_atoms) * len(cec_heavy_atoms))
        print('\n         H vector sum:', QMhydrogen_sum_vector)
        print('        wX vector sum:', QMheavy_weighted_sum_vector)
        print('           cross term:', QMhydrogen_to_QMheavy_vectors_scaled_sum)
        print('            CEC coord:', CEC)
        print('\n\tRD CALCULATION:', '(', '%.2f' % w1, '*', '%.2f' % initial_site_to_CEC_distances_sorted[0], ') + (', '%.2f' % w2,'*', '%.2f' % initial_site_to_CEC_distances_sorted[1], ')')
        print('\tRE CALCULATION:',   '(', '%.2f' % w3, '*', '%.2f' % target_site_to_CEC_distances_sorted[0],  ') + (', '%.2f' % w4,'*', '%.2f' % target_site_to_CEC_distances_sorted[1],  ')')
        print('\n\tRDE = (w1*w3*RD1E1) + (w1*w4*RD1E2) + (w2*w3*RD2E1) + (w2*w4*RD2E2)')
        print('\n\tw1: %.2f' % w1 ,'\n\tw2: %.2f' % w2 ,'\n\tw3: %.2f' % w3 ,'\n\tw4: %.2f' % w4)
        print('\n\tRD1E1: %.2f' % RD1E1 ,'\n\tRD1E2: %.2f' % RD1E2 ,'\n\tRD2E1: %.2f' % RD2E1 ,'\n\tRD2E2: %.2f' % RD2E2)
        print('\n\tRDE CALCULATION: (', '%.2f' % w1,'%.2f' % w3,'%.2f' % RD1E1, ') + (', '%.2f' % w1,'%.2f' % w4,'%.2f' % RD1E2, ') + (', '%.2f' % w2,'%.2f' % w3,'%.2f' % RD2E1, ') + (', '%.2f' % w2,'%.2f' % w4,'%.2f' % RD2E2, ')')
        print('\n\tRD:', '%.6f' % RD)
        print('\tRE:', '%.6f' % RE)
        print('\tRDE:', '%.6f' % RDE)
        print('\n')

        # DEBUGGING PLOT
        import matplotlib.pyplot as plt
        wf=[]
        for atom in cec_heavy_atoms:
            wf.append(weight_factors[atom.index])
        wf=np.array(wf)
        wf=wf/np.max(wf)
        ax = plt.axes(projection='3d')
        # plot the cec atoms
        ax.scatter3D(QMheavy_positions[:,0], QMheavy_positions[:,1], QMheavy_positions[:,2], color='red', marker='o', s=wf*100 )
        ax.scatter3D(QMhydrogen_positions[:,0], QMhydrogen_positions[:,1], QMhydrogen_positions[:,2], color='grey', marker='o')
        # now plot the initial and target site cec atoms
        ax.scatter3D(initial_cv_atoms.positions[:,0], initial_cv_atoms.positions[:,1], initial_cv_atoms.positions[:,2], color='yellow', marker='o')
        ax.scatter3D(target_cv_atoms.positions[:,0], target_cv_atoms.positions[:,1], target_cv_atoms.positions[:,2], color='yellow', marker='o')
        # label the initial and target sites at the average position of the atoms
        ax.text3D(np.mean(initial_cv_atoms.positions[:,0]), np.mean(initial_cv_atoms.positions[:,1]), np.mean(initial_cv_atoms.positions[:,2]), 'initial site', color='black')
        ax.text3D(np.mean(target_cv_atoms.positions[:,0]), np.mean(target_cv_atoms.positions[:,1]), np.mean(target_cv_atoms.positions[:,2]), 'target site', color='black')
        # connect the hydrogens to the heavy atoms with lines, colour their alpha based on the switching function
        for i in range(len(cec_hydrogen_atoms)):
            for j in range(len(cec_heavy_atoms)):
                ax.plot3D([QMhydrogen_positions[i,0], QMheavy_positions[j,0]], [QMhydrogen_positions[i,1], QMheavy_positions[j,1]], [QMhydrogen_positions[i,2], QMheavy_positions[j,2]], color='grey', alpha=(1.0 / (1.0 + np.exp((QMhydrogen_to_QMheavy_distances[i,j] - rsw) / dsw))))
        # annotate the heavy atoms with their weight factors
        for i in range(len(cec_heavy_atoms)):
            ax.text3D(QMheavy_positions[i,0], QMheavy_positions[i,1]+0.2, QMheavy_positions[i,2]+0.2, '%.1f' % weight_factors[cec_heavy_atoms[i].index], color='black')
        # plot the CEC
        ax.scatter3D(CEC[0], CEC[1], CEC[2], color='green', marker='o')
        ax.text3D(CEC[0], CEC[1], CEC[2], 'CEC', color='black')
        # plot line "RD"
        ax.plot3D([CEC[0], initial_cv_atoms.positions[0,0]], [CEC[1], initial_cv_atoms.positions[0,1]], [CEC[2], initial_cv_atoms.positions[0,2]], color='black', linestyle='--')
        ax.text3D((CEC[0] + initial_cv_atoms.positions[0,0]) / 2.0, (CEC[1] + initial_cv_atoms.positions[0,1]) / 2.0, (CEC[2] + initial_cv_atoms.positions[0,2]) / 2.0, 'RD', color='black')
        # plot line "RE"
        ax.plot3D([CEC[0], target_cv_atoms.positions[0,0]], [CEC[1], target_cv_atoms.positions[0,1]], [CEC[2], target_cv_atoms.positions[0,2]], color='black', linestyle='--')
        ax.text3D((CEC[0] + target_cv_atoms.positions[0,0]) / 2.0, (CEC[1] + target_cv_atoms.positions[0,1]) / 2.0, (CEC[2] + target_cv_atoms.positions[0,2]) / 2.0, 'RE', color='black')
        # plot line "RDE"
        ax.plot3D([initial_cv_atoms.positions[0,0], target_cv_atoms.positions[0,0]], [initial_cv_atoms.positions[0,1], target_cv_atoms.positions[0,1]], [initial_cv_atoms.positions[0,2], target_cv_atoms.positions[0,2]], color='black', linestyle='--')
        ax.text3D((initial_cv_atoms.positions[0,0] + target_cv_atoms.positions[0,0]) / 2.0, (initial_cv_atoms.positions[0,1] + target_cv_atoms.positions[0,1]) / 2.0, (initial_cv_atoms.positions[0,2] + target_cv_atoms.positions[0,2]) / 2.0, 'RDE', color='black')
        
        plt.show()

    # return Projection value and CEC position
    return CV, CEC
