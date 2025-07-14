# WRITE OUT THE PLUMED FILE ############################################################################

def write_plumed_file(cec_hydrogen_atoms, cec_heavy_atoms, initial_cv_atoms, target_cv_atoms, weight_factors, sigma=0.025, rsw=1.25, dsw=0.04, outputfile='plumed-cec.dat', debug=False):
    
    # docstrings
    """
    Write out the plumed file for the CEC colvar.
    Parameters
    ----------
    cec_hydrogen_atoms : MDanalysis selection
        MDanalysis selection of QM hydrogens.
    cec_heavy_atoms : MDanalysis selection
        MDanalysis selection of QM heavy atoms.
    initial_cv_atoms : MDanalysis selection
        List of initial CV atoms (usually the protonated O or N).
    target_cv_atoms : MDanalysis selection
        List of target CV atoms (usually the protonable O or N).
    weight_factors : dict
        Dictionary of weight factors for the heavy atoms.
    sigma : float, optional
        Width of the Gaussian in the metadynamics bias. The default is 0.025.
    rsw : float, optional
        Width of the switching function. The default is 1.25.
    dsw : float, optional
        Width of the switching function. The default is 0.04.  
    outputfile : str, optional
        Name of the output file. The default is 'plumed-cec.dat'.
    """

    cec_atoms = cec_heavy_atoms + cec_hydrogen_atoms

    with open(outputfile, 'w') as f:

        # UNITS
        f.write('#units\n\n')
        f.write('UNITS LENGTH=A\n\n')

        # CELL
        f.write('#cell\n\n')
        f.write('cell: CELL\n\n')

        #f.write('#MOLINFO STRUCTURE=%s\n\n' % reference_structure)

        # IDs are non-unique so we need to use indices (position in topology)
        # we substitute the IDs with the indices. Because PLUMED is 1-indexed we also add 1
        cec_heavy_atoms.atoms.ids = cec_heavy_atoms.atoms.indices + 1
        cec_hydrogen_atoms.atoms.ids = cec_hydrogen_atoms.atoms.indices + 1
        initial_cv_atoms.atoms.ids = initial_cv_atoms.atoms.indices + 1
        target_cv_atoms.atoms.ids = target_cv_atoms.atoms.indices + 1 

        # FIX PBC ISSUES
        f.write('#pbc\n\n')
        f.write('WHOLEMOLECULES ENTITY0=')
        for atom in cec_hydrogen_atoms:
                f.write('%s,' % atom.id)
        for atom in cec_heavy_atoms:
                if atom == cec_heavy_atoms[-1]:
                    f.write('%s\n\n' % atom.id)
                else:
                    f.write('%s,' % atom.id)

        # BASIC VARIABLES
        f.write('#basic variables\n\n')

        # COM of all heavy atoms # fixed atom at this position?
        f.write('#center of mass of all qm atoms\n')
        f.write('qm_com: CENTER ATOMS=')
        for atom in cec_heavy_atoms:
            if atom == cec_heavy_atoms[-1]:
                f.write('%s' % atom.id)
            else:
                f.write('%s,' % atom.id)
        f.write('\n\n')
        f.write('position_qm_com: POSITION ATOM=qm_com NOPBC\n\n')

        # HYDROGEN POSITIONS (rHi)
        f.write('#hydrogen positions\n')
        for atom in cec_hydrogen_atoms:
            f.write('qm_hydrogen_%s: POSITION ATOM=%s NOPBC\n' % (atom.id, atom.id))
        f.write('\n')    

        # HEAVYATOM (O/N) POSITIONS (rXj)
        f.write('#heavy atom positions\n')
        for atom in cec_heavy_atoms:
            f.write('qm_heavy_%s: POSITION ATOM=%s NOPBC\n' % (atom.id, atom.id))
        f.write('\n')

        # HEAVYATOM (O/N) WEIGHT FACTORS (WXj)
        f.write('#heavy atom weight factors\n')
        for atom in cec_heavy_atoms:
            f.write('qm_heavy_wf_%s: CONSTANT VALUE=%s \n' % (atom.id, weight_factors[atom.index])) # NB the original weight factors dictionary is zero indexed, hence the difference here
        f.write('\n')

        # INTIAL AND TARGET SITE CV ATOMS
        f.write('#initial site CV atoms (protonated O or N)\n')
        for atom in initial_cv_atoms:
            f.write('initial_site_%s: POSITION ATOM=%s NOPBC\n' % (atom.id, atom.id))
        f.write('\n')
        f.write('#target site CV atoms (protonable O or N)\n')
        for atom in target_cv_atoms:
            f.write('target_site_%s: POSITION ATOM=%s NOPBC\n' % (atom.id, atom.id))
        f.write('\n')

        # PAIRWISE VARIABLES
        f.write('\n#pairwise variables\n\n\n')

        # HYDROGEN-HEAVYATOM VECTORS (rHi - rXj)
        f.write('#hydrogen-heavy atom vectors\n')
        for component in ['x', 'y', 'z']:
            for atom in cec_hydrogen_atoms:
                for heavy in cec_heavy_atoms:    
                    f.write('qm_hydrogen_%s_qm_heavy_%s_%s: MATHEVAL ARG=qm_hydrogen_%s.%s,qm_heavy_%s.%s FUNC=x-y PERIODIC=NO \n' % (atom.id, heavy.id, component, atom.id, component, heavy.id, component))
        f.write('\n')

         # HYDROGEN-HEAVYATOM DISTANCES (dHiXj)   
        f.write('#hydrogen-heavy atom distances (dHiXj)\n')
        #for component in ['x', 'y', 'z']:
        for atom in cec_hydrogen_atoms:
            for heavy in cec_heavy_atoms:    
                f.write('qm_hydrogen_%s_qm_heavy_%s_d:   DISTANCE ATOMS=%s,%s \n'          % (atom.id, heavy.id, atom.id, heavy.id))
        f.write('\n')

        # DEBUG - SUMMED UNSCALED CROSS VECTORS
        f.write('#debug cross sum\n')
        for component in ['x', 'y', 'z']:
            f.write('debug_cross_sum_%s: COMBINE ARG=' % component)
            for atom in cec_hydrogen_atoms:
                for heavy in cec_heavy_atoms:
                    if atom == cec_hydrogen_atoms[-1] and heavy == cec_heavy_atoms[-1]:
                        f.write('qm_hydrogen_%s_qm_heavy_%s_%s' % (atom.id, heavy.id, component))
                    else:
                        f.write('qm_hydrogen_%s_qm_heavy_%s_%s,' % (atom.id, heavy.id, component))
            f.write(' COEFFICIENTS=')
            for atom in cec_hydrogen_atoms:
                for heavy in cec_heavy_atoms:
                    if atom == cec_hydrogen_atoms[-1] and heavy == cec_heavy_atoms[-1]:
                        f.write('1.0')
                    else:
                        f.write('1.0,')
            f.write(' PERIODIC=NO\n')
        f.write('\n')

        # SWITCHING FUNCTION PARAMETERS (rsw, dsw)
        f.write('#switching function parameters\n')
        f.write('rsw: CONSTANT VALUE=%s\n' % rsw)
        f.write('dsw: CONSTANT VALUE=%s\n' % dsw)
        f.write('\n')

        # CORRECTION TERM fsw(dHiXj)*(rHi - rXj)
        f.write('#hydrogen-heavy atom vectors scaled by switching function\n')
        for component in ['x', 'y', 'z']:
            for atom in cec_hydrogen_atoms:
                for heavy in cec_heavy_atoms:
                    f.write('qm_hydrogen_%s_qm_heavy_%s_%s_scaled: MATHEVAL ARG=qm_hydrogen_%s_qm_heavy_%s_d,rsw,dsw,qm_hydrogen_%s_qm_heavy_%s_%s FUNC=d*(1.0/(1.0+exp((a-b)/c))) VAR=a,b,c,d PERIODIC=NO\n' % (atom.id, heavy.id, component, atom.id, heavy.id, atom.id, heavy.id, component))
                    #f.write('qm_hydrogen_%s_qm_heavy_%s_%s_scaled: MATHEVAL ARG=qm_hydrogen_%s_qm_heavy_%s_d,rsw,dsw,qm_hydrogen_%s_qm_heavy_%s_%s FUNC=d*(1.0/(1.0+(exp(a-b)/c))) VAR=a,b,c,d PERIODIC=NO\n' % (atom.id, heavy.id, component, atom.id, heavy.id, atom.id, heavy.id, component))
        f.write('\n')


        # SUMMED TERMS 
        f.write('\n#summed terms and CEC\n\n\n')

        # SUMMED HYDROGEN VECTORS (A)
        f.write('#hydrogen vector sum (A)\n')
        for component in ['x', 'y', 'z']:
            f.write('qm_hydrogen_sum_%s: COMBINE ARG=' % component)
            for atom in cec_hydrogen_atoms:
                if atom == cec_hydrogen_atoms[-1]:
                    f.write('qm_hydrogen_%s.%s' % (atom.id, component))
                else:
                    f.write('qm_hydrogen_%s.%s,' % (atom.id, component))
            f.write(' COEFFICIENTS=')
            for atom in cec_hydrogen_atoms:
                if atom == cec_hydrogen_atoms[-1]:
                    f.write('1.0' % atom.id)
                else:
                    f.write('1.0,' % atom.id)
            f.write(' PERIODIC=NO\n')
        f.write('\n')

        # SUMMED WEIGHTED HEAVYATOM VECTORS (B)
        f.write('#heavy-atom weighted vector sum (B)\n')
        for component in ['x', 'y', 'z']:
            f.write('qm_heavy_weighted_sum_%s: COMBINE ARG=' % component)
            for atom in cec_heavy_atoms:
                if atom == cec_heavy_atoms[-1]:
                    f.write('qm_heavy_%s.%s' % (atom.id, component))
                else:
                    f.write('qm_heavy_%s.%s,' % (atom.id, component))
            f.write(' COEFFICIENTS=')
            for atom in cec_heavy_atoms:
                if atom == cec_heavy_atoms[-1]:
                    f.write('%s' % weight_factors[atom.index])
                else:
                    f.write('%s' % weight_factors[atom.index]+',')
            f.write(' PERIODIC=NO\n')
        f.write('\n')

        # SUMMED CORRECTION TERM (C)
        f.write('#scaled pairwise sum (C)\n')
        for component in ['x', 'y', 'z']:
            f.write('qm_pairwise_scaled_vectors_%s: COMBINE ARG=' % component)
            for atom in cec_hydrogen_atoms:
                for heavy in cec_heavy_atoms:
                    if atom == cec_hydrogen_atoms[-1] and heavy == cec_heavy_atoms[-1]:
                        f.write('qm_hydrogen_%s_qm_heavy_%s_%s_scaled' % (atom.id, heavy.id, component))
                    else:
                        f.write('qm_hydrogen_%s_qm_heavy_%s_%s_scaled,' % (atom.id, heavy.id, component))
            f.write(' COEFFICIENTS=')
            for atom in cec_hydrogen_atoms:
                for heavy in cec_heavy_atoms:
                    if atom == cec_hydrogen_atoms[-1] and heavy == cec_heavy_atoms[-1]:
                        f.write('1.0')
                    else:
                        f.write('1.0,')
            f.write(' PERIODIC=NO\n')
        f.write('\n')

        # CENTRE OF EXCESS CHARGE
        f.write('#centre of excess charge (A - B - C)\n')
        for component in ['x', 'y', 'z']:
            f.write('qm_centre_of_excess_charge_%s: MATHEVAL VAR=a,b,c ARG=qm_hydrogen_sum_%s,qm_heavy_weighted_sum_%s,qm_pairwise_scaled_vectors_%s FUNC=(a-b-c) PERIODIC=NO\n' % (component, component, component, component))
        f.write('\n')

        # COLVAR TERMS #####################################################################################
        f.write('\n#colvar terms\n\n')

        # INITAL AND TARGET ATOMS CEC DISTANCES 
        f.write('\n# distances between CEC and initial and target sites\n')
        for atom in initial_cv_atoms:
            for component in ['x', 'y', 'z']:
                f.write('%s-qm_centre_of_excess_charge_%s: COMBINE ARG=initial_site_%s.%s,qm_centre_of_excess_charge_%s COEFFICIENTS=1.0,-1.0 PERIODIC=NO\n' % (atom.id, component, atom.id, component, component))
        for atom in initial_cv_atoms:
            f.write('%s-qm_centre_of_excess_charge_squared: COMBINE ARG=%s-qm_centre_of_excess_charge_x,%s-qm_centre_of_excess_charge_y,%s-qm_centre_of_excess_charge_z POWERS=2.0,2.0,2.0 PERIODIC=NO\n' % (atom.id, atom.id, atom.id, atom.id))
        f.write('\n# these are rsubD(1/2) in the paper\n')
        for atom in initial_cv_atoms:
            f.write('%s-qm_cec_distance: COMBINE ARG=%s-qm_centre_of_excess_charge_squared POWERS=0.5 COEFFICIENTS=1 PERIODIC=NO\n' % (atom.id, atom.id))
        f.write('\n')
        for atom in target_cv_atoms:
            for component in ['x', 'y', 'z']:
                f.write('%s-qm_centre_of_excess_charge_%s: COMBINE ARG=target_site_%s.%s,qm_centre_of_excess_charge_%s COEFFICIENTS=1.0,-1.0 PERIODIC=NO\n' % (atom.id, component, atom.id, component, component))
        for atom in target_cv_atoms:
            f.write('%s-qm_centre_of_excess_charge_squared: COMBINE ARG=%s-qm_centre_of_excess_charge_x,%s-qm_centre_of_excess_charge_y,%s-qm_centre_of_excess_charge_z POWERS=2.0,2.0,2.0 PERIODIC=NO\n' % (atom.id, atom.id, atom.id, atom.id))
        f.write('\n# these are rsubE(1/2) in the paper\n')
        for atom in target_cv_atoms:
            f.write('%s-qm_cec_distance: COMBINE ARG=%s-qm_centre_of_excess_charge_squared POWERS=0.5 COEFFICIENTS=1 PERIODIC=NO\n' % (atom.id, atom.id))
        f.write('\n')

        # SORT THE CEC DISTANCES (RD1/RE1 < RD2/RE2)
        f.write('RD: SORT ARG=%s-qm_cec_distance,%s-qm_cec_distance\n' % (initial_cv_atoms[0].id, initial_cv_atoms[1].id))  # lowest is indexed with 1
        f.write('RE: SORT ARG=%s-qm_cec_distance,%s-qm_cec_distance\n' % (target_cv_atoms[0].id, target_cv_atoms[1].id))
        f.write('\n')

        # WEIGHTS FOR COLVAR VECTOR PROJECTION
        f.write('# weights for initial-target distance\n')
        f.write('sigma: CONSTANT VALUE=%s\n' % (sigma))
        f.write('W1: MATHEVAL ARG=RD.1,RD.2,sigma FUNC=1/((exp((x-y)/z))+1) PERIODIC=NO\n')
        f.write('W2: MATHEVAL ARG=RD.2,RD.1,sigma FUNC=1/((exp((x-y)/z))+1) PERIODIC=NO\n')
        f.write('W3: MATHEVAL ARG=RE.1,RE.2,sigma FUNC=1/((exp((x-y)/z))+1) PERIODIC=NO\n')
        f.write('W4: MATHEVAL ARG=RE.2,RE.1,sigma FUNC=1/((exp((x-y)/z))+1) PERIODIC=NO\n')

        # FINAL COLVAR TERMS (LAW OF COSINES)
        f.write('\n# first two terms for projection colvar\n')
        f.write('RD_final: MATHEVAL ARG=W1,RD.1,W2,RD.2 VAR=w,x,y,z FUNC=(w*x)+(y*z) PERIODIC=NO\n')
        f.write('RE_final: MATHEVAL ARG=W3,RE.1,W4,RE.2 VAR=w,x,y,z FUNC=(w*x)+(y*z) PERIODIC=NO\n')
        f.write('\n')

        # COMPLICATED ROUTINE TO SORT RA1B1, RA2B1, RA1B2, RA2B2 
        # such that A1 and B1 are closer to C than A2 and B2, without referencing atom names after sorting
        # C is the centre of excess charge
        f.write('# calculation of final cross-term for projection colvar\n')

        # sort A indices 
        sf=1000 # scale factor for below

        #RA1B1 = distance 
        #RA2B1 = distance
        #
        f.write('RA1_B1: DISTANCE ATOMS=%s,%s \n' % (initial_cv_atoms[0].id, target_cv_atoms[0].id))
        f.write('RA2_B1: DISTANCE ATOMS=%s,%s \n' % (initial_cv_atoms[1].id, target_cv_atoms[0].id))

        #RA1C = distance # already defined as %s-qm_cec_distance
        #RA2C = distance # already defined as %s-qm_cec_distance
        #
        f.write('RA1_C: MATHEVAL ARG=%s-qm_cec_distance FUNC=x PERIODIC=NO \n' % (initial_cv_atoms[0].id))
        f.write('RA2_C: MATHEVAL ARG=%s-qm_cec_distance FUNC=x PERIODIC=NO \n' % (initial_cv_atoms[1].id))

        #RAN_C=sort(RA1C, RA2C)
        #
        f.write('RAN_C: SORT ARG=RA1_C,RA2_C \n')

        #RA1B1_scaledAC = RA1B1 + 10*RA1C
        #RA2B1_scaledAC = RA2B1 + 10*RA2C
        #RAN_B1=sort(RA1B1_scaledAC, RA2B1_scaledAC)
        #RA1B1 = RAN_B1.1 - 10*RAN_C.1
        #RA2B1 = RAN_B1.2 - 10*RAN_C.2
        #
        f.write('RA1_B1_scaledAC: MATHEVAL ARG=RA1_B1,RA1_C FUNC=x+y*%s PERIODIC=NO \n' % sf)
        f.write('RA2_B1_scaledAC: MATHEVAL ARG=RA2_B1,RA2_C FUNC=x+y*%s PERIODIC=NO \n' % sf)
        f.write('RAN_B1: SORT ARG=RA1_B1_scaledAC,RA2_B1_scaledAC \n')
        f.write('RA1_B1_Asorted: MATHEVAL ARG=RAN_B1.1,RAN_C.1 FUNC=x-y*%s PERIODIC=NO \n' % sf)
        f.write('RA2_B1_Asorted: MATHEVAL ARG=RAN_B1.2,RAN_C.2 FUNC=x-y*%s PERIODIC=NO \n' % sf)

        # then do same for B2 (arbitraty B indices)
        #RA1B2 = distance 
        #RA2B2 = distance
        #
        f.write('RA1_B2: DISTANCE ATOMS=%s,%s \n' % (initial_cv_atoms[0].id, target_cv_atoms[1].id))
        f.write('RA2_B2: DISTANCE ATOMS=%s,%s \n' % (initial_cv_atoms[1].id, target_cv_atoms[1].id))

        #RA1B2_scaledAC = RA1B2 + 10*RA1C
        #RA2B2_scaledAC = RA2B2 + 10*RA2C
        #RAN_B2=sort(RA1B2_scaledAC, RA2B2_scaledAC)
        #RA1B2 = RAN_B2.1 - 10*RAN_C.1
        #RA2B2 = RAN_B2.2 - 10*RAN_C.2
        #
        f.write('RA1_B2_scaledAC: MATHEVAL ARG=RA1_B2,RA1_C FUNC=x+y*%s PERIODIC=NO \n' % sf)
        f.write('RA2_B2_scaledAC: MATHEVAL ARG=RA2_B2,RA2_C FUNC=x+y*%s PERIODIC=NO \n' % sf)
        f.write('RAN_B2: SORT ARG=RA1_B2_scaledAC,RA2_B2_scaledAC \n')
        f.write('RA1_B2_Asorted: MATHEVAL ARG=RAN_B2.1,RAN_C.1 FUNC=x-y*%s PERIODIC=NO \n' % sf)
        f.write('RA2_B2_Asorted: MATHEVAL ARG=RAN_B2.2,RAN_C.2 FUNC=x-y*%s PERIODIC=NO \n' % sf)

        # now A is sorted, sort B indicess

        #RB1C = distance # already defined as %s-qm_cec_distance
        #RB2C = distance # already defined as %s-qm_cec_distance
        #
        f.write('RB1_C: MATHEVAL ARG=%s-qm_cec_distance FUNC=x PERIODIC=NO \n' % (target_cv_atoms[0].id))
        f.write('RB2_C: MATHEVAL ARG=%s-qm_cec_distance FUNC=x PERIODIC=NO \n' % (target_cv_atoms[1].id))

        #RBN_C=sort(RB1C, RB2C)
        #
        f.write('RBN_C: SORT ARG=RB1_C,RB2_C \n')

        #RA1B1_scaledBC = RA1B1 + 10*RB1C
        #RA1B2_scaledBC = RA1B2 + 10*RB2C
        #RA1_BN=sort(RA1B1_scaledBC, RA1B2_scaledBC)
        #RA2B1_scaledBC = RA2B1 + 10*RB1C
        #RA2B2_scaledBC = RA2B2 + 10*RB2C
        #RA2_BN=sort(RA2B1_scaledBC, RA2B2_scaledBC)
        #
        f.write('RA1_B1_scaledBC: MATHEVAL ARG=RA1_B1_Asorted,RB1_C FUNC=x+y*%s PERIODIC=NO \n' % sf)
        f.write('RA1_B2_scaledBC: MATHEVAL ARG=RA1_B2_Asorted,RB2_C FUNC=x+y*%s PERIODIC=NO \n' % sf)
        f.write('RA1_BN: SORT ARG=RA1_B1_scaledBC,RA1_B2_scaledBC \n')
        f.write('RA2_B1_scaledBC: MATHEVAL ARG=RA2_B1_Asorted,RB1_C FUNC=x+y*%s PERIODIC=NO \n' % sf)
        f.write('RA2_B2_scaledBC: MATHEVAL ARG=RA2_B2_Asorted,RB2_C FUNC=x+y*%s PERIODIC=NO \n' % sf)
        f.write('RA2_BN: SORT ARG=RA2_B1_scaledBC,RA2_B2_scaledBC \n')

        #RA1B1 = RA1_BN.1 - 10*RBN_C.1
        #RA1B2 = RA1_BN.2 - 10*RBN_C.2
        #RA2B1 = RA2_BN.1 - 10*RBN_C.1
        #RA2B2 = RA2_BN.2 - 10*RBN_C.2
        #
        #f.write('RA1_B1_Asorted_Bsorted: MATHEVAL ARG=RA1_BN.1,RBN_C.1 FUNC=(x-y*%s) PERIODIC=NO \n' % sf) # RD1E1
        #f.write('RA1_B2_Asorted_Bsorted: MATHEVAL ARG=RA1_BN.2,RBN_C.2 FUNC=(x-y*%s) PERIODIC=NO \n' % sf) # RD1E2
        #f.write('RA2_B1_Asorted_Bsorted: MATHEVAL ARG=RA2_BN.1,RBN_C.1 FUNC=(x-y*%s) PERIODIC=NO \n' % sf) # RD2E1
        #f.write('RA2_B2_Asorted_Bsorted: MATHEVAL ARG=RA2_BN.2,RBN_C.2 FUNC=(x-y*%s) PERIODIC=NO \n' % sf) # RD2E2
        f.write('RD1E1: MATHEVAL ARG=RA1_BN.1,RBN_C.1 FUNC=(x-y*%s) PERIODIC=NO \n' % sf) # RD1E1
        f.write('RD1E2: MATHEVAL ARG=RA1_BN.2,RBN_C.2 FUNC=(x-y*%s) PERIODIC=NO \n' % sf) # RD1E2
        f.write('RD2E1: MATHEVAL ARG=RA2_BN.1,RBN_C.1 FUNC=(x-y*%s) PERIODIC=NO \n' % sf) # RD2E1
        f.write('RD2E2: MATHEVAL ARG=RA2_BN.2,RBN_C.2 FUNC=(x-y*%s) PERIODIC=NO \n' % sf) # RD2E2
                                                                 #          (a*c*i) + (a*d*f) + (b*c*g) + (b*d*h)
        f.write('\n# final cross-term for projection colvar\n')  #   RDE = (w1*w3*RD1E1) + (w1*w4*RD1E2) + (w2*w3*RD2E1) + (w2*w4*RD2E2)
        f.write('RDE_final: MATHEVAL ARG=W1,W2,W3,W4,RD1E1,RD1E2,RD2E1,RD2E2 VAR=a,b,c,d,i,f,g,h FUNC=(a*c*i)+(a*d*f)+(b*c*g)+(b*d*h) PERIODIC=NO\n')   # something is fucking dodgy about this - i think e is getting substituted
        f.write('\n')

        # FINAL CV 
        f.write('# final CV\n')
        f.write('CV: MATHEVAL ARG=RD_final,RDE_final,RE_final VAR=a,b,c FUNC=(((a*a)+(b*b)-(c*c))/(2*(b*b))) PERIODIC=NO\n')

        # FINAL PRINTOUTS     
        f.write('\n# print final CV\n')
        f.write('PRINT STRIDE=1 ARG=CV FILE=CV_PROJECTION\n')
        f.write('\n')

        f.write('\n#biasing potentials\n\n')
        # BIASING POTENTIAL ##########################################################################################

        # fc for biases
        kappa=1000
        #kappa_final=10000

        #upper and lower wall for CV
        f.write('\n# upper wall for CV\n')
        f.write('UPPER_WALLS ARG=CV AT=1.0 KAPPA=1000\n')
        f.write('LOWER_WALLS ARG=CV AT=0.0 KAPPA=1000\n')

        # # moving restraint - SMD biasing restraint - similar to paper
        # f.write('\n# moving restraint\n')
        # f.write('smd: MOVINGRESTRAINT ...\n')
        # f.write('ARG=CV VERSE=B\n')
        # f.write('STEP0=0    AT0=0.00 KAPPA0=%s\n' % kappa)
        # f.write('STEP1=800  AT1=1.00 KAPPA1=%s\n' % kappa)
        # f.write('...\n')
        # f.write('PRINT ARG=smd.bias,smd.work FILE=SMD_BIAS\n')

        # # abmd biasing restraint
        # f.write('\n# abmd restraint\n')
        # f.write('abmd: ABMD ARG=CV TO=1.0 KAPPA=%s\n' % kappa)
        # f.write('PRINT ARG=abmd.bias FILE=ABMD_BIAS\n')

        # metadynamics bias
        sigma_mtd =   0.10 #0.01    # 0.02   # 0.2 Angstroms
        height =      8.4  #16.8    # 4.2    # ~ 2*1 kcal mol-1 = ~8.4 kJ
        pace =        25   #10      # 11     # should be 25 fs - large enough for unbiased DOFs to relax
        biasfactor =  16   # 16    # based on an expected barrier of ~9 kcal mol-1 for ASP deprotonation
        temp =        310  # 310
        #
        grid_min = -5.0
        grid_max =  5.0
        #
        f.write('\n# metadynamics\n')
        f.write('METAD ARG=CV SIGMA=%s HEIGHT=%s PACE=%s FILE=HILLS GRID_MIN=%s GRID_MAX=%s TEMP=%s BIASFACTOR=%s LABEL=metad\n' % (sigma_mtd, height, pace, grid_min, grid_max, temp, biasfactor))
        f.write('\n')

        f.write('\n# print metadynamics\n')
        f.write('PRINT STRIDE=1 ARG=metad.* FILE=METAD_BIAS\n')

        if debug == True:

            # check atoms
            print('\n')
            for atom in cec_atoms:
                print('\tIndex:', atom.index, 'ID:', atom.id, ' resID:', atom.resid, 'resname', atom.resname)
            print('\n')

            # DEBUGGING PRINTOUTS ########################################################################################
            f.write('\n#debugging printouts\n\n')
            f.write('\n# variables\n')
            f.write('PRINT STRIDE=1 ARG=RD_final,RE_final,RDE_final FILE=SUMMED_TERMS_CV\n')
            f.write('PRINT STRIDE=1 ARG=qm_centre_of_excess_charge_x,qm_centre_of_excess_charge_y,qm_centre_of_excess_charge_z FILE=CEC\n')
            f.write('PRINT STRIDE=1 ARG=qm_hydrogen_sum_x,qm_hydrogen_sum_y,qm_hydrogen_sum_z FILE=SUM_HYD \n')
            f.write('PRINT STRIDE=1 ARG=qm_heavy_weighted_sum_x,qm_heavy_weighted_sum_y,qm_heavy_weighted_sum_z FILE=SUM_HEAVY \n')
            f.write('PRINT STRIDE=1 ARG=qm_pairwise_scaled_vectors_x,qm_pairwise_scaled_vectors_y,qm_pairwise_scaled_vectors_z FILE=SUM_CROSS\n')
            #f.write('PRINT STRIDE=1 ARG=W1,RA1_B1_Asorted_Bsorted,W2,RA1_B2_Asorted_Bsorted,W3,RA2_B1_Asorted_Bsorted,W4,RA2_B2_Asorted_Bsorted FILE=weights\n')
            f.write('PRINT STRIDE=1 ARG=W1,W2,W3,W4,RD1E1,RD1E2,RD2E1,RD2E2,RDE_final FILE=WEIGHTS\n')
            # COG/COM of all heavy atoms
            f.write('PRINT ARG=position_qm_com.x,position_qm_com.y,position_qm_com.z FILE=QMCOM\n\n')
            #for atom in initial_cv_atoms:
            #    if atom == initial_cv_atoms[-1]:
            #        f.write('%s-qm_cec_distance FILE=QMCEC_DIST_INITIAL\n' % (atom.id))
            #    else:
            #        f.write('%s-qm_cec_distance,' % (atom.id))
            #f.write('PRINT STRIDE=1 ARG=')
            #for atom in target_cv_atoms:
            #    if atom == target_cv_atoms[-1]:
            #        f.write('%s-qm_cec_distance FILE=QMCEC_DIST_TARGET\n' % (atom.id))
            #    else:
            #        f.write('%s-qm_cec_distance,' % (atom.id))
            #f.write('# coords\n')
            #write a file with all qm_hydrogen_%s_qm_heavy_%s for each atom
            f.write('PRINT FILE=UNWEIGHTED_CROSS_SUM ARG=debug_cross_sum_x,debug_cross_sum_y,debug_cross_sum_z')
            f.write('\n')
            # DEBUG - SUM of all the cross distances (not vectors)
            f.write('#unweighted cross sum distances\n')
            f.write('debug_cross_sum_d: COMBINE ARG=')
            for atom in cec_hydrogen_atoms:
                for heavy in cec_heavy_atoms:
                    if atom == cec_hydrogen_atoms[-1] and heavy == cec_heavy_atoms[-1]:
                        f.write('qm_hydrogen_%s_qm_heavy_%s_d' % (atom.id, heavy.id))
                    else:
                        f.write('qm_hydrogen_%s_qm_heavy_%s_d,' % (atom.id, heavy.id))
            f.write(' COEFFICIENTS=')
            for atom in cec_hydrogen_atoms:
                for heavy in cec_heavy_atoms:
                    if atom == cec_hydrogen_atoms[-1] and heavy == cec_heavy_atoms[-1]:
                        f.write('1.0')
                    else:
                        f.write('1.0,')
            f.write(' PERIODIC=NO\n')
            #write a file with all qm_hydrogen_%s_qm_heavy_%s for each atom
            f.write('PRINT FILE=UNWEIGHTED_CROSS_SUM_D ARG=debug_cross_sum_d')
            f.write('\n')

            # dump hygrogen positions heavy atom positions for debugging
            hydrogen_atoms = []
            heavy_atoms = []
            all_atoms = []
            for atom in cec_hydrogen_atoms:
                hydrogen_atoms.append(atom.id)
            for atom in cec_heavy_atoms:
                heavy_atoms.append(atom.id)
            #for atom in u.select_atoms('all'):
            #    all_atoms.append(atom.id)
            hydrogen_atoms = str(hydrogen_atoms).replace('[','').replace(']','').replace('\'','').replace(' ','')
            heavy_atoms = str(heavy_atoms).replace('[','').replace(']','').replace('\'','').replace(' ','')
            #all_atoms = str(all_atoms).replace('[','').replace(']','').replace('\'','').replace(' ','')
            f.write('DUMPATOMS STRIDE=1 FILE=hydrogens.gro ATOMS=%s\n' % ((hydrogen_atoms)))
            f.write('DUMPATOMS STRIDE=1 FILE=heavyatoms.gro ATOMS=%s\n' % ((heavy_atoms)))
            #f.write('DUMPATOMS STRIDE=1 FILE=qmatoms.gro ATOMS=%s\n' % ((all_atoms)))
