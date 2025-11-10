from imports import *

def hamiltonian_generation(n_sites, eps, gamma_L, gamma_R, F_L, F_R, J):
    """
    Returns H_re and H_im for a fermionic chain (1 qubit per site, no spin).
    
    Args:
        n_sites: Number of physical sites (equals number of qubits)
        eps: List of on-site energies [eps_1, ..., eps_n]
        gamma_L, gamma_R: Left/right reservoir couplings
        F_L, F_R: Reservoir occupation factors (0 ≤ F ≤ 1)
        t: Hopping amplitude
        
    Returns:
        H_re, H_im: Real and imaginary parts as SparsePauliOp
    """
    # Initialize Pauli lists and coefficients
    pauli_re, coeffs_re = [], []
    pauli_im, coeffs_im = [], []

    # ===== Real Part (H_re) =====

    N = 2 * n_sites
    n_sites = len(eps)
    N = 2 * n_sites
        

    # Onsite energies

    eps_index = 0
    for i in range(n_sites):
  
        z_str = ['I']* N
        z_str[i] = 'Z'
        pauli_re.append(''.join(z_str))
        coeffs_re.append(eps[eps_index]/2)

        z_str = ['I']* N
        z_str[i+n_sites] = 'Z'
        pauli_re.append(''.join(z_str))
        coeffs_re.append(-eps[eps_index]/2)

        eps_index += 1

    # Dissipation for left and right sides

    xy_str = ['I']* N
    xy_str[0] = 'X'
    xy_str[n_sites] = 'Y'
    pauli_re.append(''.join(xy_str))
    coeffs_re.append(-0.25*gamma_L*(1-2*F_L))

    yx_str = ['I']* N
    yx_str[0] = 'Y'
    yx_str[n_sites] = 'X'
    pauli_re.append(''.join(yx_str))
    coeffs_re.append(-0.25*gamma_L*(1-2*F_L))

    
    yx_str = ['I']* N
    yx_str[n_sites - 1] = 'Y'
    yx_str[N-1] = 'X'
    pauli_re.append(''.join(yx_str))
    coeffs_re.append(-0.25*gamma_R*(1-2*F_R))

    xy_str = ['I']* N
    xy_str[n_sites - 1] = 'X'
    xy_str[N-1] = 'Y'
    pauli_re.append(''.join(xy_str))
    coeffs_re.append(-0.25*gamma_R*(1-2*F_R))

    # Onsite energies

    eps_index = 0
    for i in range(n_sites):
  
        z_str = ['I']* N
        z_str[i] = 'Z'
        pauli_re.append(''.join(z_str))
        coeffs_re.append(eps[eps_index]/2)

        z_str = ['I']* N
        z_str[i+n_sites] = 'Z'
        pauli_re.append(''.join(z_str))
        coeffs_re.append(-eps[eps_index]/2)

        eps_index += 1

    # Dissipation for left and right sides

    xy_str = ['I']* N
    xy_str[0] = 'X'
    xy_str[n_sites] = 'Y'
    pauli_re.append(''.join(xy_str))
    coeffs_re.append(-0.25*gamma_L*(1-2*F_L))

    yx_str = ['I']* N
    yx_str[0] = 'Y'
    yx_str[n_sites] = 'X'
    pauli_re.append(''.join(yx_str))
    coeffs_re.append(-0.25*gamma_L*(1-2*F_L))

    
    yx_str = ['I']* N
    yx_str[n_sites - 1] = 'Y'
    yx_str[N-1] = 'X'
    pauli_re.append(''.join(yx_str))
    coeffs_re.append(-0.25*gamma_R*(1-2*F_R))

    xy_str = ['I']* N
    xy_str[n_sites - 1] = 'X'
    xy_str[N-1] = 'Y'
    pauli_re.append(''.join(xy_str))
    coeffs_re.append(-0.25*gamma_R*(1-2*F_R))


    ## Hopping terms

    for i in range(n_sites - 1):
        
        xx_str = ['I']* N
        xx_str[i], xx_str[i+1] = 'X', 'X'
        pauli_re.append(''.join(xx_str))
        coeffs_re.append(-J)

        xx_str = ['I']* N
        xx_str[n_sites+i], xx_str[n_sites+i+1] = 'X', 'X'
        pauli_re.append(''.join(xx_str))
        coeffs_re.append(J)

        yy_str = ['I']* N
        yy_str[i], yy_str[i+1] = 'Y', 'Y'
        coeffs_re.append(-J)
        pauli_re.append(''.join(yy_str))

        yy_str = ['I']* N
        yy_str[n_sites+i], yy_str[n_sites+i+1] = 'Y', 'Y'
        pauli_re.append(''.join(yy_str))
        coeffs_re.append(J)




# ===== Imaginary Part (H_im) =====

    xx_str = ['I']* N
    xx_str[0] = 'X'
    xx_str[n_sites] = 'X'
    pauli_im.append(''.join(xx_str))
    coeffs_im.append(-gamma_L/4)

    ## Hopping terms

    for i in range(n_sites - 1):
        
        xx_str = ['I']* N
        xx_str[i], xx_str[i+1] = 'X', 'X'
        pauli_re.append(''.join(xx_str))
        coeffs_re.append(-J)

        xx_str = ['I']* N
        xx_str[n_sites+i], xx_str[n_sites+i+1] = 'X', 'X'
        pauli_re.append(''.join(xx_str))
        coeffs_re.append(J)

        yy_str = ['I']* N
        yy_str[i], yy_str[i+1] = 'Y', 'Y'
        coeffs_re.append(-J)
        pauli_re.append(''.join(yy_str))

        yy_str = ['I']* N
        yy_str[n_sites+i], yy_str[n_sites+i+1] = 'Y', 'Y'
        pauli_re.append(''.join(yy_str))
        coeffs_re.append(J)




# ===== Imaginary Part (H_im) =====

    xx_str = ['I']* N
    xx_str[0] = 'X'
    xx_str[n_sites] = 'X'
    pauli_im.append(''.join(xx_str))
    coeffs_im.append(-gamma_L/4)

    yy_str = ['I']* N
    yy_str[0] = 'Y'
    yy_str[n_sites]= 'Y'
    pauli_im.append(''.join(yy_str))
    coeffs_im.append(0.25*gamma_L)

    ## I terms
    I_str = ['I'] * N
    pauli_im.append(''.join(I_str))
    coeffs_im.append(gamma_L/2)

    ## ZZ terms Left
    z_str = ['I']* N
    z_str[0] = 'Z'
    pauli_im.append(''.join(z_str))
    coeffs_im.append(-0.25*gamma_L*(1-2*F_L))

    z_str = ['I']* N
    z_str[n_sites] = 'Z'
    pauli_im.append(''.join(z_str))
    coeffs_im.append(-0.25*gamma_L*(1-2*F_L))

    ## Right reservior imaginary terms

    ##XX term right
    xx_str = ['I']* N
    xx_str[n_sites-1] = 'X'
    xx_str[N-1] = 'X'
    pauli_im.append(''.join(xx_str))
    coeffs_im.append(-gamma_R/4)

    yy_str = ['I']* N
    yy_str[n_sites-1] = 'Y'
    yy_str[N-1]= 'Y'
    pauli_im.append(''.join(yy_str))
    coeffs_im.append(0.25*gamma_R)

    ##II term Right
    I_str = ['I'] * N
    pauli_im.append(''.join(I_str))
    coeffs_im.append(gamma_R/2)

    ## ZZ terms Right
    z_str = ['I']* N
    z_str[n_sites-1] = 'Z'
    pauli_im.append(''.join(z_str))
    coeffs_im.append(-0.25*gamma_R*(1-2*F_R))

    z_str = ['I']* N
    z_str[N-1] = 'Z'
    pauli_im.append(''.join(z_str))
    coeffs_im.append(-0.25*gamma_R*(1-2*F_R))

    ##YY Terms Right



    yy_str = ['I']* N
    yy_str[n_sites-1] = 'Y'
    yy_str[N-1]= 'Y'
    pauli_im.append(''.join(yy_str))
    coeffs_im.append(0.25*gamma_R)

    ##II term Right
    I_str = ['I'] * N
    pauli_im.append(''.join(I_str))
    coeffs_im.append(gamma_R/2)

     ## ZZ terms Right
    z_str = ['I']* N
    z_str[n_sites-1] = 'Z'
    pauli_im.append(''.join(z_str))
    coeffs_im.append(-0.25*gamma_R*(1-2*F_R))

    z_str = ['I']* N
    z_str[N-1] = 'Z'
    pauli_im.append(''.join(z_str))
    coeffs_im.append(-0.25*gamma_R*(1-2*F_R))
    
    return SparsePauliOp(pauli_re, coeffs=np.array(coeffs_re)), SparsePauliOp(pauli_im, coeffs=np.array(coeffs_im))


def statevector_to_densitymatrix(v):
    """
    Converts a Statevector to a density matrix.

    Inputs:
        v: (numpy.ndarray) The state vector to be converted.
    
    Returns:
        (numpy.ndarray) The corresponding density matrix.
    """
    m = int(np.sqrt(len(v)))
    return np.reshape(v, (m, m), order='F')

def perform_vqte(ham_real, ham_imag, init_state, dt, nt, ansatz, init_param_values):
    """
    Performs the VQTE simulation with corrected calculations.
    """
    # Define the variational principles
    real_var_principle = RealMcLachlanPrinciple(qgt=ReverseQGT(), gradient=ReverseEstimatorGradient(derivative_type=DerivativeType.IMAG))
    imag_var_principle = ImaginaryMcLachlanPrinciple(qgt=ReverseQGT(), gradient=ReverseEstimatorGradient())

    # Construct number operator
    num_op = 0.5 * SparsePauliOp("III") - 0.5 * SparsePauliOp("IIZ")
    
    # Get initial expectation value
    initial_exp_val = init_state.expectation_value(num_op).real
    num_op_list = [initial_exp_val]
    vqte_fidelity = [statevector_to_densitymatrix(init_state.data)]

    # --- Perform time evolution ---
    for t in range(nt):
      
        print("Step", t , "out of", nt)
        # Real evolution
        evolution_problem_re = TimeEvolutionProblem(ham_real, dt )
        var_qrte = VarQRTE(ansatz, init_param_values, real_var_principle, num_timesteps=1)
        evolution_result_re = var_qrte.evolve(evolution_problem_re)
        init_param_values = evolution_result_re.parameter_values[-1]
        
        norm_squared = 1.0 
        
        psi_after_re = Statevector(ansatz.assign_parameters(init_param_values))
        exp_val_H_imag = psi_after_re.expectation_value(ham_imag).real
        norm_squared *= (1 + exp_val_H_imag * dt)

        # Imaginary evolution
        evolution_problem_im = TimeEvolutionProblem(ham_imag, dt )
        var_qite = VarQITE(ansatz, init_param_values, imag_var_principle, num_timesteps=1)
        evolution_result_im = var_qite.evolve(evolution_problem_im)
        init_param_values = evolution_result_im.parameter_values[-1]

         # Normalized so the trace is always 1
        final_psi_normalized = Statevector(ansatz.assign_parameters(init_param_values))
        
        # Create the physically correct, unnormalized density matrix by scaling with the tracked norm
        rho_unnormalized_vec = np.sqrt(norm_squared) * final_psi_normalized.data
        rho_matrix = statevector_to_densitymatrix(rho_unnormalized_vec)

        # Extract expectation values
        true_trace = np.trace(rho_matrix)
       
        
        rho_matrix /= true_trace

        exp_val = np.trace(rho_matrix @ np.array([[0, 0], [0, 1]]))
        
        
        num_op_list.append(exp_val.real)
        vqte_fidelity.append(rho_matrix)


    return num_op_list