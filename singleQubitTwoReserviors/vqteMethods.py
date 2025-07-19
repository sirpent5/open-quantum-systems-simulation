from imports import *

def hamiltonian_generation(eps, gamma_L, gamma_R, F_R,F_L):
#     """
#     Generates the Hamiltonian for the system of a single qubit coupled to a reservoir.
    
#     Inputs:
#         eps: (float) Coupling strength of the qubit to the reservoir.
#         gamma: (float) Coupling strength of the qubit to the environment.
#         mu: (float) Chemical potential of the reservoir.
#         T: (float) Temperature of the reservoir.
#     Returns:
#         hamiltonian_re: SparsePauliOp representing the real part of the Hamiltonian of the system.
#         hamiltonian_im: SparsePauliOp representing the imaginary part of the Hamiltonian of the system.
#     """

  
    # hamiltonian_re = SparsePauliOp( ["IZ", "ZI", "XY", "YX"],
    #     coeffs=[(-eps/2),(eps/2), ((-1/4)*((gamma_L*(1-2*F_L))+(gamma_R*(1-2*F_R)))),((-1/4)*((gamma_L*(1-2*F_L))+(gamma_R*(1-2*F_R))))]
    # )


    # hamiltonian_im = SparsePauliOp( ["XX", "YY", "II", "IZ", "ZI"],
    #     coeffs=[-(gamma_L+gamma_R)/4, (gamma_L+gamma_R)/4, (gamma_L+gamma_R)/2,
    #              ((-1/4)*((gamma_L*(1-2*F_L))+(gamma_R*(1-2*F_R)))), ((-1/4)*((gamma_L*(1-2*F_L))+(gamma_R*(1-2*F_R)))) ] )

    
    # return hamiltonian_re, hamiltonian_im
    gamma =1
    hamiltonian_re_coeffs = {
        # System energy influenced by both reservoirs
        "IIZ": -eps,
        # Energy of Reservoir 1
        "IZI": eps/ 2,
        # Energy of Reservoir 2
        "ZII": eps / 2,
        # System <-> Reservoir 1 coupling
        "XYI": -(gamma * F_L) / 4,
        "YXI": -(gamma * F_L) / 4,
        # System <-> Reservoir 2 coupling
        "XIY": -(gamma * F_R) / 4,
        "YIX": -(gamma * F_R) / 4, }

    hamiltonian_im_coeffs = {
        # Global term from both reservoirs
        "III": -gamma,
        # System <-> Res 1 dissipative interaction
        "XXI": gamma / 4,
        "YYI": -gamma / 4,
        # System <-> Res 2 dissipative interaction
        "XIX": gamma / 4,
        "YIY": -gamma / 4,
        # Bias terms for System-Res1
        "IIZ": (gamma * F_L) / 4 + (gamma * F_R) / 4, # System bias from both
        "IZI": (gamma * F_L) / 4,                     # Res 1 bias
        "ZII": (gamma * F_R) / 4,                     # Res 2 bias
    }
    hamiltonian_re = SparsePauliOp(list(hamiltonian_re_coeffs.keys()), coeffs=list(hamiltonian_re_coeffs.values()))
    hamiltonian_im = -1 * SparsePauliOp(list(hamiltonian_im_coeffs.keys()), coeffs=list(hamiltonian_im_coeffs.values()))

    return hamiltonian_re, hamiltonian_im


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
        exp_val = np.trace(rho_matrix @ np.array([[0, 0], [0, 1]])) / true_trace
        
        num_op_list.append(exp_val.real)



    return num_op_list