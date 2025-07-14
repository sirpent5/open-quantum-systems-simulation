
# from qiskit.quantum_info import SparsePauliOp
# import numpy as np
from imports import *
def hamiltonian_generation(eps, gamma, F_R,F_L,mu_L,mu_R):
    """
    Generates the Hamiltonian for the system of a single qubit coupled to a reservoir.
    
    Inputs:
        eps: (float) Coupling strength of the qubit to the reservoir.
        gamma: (float) Coupling strength of the qubit to the environment.
        mu: (float) Chemical potential of the reservoir.
        T: (float) Temperature of the reservoir.
    Returns:
        hamiltonian_re: SparsePauliOp representing the real part of the Hamiltonian of the system.
        hamiltonian_im: SparsePauliOp representing the imaginary part of the Hamiltonian of the system.
    """

    # hamiltonian_re = SparsePauliOp(["IZ", "ZI", "XY", "YX", "XY", "YX"], coeffs=[-eps / 2, eps / 2, -(gamma * (1 - 2*F_L)) / 4, -(gamma * (1 - 2*F_L)) / 4, -(gamma * (1 - 2*F_R)) / 4, -(gamma * (1 - 2*F_R)) / 4])
 
    # hamiltonian_im_L = -1 * SparsePauliOp(["XX", "YY", "II", "IZ", "ZI"], coeffs=[gamma / 4, -gamma / 4, -gamma / 2, (gamma * (1 - 2*F_L)) / 4, (gamma * (1 - 2*F_L)) / 4])
    # hamiltonian_im_R = -1 * SparsePauliOp(["XX", "YY", "II", "IZ", "ZI"], coeffs=[gamma / 4, -gamma / 4, -gamma / 2, (gamma * (1 - 2*F_R)) / 4, (gamma * (1 - 2*F_R)) / 4])
    # hamiltonian_im = hamiltonian_im_L + hamiltonian_im_R
    

    common_coeff = -gamma / 2 * (1 - F_L - F_R)

    hamiltonian_re = SparsePauliOp(
        ["IZ", "ZI", "XY", "YX"],
        coeffs=[-eps / 2, eps / 2, common_coeff, common_coeff]
    )

    hamiltonian_im = SparsePauliOp(
        ["XX", "YY", "II", "IZ", "ZI"],
        coeffs=[-gamma / 2, gamma / 2, gamma, common_coeff, common_coeff]
    )
    return hamiltonian_re, hamiltonian_im

def hamiltonian_generation_simple():
    """
    Generates a simple Hamiltonian for a single qubit system.

    Returns:
        hamiltonian_re: SparsePauliOp representing the Hamiltonian.
    """
    return SparsePauliOp(["IX", "XI"], coeffs=[1, -1])  # Example coefficients

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

    trace_list = [1.0]
    num_op = 0.5 * SparsePauliOp("III") - 0.5 * SparsePauliOp("IIZ")
    
    initial_exp_val = init_state.expectation_value(num_op).real
    num_op_list = [initial_exp_val]
    
    print(f"Initial expectation value of number operator: {initial_exp_val:.4f}")

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

        true_trace = np.trace(rho_matrix)

        exp_val = np.trace(rho_matrix @ np.array([[0, 0], [0, 1]])) / true_trace
        
        num_op_list.append(exp_val.real)
        trace_list.append(1.0)


    return num_op_list, trace_list

def output_vqte_results(vqte_results, time, nt, eps, mu_L,mu_R,T_L, T_R):

    plt.figure(figsize=(10, 6))
    time_axis = np.linspace(0, time, nt+1)
    mu_effective = (mu_L + mu_R) / 2
    T_effective = (T_L + T_R) / 2
    time_axis = np.linspace(0, time, nt+1)
    plt.plot(time_axis, [1 / (1 + np.exp((eps - mu_effective) / T_effective))] * (nt+1), label='Steady State Expectation Value', linestyle='solid')
    plt.plot(time_axis/2, vqte_results,marker='', linestyle='dashed', label='VQTE Result', color='blue')

    plt.title("VQTE Results of two reserviors coupled to a single qubit")
    plt.xlabel("Time (t)")
    plt.ylabel("⟨n⟩ (Expectation Value)")
    plt.grid(True)
    plt.legend()
    
    plt.show()