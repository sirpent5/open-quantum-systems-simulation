
# from qiskit.quantum_info import SparsePauliOp
# import numpy as np
from imports import *

def hamiltonian_generation(N, eps, gamma, F_L, F_R,j):
    """
    Generates the real and imaginary parts of an effective Hamiltonian for an N-qubit chain.

    The real part corresponds to on-site energy terms. The imaginary part corresponds 
    to dissipation from reservoirs coupled to the first and last qubits.

    Args:
        N (int): The number of qubits in the chain.
        eps (float): The on-site energy for each qubit (coefficient for Z terms).
        gamma (float): The coupling strength of the reservoirs.
        F_L (float): The Fermi-Dirac distribution factor for the left reservoir.
        F_R (float): The Fermi-Dirac distribution factor for the right reservoir.

    Returns:
        (SparsePauliOp, SparsePauliOp): A tuple containing the real (Hermitian) and
                                         imaginary (dissipative) parts of the Hamiltonian.
   """
    n = 2*N
    hopping_coeff = j * (1 - F_L - F_R)
    # 1. Build the real part (H_re): The system Hamiltonian
    # Coefficients for reservoir-coupled Z terms
    coeff_left = -gamma / 2 * (1 - F_L)   # Z on first qubit
    coeff_right = -gamma / 2 * (1 - F_R)   # Z on last qubit

    # Initialize empty lists for Pauli strings and coefficients
    pauli_list_re = []
    coeffs_re = []
    pauli_list_im = []
    coeffs_im = []

    # Single-qubit Z terms 
    for i in range(n):
        sign = (-1) ** i
        z_term = ["I"] * n
        z_term[i] = "Z"
        pauli_list_re.append("".join(z_term))
        coeffs_re.append(sign * eps / 2)

    # Nearest-neighbor XY and YX terms (scaled by j)
    for i in range(n - 1):
        # XY term
        xy_term = ["I"] * n
        xy_term[i] = "X"
        xy_term[i + 1] = "Y"
        pauli_list_re.append("".join(xy_term))
        coeffs_re.append(hopping_coeff)

        # YX term
        yx_term = ["I"] * n
        yx_term[i] = "Y"
        yx_term[i + 1] = "X"
        pauli_list_re.append("".join(yx_term))
        coeffs_re.append(hopping_coeff)

    # Hamiltonian_im terms 
    # Nearest-neighbor XX and YY terms
    for i in range(n - 1):
        # XX term
        xx_term = ["I"] * n
        xx_term[i] = "X"
        xx_term[i + 1] = "X"
        pauli_list_im.append("".join(xx_term))
        coeffs_im.append(-gamma / 2)

        # YY term
        yy_term = ["I"] * n
        yy_term[i] = "Y"
        yy_term[i + 1] = "Y"
        pauli_list_im.append("".join(yy_term))
        coeffs_im.append(gamma / 2)

    # Global II...I term (gamma)
    pauli_list_im.append("I" * n)
    coeffs_im.append(gamma)

    # Left reservoir (Z on first qubit)
    z_left = ["I"] * n
    z_left[0] = "Z"
    pauli_list_im.append("".join(z_left))
    coeffs_im.append(coeff_left)

    # Right reservoir (Z on last qubit)
    z_right = ["I"] * n
    z_right[-1] = "Z"
    pauli_list_im.append("".join(z_right))
    coeffs_im.append(coeff_right)

    # Construct the Hamiltonians
    hamiltonian_re = SparsePauliOp(pauli_list_re, coeffs_re)
    hamiltonian_im = SparsePauliOp(pauli_list_im, coeffs_im)
    return hamiltonian_re, hamiltonian_im


def create_number_operator(N, qubit_index):
    """
    Creates the number operator for a specific qubit in an N-qubit system.

    The number operator is defined as n = (I - Z) / 2, which acts as a
    projector onto the |1⟩ state for the specified qubit.

    Args:
        N (int): The total number of qubits in the system.
        qubit_index (int): The index of the qubit for the operator (0 to N-1).

    Returns:
        SparsePauliOp: The number operator for the specified qubit.
    """
    if not 0 <= qubit_index < N:
        raise ValueError("qubit_index must be between 0 and N-1.")

    # Create the two necessary Pauli strings
    identity_str = 'I' * N
    z_pauli_list = ['I'] * N
    z_pauli_list[qubit_index] = 'Z'
    z_str = "".join(z_pauli_list)

    # The number operator is n = 0.5 * I - 0.5 * Z
    number_op = SparsePauliOp([identity_str, z_str], coeffs=[0.5, -0.5])
    
    return number_op
def hamiltonian_generation_simple():
    """
    Generates a simple Hamiltonian for a single qubit system.

    Returns:
        hamiltonian_re: SparsePauliOp representing the Hamiltonian.
    """
    return SparsePauliOp(["IX", "XI"], coeffs=[1, -1])  # Example coefficients

def statevector_to_densitymatrix(state_vector):
    """
    Converts a Statevector to a density matrix.

    Inputs:
        v: (numpy.ndarray) The state vector to be converted.
    
    Returns:
        (numpy.ndarray) The corresponding density matrix.
    """
    
    return np.outer(state_vector, np.conj(state_vector))

def perform_vqte(ham_real, ham_imag, init_state,dt, nt, ansatz, init_param_values,N):
    real_var_principle = RealMcLachlanPrinciple(qgt=ReverseQGT(), gradient=ReverseEstimatorGradient(derivative_type=DerivativeType.IMAG))
    imag_var_principle = ImaginaryMcLachlanPrinciple(qgt=ReverseQGT(), gradient=ReverseEstimatorGradient())


# Initialize lists to store results
#second is the is expectation value of the number operator
   
    number_operators = [create_number_operator(N, i) for i in range(N)]
    # Perform time evolution
    results_history = [[] for _ in range(N)]

    print("Initial expectation values:")
    for i, op in enumerate(number_operators):
        initial_exp_val = init_state.expectation_value(op).real
        results_history[i].append(initial_exp_val)
        print(f"  Qubit {i}: {initial_exp_val:.4f}")

    time_points = [0.0]
    plt.ion()
    fig = plt.subplots(figsize=(10, 6))
    #update_live_plot(results_history, time_points, N)

    plot_interval = 100



    # --- Time Evolution Loop ---
    for t in range(nt):
        print(f"Step {t+1} out of {nt}")
        # Real and Imaginary evolution steps
        evolution_problem_re = TimeEvolutionProblem(ham_real, dt)
        var_qrte = VarQRTE(ansatz, init_param_values, real_var_principle, num_timesteps=1)
        evolution_result_re = var_qrte.evolve(evolution_problem_re)
        init_param_values = evolution_result_re.parameter_values[-1]

        evolution_problem_im = TimeEvolutionProblem(ham_imag, dt)
        var_qite = VarQITE(ansatz, init_param_values, imag_var_principle, num_timesteps=1)
        evolution_result_im = var_qite.evolve(evolution_problem_im)
        init_param_values = evolution_result_im.parameter_values[-1]

        # --- Measurement ---
        current_psi = Statevector(ansatz.assign_parameters(init_param_values))
        norm = np.linalg.norm(current_psi.data)
        normalized_psi = current_psi / norm if norm != 0 else current_psi

        # Calculate expectation value for each qubit's number operator
        for i, op in enumerate(number_operators):
            exp_val = normalized_psi.expectation_value(op).real
            results_history[i].append(exp_val)


        time_points.append((t + 1) * dt)
     # MODIFICATION: Call the live plot function
        # plt.clf()
        #if (t + 1) % plot_interval == 0 or t == nt - 1:
            #update_live_plot(results_history, time_points, N)

    # MODIFICATION: Finalize plotting
    plt.ioff()
    plt.title("Final VQTE Results")
    plt.show()


    return results_history



def update_live_plot(expectation_value_history, time_points, N):
    """
    Clears the current plot and redraws it with the updated data.
    """
    plt.clf()  # Clear the current figure
    
    for site_idx in range(N):
        plt.plot(time_points, expectation_value_history[site_idx], label=f'Qubit {site_idx}', marker='o', markersize=3, linestyle='-')

    plt.title("Live VQTE: Qubit Occupation")
    plt.xlabel("Time (t)")
    plt.ylabel("⟨n⟩ (Expectation Value)")
    plt.grid(True)
    plt.legend()
    plt.ylim(-0.05, 1.05) # Keep axes stable
    
    plt.draw()
    plt.pause(0.01)