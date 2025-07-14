
# from qiskit.quantum_info import SparsePauliOp
# import numpy as np
from imports import *



def hamiltonian_generation(N, eps, gamma, F_L, F_R, hopping_coeff):
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



    # Initialize empty lists for Pauli strings and coefficients
    pauli_list_re = []
    coeffs_re = []
    pauli_list_im = []
    coeffs_im = []

    # Single-qubit Z terms 
    for i in range(n):
        z_term = ["I"] * n
        z_term[i] = "Z"
        pauli_list_re.append("".join(z_term))
        coeffs_re.append(eps/2)

    for i in range(n - 1):
        # XX term
        xx_term = ["I"] * n
        xx_term[i] = "X"
        xx_term[i + 1] = "X"  
        pauli_list_re.append("".join(xx_term))
        coeffs_re.append(hopping_coeff)

        # YY term
        yy_term = ["I"] * n
        yy_term[i] = "Y"
        yy_term[i + 1] = "Y"  # Should be 'Y'
        pauli_list_re.append("".join(yy_term))
        coeffs_re.append(hopping_coeff)

    # Hamiltonian_im terms 

   
    # pauli_list_im.append("I" * n)
    # coeffs_im.append(gamma / 2)

    x_left = ["I"] * n
    x_left[0] = "X"
    pauli_list_im.append("".join(x_left))
    coeffs_im.append(-gamma * (1 - F_L) / 4)  # X part

    y_left = ["I"] * n
    y_left[0] = "Y"
    pauli_list_im.append("".join(y_left))
    coeffs_im.append(-gamma * (1 - F_L) / 4)

    x_left_gain = ["I"] * n
    x_left_gain[0] = "X"
    pauli_list_im.append("".join(x_left_gain))
    coeffs_im.append(-gamma * (1 - F_L))

    y_left_gain = ["I"] * n
    y_left_gain[0] = "Y"
    pauli_list_im.append("".join(y_left_gain))
    coeffs_im.append(gamma * F_L/4)  # Note sign flip for -iY


    x_right = ["I"] * n
    x_right[-1] = "X"
    pauli_list_im.append("".join(x_right))
    coeffs_im.append(-gamma * (1 - F_R) / 4)

    y_right = ["I"] * n
    y_right[-1] = "Y"
    pauli_list_im.append("".join(y_right))
    coeffs_im.append(-gamma * (1 - F_R) / 4)

    x_right_gain = ["I"] * n
    x_right_gain[-1] = "X"
    pauli_list_im.append("".join(x_right_gain))
    coeffs_im.append(-gamma * F_R / 4)

    y_right_gain = ["I"] * n
    y_right_gain[-1] = "Y"
    pauli_list_im.append("".join(y_right_gain))
    coeffs_im.append(gamma * F_R / 4)

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

    number_operators = [create_number_operator(N, i) for i in range(N)]
    # Perform time evolution
    results_history = [[] for _ in range(N)]

    print("Initial expectation values:")
    for i, op in enumerate(number_operators):
        initial_exp_val = init_state.expectation_value(op).real
        results_history[i].append(initial_exp_val)
        print(f"  Qubit {i}: {initial_exp_val:.4f}")


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

    return results_history



    
