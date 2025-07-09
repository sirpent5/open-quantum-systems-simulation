
# from qiskit.quantum_info import SparsePauliOp
# import numpy as np
from imports import *
# def hamiltonian_generation(eps, gamma, F_R,F_L):
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
#     F = -((1- 2*F_L)+(1- 2*F_R))

#     hamiltonian_re = SparsePauliOp(["IZ", "ZI", "XY", "YX"], coeffs=[-eps / 2, eps / 2, -(gamma * (1 - 2*F)) / 4, -(gamma * (1 - 2*F)) / 4])
#     hamiltonian_im = -1 * SparsePauliOp(["XX", "YY", "II", "IZ", "ZI"], coeffs=[gamma / 4, -gamma / 4, -gamma / 2, (gamma * (1 - 2*F)) / 4, (gamma * (1 - 2*F)) / 4])
    
#     return hamiltonian_re, hamiltonian_im
def hamiltonian_generation(N, eps, gamma, F_L, F_R):
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

    re_paulis = []
    re_coeffs = []
    im_paulis = []
    im_coeffs = []

    # 1. Build the real part (H_re): The system Hamiltonian
    # On-site Z terms (energy) for each qubit
    for i in range(N):
        pauli_str = ['I'] * N
        pauli_str[i] = 'Z'
        re_paulis.append("".join(pauli_str))
        re_coeffs.append(eps)

    # NOTE: Nearest-neighbor XX and YY interaction terms have been removed.

    # 2. Build the imaginary part (H_im): The dissipative terms
    # Identity term from both reservoirs
    im_paulis.append('I' * N)
    im_coeffs.append(gamma / 2.0)

    # Left reservoir term acting on the first qubit (qubit 0)
    pauli_zl = ['I'] * N
    pauli_zl[0] = 'Z'
    im_paulis.append("".join(pauli_zl))
    im_coeffs.append(-gamma / 4.0 * (1 - 2 * F_L))

    # Right reservoir term acting on the last qubit (qubit N-1)
    pauli_zr = ['I'] * N
    pauli_zr[N-1] = 'Z'
    im_paulis.append("".join(pauli_zr))
    im_coeffs.append(-gamma / 4.0 * (1 - 2 * F_R))

    # Create the SparsePauliOp objects
    # Qiskit automatically sums coefficients for identical Pauli strings (e.g., if N=1)
    hamiltonian_re = SparsePauliOp(re_paulis, coeffs=re_coeffs)
    hamiltonian_im = SparsePauliOp(im_paulis, coeffs=im_coeffs)

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