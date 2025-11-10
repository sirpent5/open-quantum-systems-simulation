
# from qiskit.quantum_info import SparsePauliOp
# import numpy as np
from imports import *

import numpy as np
from qiskit.quantum_info import SparsePauliOp

from qiskit.quantum_info import SparsePauliOp
import numpy as np


        
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

# def hamiltonian_generation(N, eps, gamma_L, gamma_R, F_L, F_R, J):
#     """
#     Generates the real and imaginary parts of an effective Hamiltonian for an N-qubit chain.

#     The real part corresponds to on-site energy terms. The imaginary part corresponds 
#     to dissipation from reservoirs coupled to the first and last qubits.

#     Args:
#         N (int): The number of qubits in the chain.
#         eps (float): The on-site energy for each qubit (coefficient for Z terms).
#         gamma (float): The coupling strength of the reservoirs.
#         F_L (float): The Fermi-Dirac distribution factor for the left reservoir.
#         F_R (float): The Fermi-Dirac distribution factor for the right reservoir.

#     Returns:
#         (SparsePauliOp, SparsePauliOp): A tuple containing the real (Hermitian) and
#                                          imaginary (dissipative) parts of the Hamiltonian.
#    """
  

#     # Initialize empty lists for Pauli strings and coefficients
#     pauli_list_re = []
#     coeffs_re = []
#     pauli_list_im = []
#     coeffs_im = []

#     for site in range(N):
#             hamiltonian_on_site = SparsePauliOp( ["IZ", "ZI", "XY", "YX"],
#             coeffs=[(-eps/2),(eps/2), ((-1/4)*((gamma_L*(1-2*F_L))+(gamma_R*(1-2*F_R)))),((-1/4)*((gamma_L*(1-2*F_L))+(gamma_R*(1-2*F_R))))]
#         )


#     hamiltonian_im = SparsePauliOp( ["XX", "YY", "II", "IZ", "ZI"],
#         coeffs=[-(gamma_L+gamma_R)/4, (gamma_L+gamma_R)/4, (gamma_L+gamma_R)/2,
#                  ((-1/4)*((gamma_L*(1-2*F_L))+(gamma_R*(1-2*F_R)))), ((-1/4)*((gamma_L*(1-2*F_L))+(gamma_R*(1-2*F_R)))) ] )





#     # Construct the Hamiltonians
#     hamiltonian_re = SparsePauliOp(pauli_list_re, coeffs_re)
#     hamiltonian_im = SparsePauliOp(pauli_list_im, coeffs_im)
#     return hamiltonian_re, hamiltonian_im


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

    print("Oh my gosh it is starting.")


    for i, op in enumerate(number_operators):
        op_matrix = op.to_matrix()
        density_matrix = statevector_to_densitymatrix(init_state.data)
  
        initial_exp_val = np.trace(density_matrix @ op_matrix).real
    
        results_history[i].append(initial_exp_val)
    print("Ham_real terms and coefficients:")
    for op, coeff in zip(ham_real.paulis, ham_real.coeffs):
        print(f"{op}: {coeff}")
    print("Oh my gosh it is starting.")


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

    for op, coeff in zip(ham_real.paulis, ham_real.coeffs):
        print(f"{op}: {coeff}")
    return results_history



    