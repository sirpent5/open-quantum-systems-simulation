
# from qiskit.quantum_info import SparsePauliOp
# import numpy as np
from imports import *

import numpy as np
from qiskit.quantum_info import SparsePauliOp

from qiskit.quantum_info import SparsePauliOp
import numpy as np

def jordan_wigner_hamiltonian(N, eps, gamma_L, gamma_R, F_L, F_R, t):
    """
    N-site fermionic chain with reservoirs (2N qubits, Jordan-Wigner).
    
    Args:
        N: Number of physical sites
        eps: On-site energy (list of length N)
        gamma_L/R: Reservoir couplings
        F_L/R: Reservoir occupations (0 ≤ F ≤ 1)
        t: Hopping amplitude
    """
    num_qubits = 2 * N
    pauli_re, coeffs_re = [], []  # Real part (coherent dynamics)
    pauli_im, coeffs_im = [], []  # Imaginary part (dissipative)

    # === Real Part ===
    # 1. On-site terms (Z on occupation qubits)
    for i in range(N):
        op = ['I'] * num_qubits
        op[2*i] = 'Z'  # Occupation qubit
        pauli_re.append(''.join(op))
        coeffs_re.append(0.5 * eps[i])

    # 2. Hopping terms (XX + YY between adjacent sites)
    for i in range(N-1):
        # Hopping between site i (qubit 2i+1) and i+1 (qubit 2i+2)
        xx = ['I'] * num_qubits
        xx[2*i+1] = 'X'  # Phase qubit of site i
        xx[2*i+2] = 'X'  # Occupation qubit of site i+1
        pauli_re.append(''.join(xx))
        coeffs_re.append(t)

        yy = ['I'] * num_qubits
        yy[2*i+1] = 'Y'
        yy[2*i+2] = 'Y'
        pauli_re.append(''.join(yy))
        coeffs_re.append(t)

    # 3. Reservoir couplings (left/right edges)
    # Left reservoir (site 0: qubits 0,1)
    xy = ['I'] * num_qubits
    xy[0], xy[1] = 'X', 'Y'
    pauli_re.append(''.join(xy))
    coeffs_re.append(-0.25 * gamma_L * (1 - 2*F_L))

    yx = ['I'] * num_qubits
    yx[0], yx[1] = 'Y', 'X'
    pauli_re.append(''.join(yx))
    coeffs_re.append(-0.25 * gamma_L * (1 - 2*F_L))

    # Right reservoir (site N-1: qubits 2N-2, 2N-1)
    xy = ['I'] * num_qubits
    xy[-2], xy[-1] = 'X', 'Y'
    pauli_re.append(''.join(xy))
    coeffs_re.append(-0.25 * gamma_R * (1 - 2*F_R))

    yx = ['I'] * num_qubits
    yx[-2], yx[-1] = 'Y', 'X'
    pauli_re.append(''.join(yx))
    coeffs_re.append(-0.25 * gamma_R * (1 - 2*F_R))

    # === Imaginary Part ===
    # 1. Dissipators (acting on reservoir edges)
    # Left (site 0)
    xx = ['I'] * num_qubits
    xx[0], xx[1] = 'X', 'X'
    pauli_im.append(''.join(xx))
    coeffs_im.append(-0.25 * gamma_L)

    yy = ['I'] * num_qubits
    yy[0], yy[1] = 'Y', 'Y'
    pauli_im.append(''.join(yy))
    coeffs_im.append(0.25 * gamma_L)

    # Right (site N-1)
    xx = ['I'] * num_qubits
    xx[-2], xx[-1] = 'X', 'X'
    pauli_im.append(''.join(xx))
    coeffs_im.append(-0.25 * gamma_R)

    yy = ['I'] * num_qubits
    yy[-2], yy[-1] = 'Y', 'Y'
    pauli_im.append(''.join(yy))
    coeffs_im.append(0.25 * gamma_R)

    # 2. Boundary chemical potentials
    # Left (site 0)
    z_left = ['I'] * num_qubits
    z_left[0] = 'Z'
    pauli_im.append(''.join(z_left))
    coeffs_im.append(-0.25 * gamma_L * (1 - 2*F_L))

    # Right (site N-1)
    z_right = ['I'] * num_qubits
    z_right[-2] = 'Z'  # Occupation qubit of last site
    pauli_im.append(''.join(z_right))
    coeffs_im.append(-0.25 * gamma_R * (1 - 2*F_R))

    return SparsePauliOp(pauli_re, coeffs=np.array(coeffs_re)), \
           SparsePauliOp(pauli_im, coeffs=np.array(coeffs_im))
# def hamiltonian_generation(n_sites, eps, gamma_L, gamma_R, F_L, F_R, t):
#     """
#     Returns H_re and H_im for a fermionic chain (1 qubit per site, no spin).
    
#     Args:
#         n_sites: Number of physical sites (equals number of qubits)
#         eps: List of on-site energies [eps_1, ..., eps_n]
#         gamma_L, gamma_R: Left/right reservoir couplings
#         F_L, F_R: Reservoir occupation factors (0 ≤ F ≤ 1)
#         t: Hopping amplitude
        
#     Returns:
#         H_re, H_im: Real and imaginary parts as SparsePauliOp
#     """
#     # Initialize Pauli lists and coefficients
#     pauli_re, coeffs_re = [], []
#     pauli_im, coeffs_im = [], []

#     # ===== Real Part (H_re) =====
#     #I like this 
#     N = 2 * n_sites

#     for i in range(n_sites):
#         j = 2*i
#         z_str = ['I']* N
#         z_str[j] = 'Z'
#         pauli_re.append(''.join(z_str))
#         coeffs_re.append(eps[i]/2)

#         z_str = ['I']* N
#         z_str[j+1] = 'Z'
#         pauli_re.append(''.join(z_str))
#         coeffs_re.append(-eps[i]/2)
#     #I like this too
#     # 2. Hopping terms (t)

#     for i in range(n_sites-1):
#         # XX term
#         xx = ['I']*N
#         xx[2*i+1], xx[2*i+2] = 'X', 'X'
#         pauli_re.append(''.join(xx))
#         coeffs_re.append(t)
        
#         # YY term
#         yy = ['I']*N
#         yy[2*i+1], yy[2*i+2] = 'Y', 'Y'
#         pauli_re.append(''.join(yy))
#         coeffs_re.append(t)
    
#     # 3. Reservoir-induced XY/YX terms (only for edge pairs)
#     # Left edge (sites 0-1)

#     xy = ['I']*N
#     xy[0], xy[1] = 'X', 'Y'
#     pauli_re.append(''.join(xy))
#     coeffs_re.append(-0.25*gamma_L*(1-2*F_L))
    
#     yx = ['I']*N
#     yx[0], yx[1] = 'Y', 'X'
#     pauli_re.append(''.join(yx))
#     coeffs_re.append(-0.25*gamma_L*(1-2*F_L))

# # Right edge (sites n-2, n-1)

#     xy = ['I']*N
#     xy[-2], xy[-1] = 'X', 'Y'
#     pauli_re.append(''.join(xy))
#     coeffs_re.append(-0.25*gamma_R*(1-2*F_R))
        
#     yx = ['I']*N
#     yx[-2], yx[-1] = 'Y', 'X'
#     pauli_re.append(''.join(yx))
#     coeffs_re.append(-0.25*gamma_R*(1-2*F_R))

# # ===== Imaginary Part (H_im) =====
# # 1. Dissipative XX/YY terms (edges only)
# # Left edge
#     xx = ['I']*N
#     xx[0], xx[1] = 'X', 'X'
#     pauli_im.append(''.join(xx))
#     coeffs_im.append(-gamma_L/4)
    
#     yy = ['I']*N
#     yy[0], yy[1] = 'Y', 'Y'
#     pauli_im.append(''.join(yy))
#     coeffs_im.append(gamma_L/4)
    
# # Right edge

#     xx = ['I']*N
#     xx[-2], xx[-1] = 'X', 'X'
#     pauli_im.append(''.join(xx))
#     coeffs_im.append(-gamma_R/4)
    
#     yy = ['I']*N
#     yy[-2], yy[-1] = 'Y', 'Y'
#     pauli_im.append(''.join(yy))
#     coeffs_im.append(gamma_R/4)
    
#     # 2. Global identity term
#     pauli_im.append('I'*N)
#     coeffs_im.append((gamma_L + gamma_R)/2)
    
#     # 3. Boundary Z terms
#     # Left reservoir
#     z_left = ['I']*N
#     z_left[0] = 'Z'
#     pauli_im.append(''.join(z_left))
#     coeffs_im.append(-0.25*gamma_L*(1-2*F_L))
    
#     # Right reservoir
#     z_right = ['I']*N
#     z_right[-1] = 'Z'
#     pauli_im.append(''.join(z_right))
#     coeffs_im.append(-0.25*gamma_R*(1-2*F_R))

#     return SparsePauliOp(pauli_re, coeffs=np.array(coeffs_re)), SparsePauliOp(pauli_im, coeffs=np.array(coeffs_im))
        

# def hamiltonian_generation(n_sites, eps, gamma_L, gamma_R, F_L, F_R, t):
#     """
#     Returns H_re and H_im for a fermionic chain (1 qubit per site, no spin).
    
#     Args:
#         n_sites: Number of physical sites (equals number of qubits)
#         eps: List of on-site energies [eps_1, ..., eps_n]
#         gamma_L, gamma_R: Left/right reservoir couplings
#         F_L, F_R: Reservoir occupation factors (0 ≤ F ≤ 1)
#         t: Hopping amplitude
        
#     Returns:
#         H_re, H_im: Real and imaginary parts as SparsePauliOp
#     """
#     # Initialize Pauli lists and coefficients

#     pauli_re, coeffs_re = [], []
#     pauli_im, coeffs_im = [], []
#     #     hamiltonian_re = SparsePauliOp( ["IZ", "ZI", "XY", "YX"],
#     #     coeffs=[(-eps/2),(eps/2), ((-1/4)*((gamma_L*(1-2*F_L))+(gamma_R*(1-2*F_R)))),((-1/4)*((gamma_L*(1-2*F_L))+(gamma_R*(1-2*F_R))))]
#     # )


#     # hamiltonian_im = SparsePauliOp( ["XX", "YY", "II", "IZ", "ZI"],
#     #     coeffs=[-(gamma_L+gamma_R)/4, (gamma_L+gamma_R)/4, (gamma_L+gamma_R)/2,
#     #              ((-1/4)*((gamma_L*(1-2*F_L))+(gamma_R*(1-2*F_R)))), ((-1/4)*((gamma_L*(1-2*F_L))+(gamma_R*(1-2*F_R)))) ] )



#     # ===== Real Part (H_re) =====
#     #I like this 
    
#     N  = n_sites * 2

#     for n in range(n_sites):
#         eps.append(1)
#     # hamiltonian_re = SparsePauliOp( ["IZ", "ZI", "XY", "YX"],
#     #     coeffs=[(-eps/2),(eps/2), ((-1/4)*((gamma_L*(1-2*F_L))+(gamma_R*(1-2*F_R)))),((-1/4)*((gamma_L*(1-2*F_L))+(gamma_R*(1-2*F_R))))]
#     # )


#     # hamiltonian_im = SparsePauliOp( ["XX", "YY", "II", "IZ", "ZI"],
#     #     coeffs=[-(gamma_L+gamma_R)/4, (gamma_L+gamma_R)/4, (gamma_L+gamma_R)/2,
#     #              ((-1/4)*((gamma_L*(1-2*F_L))+(gamma_R*(1-2*F_R)))), ((-1/4)*((gamma_L*(1-2*F_L))+(gamma_R*(1-2*F_R)))) ] )

#     for i in range(n_sites):
#         j = 2*i
#         z_str = ['I']* N
#         z_str[j] = 'Z'
#         pauli_re.append(''.join(z_str))
#         coeffs_re.append(eps[i]/2)

#         z_str = ['I']* N
#         z_str[j+1] = 'Z'
#         pauli_re.append(''.join(z_str))
#         coeffs_re.append(-eps[i]/2)

#     #I like this too
#     # 2. Hopping terms (t)

#     for i in range(n_sites-1):
#         # XX term

#         xx = ['I']*N
#         xx[i], xx[2*i+1] = 'X', 'X'
#         pauli_re.append(''.join(xx))
#         coeffs_re.append(t)
        
#         # YY term
#         yy = ['I']*N
#         yy[i], yy[2*i+1] = 'Y', 'Y'
#         pauli_re.append(''.join(yy))
#         coeffs_re.append(t)
    
#     # 3. Reservoir-induced XY/YX terms (only for edge pairs)
#     # Left edge (sites 0-1)
  
#     xy = ['I']* N
#     xy[0], xy[1] = 'X', 'Y'
#     pauli_re.append(''.join(xy))
#     coeffs_re.append(-0.25*gamma_L*(1-2*F_L))
    
#     yx = ['I']* N
#     yx[0], yx[1] = 'Y', 'X'
#     pauli_re.append(''.join(yx))
#     coeffs_re.append(-0.25*gamma_L*(1-2*F_L))

# # Right edge (sites n-2, n-1)

#     xy = ['I']* N
#     xy[-2], xy[-1] = 'X', 'Y'
#     pauli_re.append(''.join(xy))
#     coeffs_re.append(-0.25*gamma_R*(1-2*F_R))
        
#     yx = ['I']* N
#     yx[-2], yx[-1] = 'Y', 'X'
#     pauli_re.append(''.join(yx))
#     coeffs_re.append(-0.25*gamma_R*(1-2*F_R))

#         # print(SparsePauliOp(pauli_re, coeffs=np.array(coeffs_re)))

#     # ===== Imaginary Part (H_im) =====
#     # 1. Dissipative XX/YY terms (edges only)
#     # Left edge

#     xx = ['I']*N
#     xx[0], xx[1] = 'X', 'X'
#     pauli_im.append(''.join(xx))
#     coeffs_im.append(-gamma_L/4)
    
#     yy = ['I']*N
#     yy[0], yy[1] = 'Y', 'Y'
#     pauli_im.append(''.join(yy))
#     coeffs_im.append(gamma_L/4)
        
#     # Right edge

#     xx = ['I']*N
#     xx[-2], xx[-1] = 'X', 'X'
#     pauli_im.append(''.join(xx))
#     coeffs_im.append(-gamma_R/4)
    
#     yy = ['I']*N
#     yy[-2], yy[-1] = 'Y', 'Y'
#     pauli_im.append(''.join(yy))
#     coeffs_im.append(gamma_R/4)

#     # 2. Global identity term
#     pauli_im.append('I'*N)
#     coeffs_im.append((gamma_L + gamma_R)/2)
    
#     # 3. Boundary Z terms
#     # Left reservoir
#     z_left = ['I']*N
#     z_left[0] = 'Z'
#     pauli_im.append(''.join(z_left))
#     coeffs_im.append(-0.25*gamma_L*(1-2*F_L))
    
#     # Right reservoir
#     z_right = ['I']*N
#     z_right[-1] = 'Z'
#     pauli_im.append(''.join(z_right))
#     coeffs_im.append(-0.25*gamma_R*(1-2*F_R))

#     return SparsePauliOp(pauli_re, coeffs=np.array(coeffs_re)), SparsePauliOp(pauli_im, coeffs=np.array(coeffs_im))
        

    

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

# def perform_vqte(ham_real, ham_imag, init_state,dt, nt, ansatz, init_param_values,N):
#     real_var_principle = RealMcLachlanPrinciple(qgt=ReverseQGT(), gradient=ReverseEstimatorGradient(derivative_type=DerivativeType.IMAG))
#     imag_var_principle = ImaginaryMcLachlanPrinciple(qgt=ReverseQGT(), gradient=ReverseEstimatorGradient())

#     number_operators = [create_number_operator(N, i) for i in range(N)]

#     print("Number op: ",number_operators[0])

#     # Perform time evolution
#     results_history = [[] for _ in range(N)]
   
#    # num_op_list = [np.trace(statevector_to_densitymatrix(init_state.data)
#    #                          @ np.array([[0, 0], [0, 1]])) / np.trace(statevector_to_densitymatrix(init_state.data))]
#     # for i, op in enumerate(number_operators):
#     #         op_matrix = op.to_matrix()

#     #         initial_exp_val = [np.trace(statevector_to_densitymatrix(init_state.data)
#     #                                     @ op_matrix) / np.trace(statevector_to_densitymatrix(init_state.data))

#     #         ]
            
#     #         init_state.expectation_value(op).real
#     #         results_history[i].append(initial_exp_val)

#     for i, op in enumerate(number_operators):
        
#         op_matrix = op.to_matrix()
#         density_matrix = statevector_to_densitymatrix(init_state.data)


#         print("Number Op: ", op_matrix)
#         print("Density Matrix: ", density_matrix)
#         #density_matrix = np.where(np.abs(density_matrix) < 1e-10, 0, density_matrix)
#         initial_exp_val = np.trace(density_matrix @ op_matrix).real
    
#         results_history[i].append(initial_exp_val)
#         # initial_exp_val = (np.trace(density_matrix @ op_matrix) / np.trace(density_matrix)).real
        

   
#         print("This is the init valie: ",initial_exp_val)


#     #trace = np.trace(statevector_to_densitymatrix(Statevector(ansatz.assign_parameters(init_param_values)).data))
#     print(results_history[0])
#     # --- Time Evolution Loop ---
#     for t in range(nt):
#         print(f"Step {t+1} out of {nt}")
#         # Real and Imaginary evolution steps
#         evolution_problem_re = TimeEvolutionProblem(ham_real, dt)
#         var_qrte = VarQRTE(ansatz, init_param_values, real_var_principle, num_timesteps=1)
#         evolution_result_re = var_qrte.evolve(evolution_problem_re)
#         init_param_values = evolution_result_re.parameter_values[-1]

#         evolution_problem_im = TimeEvolutionProblem(ham_imag, dt)
#         var_qite = VarQITE(ansatz, init_param_values, imag_var_principle, num_timesteps=1)
#         evolution_result_im = var_qite.evolve(evolution_problem_im)
#         init_param_values = evolution_result_im.parameter_values[-1]



#         # --- Measurement ---
#         current_psi = Statevector(ansatz.assign_parameters(init_param_values))
#         norm = np.linalg.norm(current_psi.data)
#         current_psi = current_psi / norm if norm != 0 else current_psi
#         #current_psi = current_psi / np.linalg.norm(current_psi.data)
#         # Calculate expectation value for each qubit's number operator
#         for i, op in enumerate(number_operators):
#             op_matrix = op.to_matrix()
         
#             density_matrix = statevector_to_densitymatrix(current_psi.data)
#             #density_matrix = np.where(np.abs(density_matrix) < 1e-10, 0, density_matrix)
#             exp_val = np.trace(density_matrix @ op_matrix).real
#             results_history[i].append(exp_val)
#             # exp_val = (np.trace(statevector_to_densitymatrix(
#             # Statevector(ansatz.assign_parameters(init_param_values)).data) @ op_matrix) / trace)

            
      
#     return results_history
def perform_vqte(ham_real, ham_imag, init_state, dt, nt, ansatz, init_param_values, N):
    real_var_principle = RealMcLachlanPrinciple(qgt=ReverseQGT(), gradient=ReverseEstimatorGradient(derivative_type=DerivativeType.IMAG))
    imag_var_principle = ImaginaryMcLachlanPrinciple(qgt=ReverseQGT(), gradient=ReverseEstimatorGradient())

    number_operators = [create_number_operator(N, i) for i in range(N)]
    results_history = [[] for _ in range(N)]


    for i, op in enumerate(number_operators):
        op_matrix = op.to_matrix()
        density_matrix = statevector_to_densitymatrix(init_state.data)


        print("Number Op: ", op_matrix)
        print("Density Matrix: ", density_matrix)
  
        initial_exp_val = np.trace(density_matrix @ op_matrix).real
    
        results_history[i].append(initial_exp_val)

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
    
        for i, op in enumerate(number_operators):
            op_matrix = op.to_matrix()
            exp_val = np.trace(rho_matrix @ op_matrix)/true_trace
            results_history[i].append(exp_val)

    return results_history

