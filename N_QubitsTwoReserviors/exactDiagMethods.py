from imports import *
from defs import numberop, Sigma_minus, Sigma_plus, Sigma_x, Sigma_y, Sigma_z
from globalMethods import verify_density_matrix

def Liouvillian(H, Ls, hbar = 1):
    d = len(H) # dimension of the system
    superH = -1j/hbar * ( np.kron(np.eye(d),H)-np.kron(H.T,np.eye(d)) ) # Hamiltonian part
    superL = sum([np.kron(L.conjugate(),L) 
                  - 1/2 * ( np.kron(np.eye(d),L.conjugate().T.dot(L)) +
                            np.kron(L.T.dot(L.conjugate()),np.eye(d)) 
                          ) for L in Ls])
    return superH + superL

def Enlarge_Matrix_site_j(j, N, matrix):
    # I⊗...⊗I⊗M⊗I...⊗I: Convert local operators into global operators.
    # j: site, starts in 0.
    
    M = np.eye(len(matrix))
    if j == 0: M = matrix
    
    for i in range(1,N):
        if i == j: M = np.kron(M, matrix)
        else: M = np.kron(M, np.eye(len(matrix)))        

    return M

def Correlation_Matrix_i_Matrix_j(i,j,N,matrix_i, matrix_j):
    # I⊗...⊗I⊗M⊗I...⊗I⊗M⊗I⊗I...⊗I

    M = np.eye(len(matrix_i))
    
    if j == 0: M = matrix_j
    elif i == 0: M = matrix_i
    
    for k in range(1,N):
        if k == j: M = np.kron(M, matrix_j)
        elif k == i: M = np.kron(M, matrix_i)
        else: M = np.kron(M, np.eye(len(matrix_i)))        

    return M

def S_Term(N, cte_list, SigmaMatrix):
    # I⊗...⊗I⊗ΔSigma⊗Sigma⊗I...⊗I: Sigma can be Sigmax, Sigmay or Sigmaz.
    # cte_list = [cte_L, cte_M, cte_R], can be J or ∆.

    Matrix_Sigma = np.zeros((2**N, 2**N))

    cte_L = cte_list[0]
    cte_M = cte_list[1]
    cte_R = cte_list[2]    
    
    for i in range(0,N-1):
        M = np.eye(len(SigmaMatrix))
        
        if i == 0: M = SigmaMatrix
       
        for j in range(1,N):
            if j == i or j == i + 1: M = np.kron(M, SigmaMatrix)
            else: M = np.kron(M, np.eye(len(SigmaMatrix)))        

        if i < N/2 - 1: cte = cte_L
        elif i > N/2 - 1: cte = cte_R
        else: cte = cte_M

        Matrix_Sigma = Matrix_Sigma + M*cte #cte can be ∆_i or J_i

    return Matrix_Sigma #∑ I⊗...⊗I⊗ΔSigma⊗Sigma⊗I...⊗I

def build_number_op_list(N):
    """
    Builds a list of number operators, where each element corresponds
    to the number operator acting on a specific site.
    """
    return [Enlarge_Matrix_site_j(j, N, numberop) for j in range(N)]

# def perform_exact_diag(gamma, F_L, F_R, dt, nt, initial_state,H,N,eps):

    
#     L_K = []
    
#     L_K.append(np.sqrt(gamma*(1-F_L))*Enlarge_Matrix_site_j(0, N, Sigma_minus))  
#     L_K.append(np.sqrt(gamma*F_L)*Enlarge_Matrix_site_j(0, N, Sigma_plus)) 
#     L_K.append(np.sqrt(gamma *(1-F_R))*Enlarge_Matrix_site_j(N-1, N, Sigma_minus))
#     L_K.append(np.sqrt(gamma * F_R)*Enlarge_Matrix_site_j(N-1, N, Sigma_plus))
     
#     Superoperator = Liouvillian(H, L_K)



#     # Create time evolution operator
#     d = len(H)
#     U = scipy.linalg.expm(Superoperator * dt)
#     rho_t = initial_state.reshape(d**2,1)  # Vectorized  state
#     number_ops = build_number_op_list(N)

#     expectation_value_history= [[] for qubit in range(N)]

#     time_points = [0]

#     rho_matrix = initial_state / np.trace(initial_state)
#     for site in range(N):
#         expectation_value_history[site].append(np.trace(number_ops[site] @ rho_matrix))



#     print("Initial expectation value of number operator:", expectation_value_history[0])
#     # Time evolution loop
#     for step in range(1,nt+1):
        
#         rho_t = U @ rho_t
#         rho_matrix = rho_t.reshape(d ,d)
#         rho_matrix = rho_matrix / np.trace(rho_matrix)

#         for site in range(N):
#             expectation_value_history[site].append(np.trace(number_ops[site] @ rho_matrix))
           
#         time_points.append(step * dt)
#     return expectation_value_history, time_points

def perform_exact_diag(gamma, F_L, F_R, dt, nt, initial_state, H, N, eps):
    try:
        # Debug: Check initial state validity
        print("\n=== Debug: Initial Checks ===")
        print(f"Initial state trace: {np.trace(initial_state)}")
        assert np.isclose(np.trace(initial_state), 1.0, atol=1e-6), "Initial state not trace 1!"
        assert np.allclose(initial_state, initial_state.conj().T), "Initial state not Hermitian!"
        eigvals = np.linalg.eigvalsh(initial_state)
        assert np.all(eigvals >= -1e-8), "Initial state has negative eigenvalues!"

        # Build Lindblad operators
        L_K = [
            np.sqrt(gamma*(1-F_L)) * Enlarge_Matrix_site_j(0, N, Sigma_minus),
            np.sqrt(gamma*F_L) * Enlarge_Matrix_site_j(0, N, Sigma_plus),
            np.sqrt(gamma*(1-F_R)) * Enlarge_Matrix_site_j(N-1, N, Sigma_minus),
            np.sqrt(gamma*F_R) * Enlarge_Matrix_site_j(N-1, N, Sigma_plus)
        ]

        # Debug: Check Hamiltonian Hermiticity
        assert np.allclose(H, H.conj().T), "Hamiltonian not Hermitian!"

        # Construct superoperator
        Superoperator = Liouvillian(H, L_K)
        d = len(H)
        
        # Debug: Check superoperator shape
        print(f"Superoperator shape: {Superoperator.shape} (expected: {d**2}x{d**2})")

        # Time evolution operator
        U = scipy.linalg.expm(Superoperator * dt)
        
        # Initialize state (column-major vectorization)
        rho_t = initial_state.reshape(d**2, 1, order='F')
        number_ops = build_number_op_list(N)
        expectation_value_history = [[] for _ in range(N)]
        time_points = [0]

        # Initial measurements
        rho_matrix = initial_state.copy()
        for site in range(N):
            expectation_value_history[site].append(np.real(np.trace(number_ops[site] @ rho_matrix)))

        print("\n=== Debug: Starting Time Evolution ===")
        
        # Time evolution loop
        for step in range(1, nt+1):
            try:
                rho_t = U @ rho_t
                rho_matrix = rho_t.reshape(d, d, order='F')
                
                # Debug: Check trace and properties at each step
                current_trace = np.trace(rho_matrix)
                if not np.isclose(current_trace, 1.0, atol=1e-6):
                    print(f"WARNING: Step {step} trace = {current_trace:.6f} (renormalizing)")
                    rho_matrix /= current_trace
                
                if not np.allclose(rho_matrix, rho_matrix.conj().T):
                    print(f"WARNING: Step {step} density matrix not Hermitian! Applying Hermitization.")
                    rho_matrix = 0.5 * (rho_matrix + rho_matrix.conj().T)
                
                eigvals = np.linalg.eigvalsh(rho_matrix)
                if np.any(eigvals < -1e-8):
                    print(f"WARNING: Step {step} has negative eigenvalues! Clipping to 0.")
                    eigvals = np.clip(eigvals, 0, None)
                    eigvecs = np.linalg.eigh(rho_matrix)[1]
                    rho_matrix = (eigvecs * eigvals) @ eigvecs.conj().T
                    rho_matrix /= np.trace(rho_matrix)  # Renormalize after clipping

                # Store results
                for site in range(N):
                    expectation_value_history[site].append(np.real(np.trace(number_ops[site] @ rho_matrix)))
                
                time_points.append(step * dt)
                
                # Progress reporting
                if step % max(1, nt//10) == 0:
                    print(f"Completed step {step}/{nt} (t = {step*dt:.2f}), max eigval: {np.max(eigvals):.4f}")
            
            except Exception as step_error:
                print(f"\n!!! Error at step {step} !!!")
                print(f"rho_matrix trace: {np.trace(rho_matrix)}")
                print(f"rho_matrix Hermitian check: {np.allclose(rho_matrix, rho_matrix.conj().T)}")
                if 'eigvals' in locals():
                    print(f"Eigenvalues: {np.sort(eigvals)}")
                raise step_error

        print("\n=== Debug: Time Evolution Completed ===")
        print(f"Final state trace: {np.trace(rho_matrix)}")
        print(f"Final state eigenvalues: {np.linalg.eigvalsh(rho_matrix)}")
        
        return expectation_value_history, np.array(time_points)

    except Exception as main_error:
        print("\n!!! Main Procedure Error !!!")
        print(str(main_error))
        print("\n=== Debug Info ===")
        if 'H' in locals():
            print(f"Hamiltonian shape: {H.shape}")
            print(f"H Hermitian check: {np.allclose(H, H.conj().T)}")
        if 'Superoperator' in locals():
            print(f"Superoperator shape: {Superoperator.shape}")
        if 'rho_matrix' in locals():
            print(f"Last rho_matrix trace: {np.trace(rho_matrix)}")
        raise
# def build_exact_diag_hamiltonian(N, j, eps):
#     """
#     Builds the exact diagonalization Hamiltonian for a chain of N qubits.
    
#     Parameters:
#     -----------
#     N : int
#         Number of qubits in the chain
#     j : float
#         Coupling strength between neighboring qubits
#     eps : float
#         On-site energy parameter
        
#     Returns:
#     --------
#     np.ndarray
#         The Hamiltonian matrix (2^N x 2^N)
        
#     Raises:
#     -------
#     ValueError
#         If inputs are invalid or Hamiltonian construction fails
#     """
#     H = np.zeros((2**N, 2**N), dtype=complex)

#     # for i in range(N - 1):
#     #     # Use the correlation matrix helper for a cleaner implementation
#     #     term1 = Correlation_Matrix_i_Matrix_j(i, i + 1, N, Sigma_plus, Sigma_minus)
#     #     term2 = Correlation_Matrix_i_Matrix_j(i, i + 1, N, Sigma_minus, Sigma_plus)
#     #     H += j * (term1 + term2)
        

#     # for i in range(N - 1):
#     #     term1 = Correlation_Matrix_i_Matrix_j(i, i + 1, N, Sigma_plus, Sigma_minus)
#     #     term2 = Correlation_Matrix_i_Matrix_j(i, i + 1, N, Sigma_minus, Sigma_plus)
#     #     H += j * (term1 + term2)
    
#     # # On-site energies (alternating signs)
#     # for i in range(N):
#     #     sign = (-1)**i
#     #     Z_i = Correlation_Matrix_i_Matrix_j(i, i, N, Sigma_plus, Sigma_minus) - \
#     #           Correlation_Matrix_i_Matrix_j(i, i, N, Sigma_minus, Sigma_plus)
#     #     H += sign * (eps / 2) * Z_i
#     try:
#         # Input validation
#         if not isinstance(N, int) or N <= 0:
#             raise ValueError(f"N must be positive integer, got {N}")
#         if not isinstance(j, (int, float)):
#             raise ValueError(f"j must be numeric, got {type(j)}")
#         if not isinstance(eps, (int, float)):
#             raise ValueError(f"eps must be numeric, got {type(eps)}")
            
#         print(f"\nBuilding Hamiltonian for N={N} qubits")
#         print(f"Parameters: j={j}, ε={eps}")
        
#         dim = 2**N
#         H = np.zeros((dim, dim), dtype=complex)
        
#         # Helper function to validate Hamiltonian terms
#         def validate_term(term, description):
#             if not term.shape == (dim, dim):
#                 raise ValueError(f"{description} has wrong shape {term.shape}, expected {(dim, dim)}")
#             if not np.allclose(term, term.conj().T):
#                 print(f"WARNING: {description} not Hermitian! Applying Hermitization.")
#                 term = 0.5 * (term + term.conj().T)
#             return term
        
#         # Nearest-neighbor coupling terms
#         print("\nConstructing coupling terms...")
#         for i in range(N - 1):
#             try:
#                 term1 = Correlation_Matrix_i_Matrix_j(i, i+1, N, Sigma_plus, Sigma_minus)
#                 term2 = Correlation_Matrix_i_Matrix_j(i, i+1, N, Sigma_minus, Sigma_plus)
                
#                 term1 = validate_term(term1, f"Term1 for sites ({i},{i+1})")
#                 term2 = validate_term(term2, f"Term2 for sites ({i},{i+1})")
                
#                 H += j * (term1 + term2)
#                 print(f"Added coupling terms for sites {i}-{i+1}")
                
#             except Exception as term_error:
#                 print(f"!!! Error constructing terms for sites {i}-{i+1} !!!")
#                 raise term_error
        
#         # On-site energy terms
#         print("\nConstructing on-site energy terms...")
#         for i in range(N):
#             try:
#                 sign = (-1)**i
#                 Z_i = Correlation_Matrix_i_Matrix_j(i, i, N, Sigma_plus, Sigma_minus) - \
#                       Correlation_Matrix_i_Matrix_j(i, i, N, Sigma_minus, Sigma_plus)
                
#                 Z_i = validate_term(Z_i, f"On-site term {i}")
                
#                 H += sign * (eps / 2) * Z_i
#                 print(f"Added on-site term for site {i} (sign: {'+' if sign > 0 else '-'})")
                
#             except Exception as onsite_error:
#                 print(f"!!! Error constructing on-site term for site {i} !!!")
#                 raise onsite_error
        
#         # Final Hamiltonian validation
#         print("\nValidating final Hamiltonian...")
#         if not np.allclose(H, H.conj().T):
#             print("WARNING: Final Hamiltonian not Hermitian! Applying Hermitization.")
#             H = 0.5 * (H + H.conj().T)
        
#         # Check if Hamiltonian is all zeros (construction likely failed)
#         if np.allclose(H, 0):
#             raise ValueError("Constructed Hamiltonian is all zeros - likely construction error")
        
#         print("✓ Hamiltonian construction successful")
#         print(f"Final Hamiltonian shape: {H.shape}")
#         print(f"Norm: {np.linalg.norm(H):.4f}")
#         print(f"Trace: {np.trace(H):.4f}")
        
#         return H
        
#     except Exception as e:
#         print("\n!!! Hamiltonian construction failed !!!")
#         print(f"Error type: {type(e).__name__}")
#         print(f"Error message: {str(e)}")
        
#         # Debug info
#         if 'H' in locals():
#             print("\nPartial Hamiltonian info:")
#             print(f"Shape: {H.shape}")
#             print(f"Norm: {np.linalg.norm(H)}")
#             print("Matrix elements:")
#             print(H)
        
#         raise ValueError("Failed to construct Hamiltonian") from e

# def build_exact_diag_hamiltonian(N, j, eps):

#     H = np.zeros((2**N, 2**N), dtype=complex)

#     # for i in range(N - 1):
#     #     # Use the correlation matrix helper for a cleaner implementation
#     #     term1 = Correlation_Matrix_i_Matrix_j(i, i + 1, N, Sigma_plus, Sigma_minus)
#     #     term2 = Correlation_Matrix_i_Matrix_j(i, i + 1, N, Sigma_minus, Sigma_plus)
#     #     H += j * (term1 + term2)
        

#     for i in range(N - 1):
#         term1 = Correlation_Matrix_i_Matrix_j(i, i + 1, N, Sigma_plus, Sigma_minus)
#         term2 = Correlation_Matrix_i_Matrix_j(i, i + 1, N, Sigma_minus, Sigma_plus)
#         H += j * (term1 + term2)
    
#     # On-site energies (alternating signs)
#     for i in range(N):
#         sign = (-1)**i
#         Z_i = Correlation_Matrix_i_Matrix_j(i, i, N, Sigma_plus, Sigma_minus) - \
#               Correlation_Matrix_i_Matrix_j(i, i, N, Sigma_minus, Sigma_plus)
#         H += sign * (eps / 2) * Z_i
     
    
#     return H
def build_exact_diag_hamiltonian(N, j, eps):

    dim = 2**N
    H = np.zeros((dim, dim), dtype=complex)
    
    # Nearest-neighbor coupling terms
    for i in range(N - 1):
        term1 = Correlation_Matrix_i_Matrix_j(i, i+1, N, Sigma_plus, Sigma_minus)
        term2 = Correlation_Matrix_i_Matrix_j(i, i+1, N, Sigma_minus, Sigma_plus)
        H += j * (term1 + term2)
    
    # On-site energy terms (alternating signs)
    for i in range(N):
        sign = (-1)**i
        Z_i = Correlation_Matrix_i_Matrix_j(i, i, N, Sigma_plus, Sigma_minus) - \
              Correlation_Matrix_i_Matrix_j(i, i, N, Sigma_minus, Sigma_plus)
        H += sign * (eps / 2) * Z_i
    
    # Ensure Hermiticity
    H = 0.5 * (H + H.conj().T)
    
    return H
def output_exact_diag_results(exact_diag_results, time, nt, eps, mu_L,mu_R,T_L, T_R, time_points):

    plt.figure(figsize=(10, 6))

    for site_idx in range(len(exact_diag_results)): # Iterate through each site's data
        plt.plot(time_points, exact_diag_results[site_idx], label=f'Site {site_idx} Occupation', marker='', linestyle='solid')

    plt.title("Exact Diagonalization Results of two reserviors coupled to N qubits")
    plt.xlabel("Time (t)")
    plt.ylabel("⟨n⟩ (Expectation Value)")
    plt.grid(True)
    plt.legend()
    
    plt.show()


