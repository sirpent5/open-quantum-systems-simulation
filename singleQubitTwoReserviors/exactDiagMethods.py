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




def perform_exact_diag(gamma_L, F_L, gamma_R, F_R, dt, nt, initial_state, H):
    """
    Performs exact diagonalization for a quantum system coupled to two reservoirs.

    Args:
        gamma_L (float): Coupling strength for the left reservoir.
        F_L (float): Fermi-Dirac distribution parameter for the left reservoir.
        gamma_R (float): Coupling strength for the right reservoir.
        F_R (float): Fermi-Dirac distribution parameter for the right reservoir.
        dt (float): Time step for the evolution.
        nt (int): Number of time steps.
        initial_state (np.ndarray): The initial density matrix (2x2).
        H (np.ndarray): The system Hamiltonian (2x2).
    
    Returns:
        tuple: A tuple containing:
            - expectation_value_history (list): The expectation value of the number operator at each time step.
            - time_points (list): The list of time points.
    """
    # Define Lindblad operators for the Left Reservoir
    L_plus_L = np.sqrt(gamma_L * (1 - F_L)) * Sigma_plus
    L_minus_L = np.sqrt(gamma_L * F_L) * Sigma_minus

 
    L_plus_R = np.sqrt(gamma_R * (1 - F_R)) * Sigma_plus
    L_minus_R = np.sqrt(gamma_R * F_R) * Sigma_minus

   
    L_K = [L_minus_L, L_plus_L, L_minus_R, L_plus_R]
    
    # Construct the Liouvillian superoperator
    Superoperator = Liouvillian(H, L_K)
    
    # Calculate the steady state (null space of the Liouvillian)
    null = null_space(Superoperator)
    NULL = null[:, 0]
    rho_ss = NULL.reshape(2, 2)
    rho_ss = rho_ss / np.trace(rho_ss)

    verify_density_matrix(initial_state)

    # Create time evolution operator
    U = scipy.linalg.expm(Superoperator * dt)
    rho_t = initial_state.reshape(4, 1)  # Vectorized state

    # Initialize history lists
    expectation_value_history = [np.trace(numberop @ initial_state) / np.trace(initial_state)]
    print(f"Initial expectation value of number operator: {expectation_value_history[0]}")
    time_points = [0]

    # Time evolution loop
    for step in range(1, nt + 1):
        rho_t = U @ rho_t
        rho_matrix = rho_t.reshape(2, 2)
        
        # Ensure the trace is 1 at each step to prevent numerical drift
        rho_matrix = rho_matrix / np.trace(rho_matrix)
        
        expectation_value_history.append(np.trace(numberop @ rho_matrix))
        time_points.append(step * dt)
        
    return expectation_value_history, time_points







def build_exact_diag_hamiltonian(J, epsilon):
    N = len(epsilon)

    dim = 2**N
    H = np.zeros((dim, dim), dtype=complex)
    
    # On-site energy terms (ε a_j^† a_j)
    for j in range(N):
        # a_j^† a_j = (σ_j^- σ_j^+) = (1 - σ_j^z)/2
        Z_j = Enlarge_Matrix_site_j(j, N, Sigma_z)
        H += epsilon[j] * 0.5 * (np.eye(dim) - Z_j)
    
    # Hopping terms 
    for j in range(N-1):
        H += J*Correlation_Matrix_i_Matrix_j(j,j+1,N, Sigma_x, Sigma_x)
        H += J*Correlation_Matrix_i_Matrix_j(j,j+1,N, Sigma_y, Sigma_y) 
    return H
 



