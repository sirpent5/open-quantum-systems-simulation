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

def perform_exact_diag(gamma, F_L,F_R, dt, nt, initial_state, H,N,eps):
    num_ops_per_site = build_number_op_list(N)
    total_num_op = sum(num_ops_per_site)
    L_K = []
    
    L_K.append(np.sqrt(gamma*(1-F_L))*Enlarge_Matrix_site_j(0, N, Sigma_minus))  
    L_K.append(np.sqrt(gamma*F_L)*Enlarge_Matrix_site_j(0, N, Sigma_plus)) 
    L_K.append(np.sqrt(gamma *(1-F_R))*Enlarge_Matrix_site_j(N-1, N, Sigma_minus))
    L_K.append(np.sqrt(gamma * F_R)*Enlarge_Matrix_site_j(N-1, N, Sigma_plus))
     
    Superoperator = Liouvillian(H, L_K)

    null = null_space(Superoperator)
    NULL = null[:, 0]
    rho_ss = NULL.reshape(2**N, 2**N)
    rho_ss = rho_ss / np.trace(rho_ss)

    referenceN = np.trace((eps * total_num_op) @ rho_ss)
    print(f"Reference number operator expectation value: {referenceN}")

    # verify_density_matrix(rho_ss)
    verify_density_matrix(initial_state)

    # Create time evolution operator
    U = scipy.linalg.expm(Superoperator * dt)

    rho_t = initial_state.reshape(2**N * 2**N, 1)  # Vectorized  state
    time_points = [0.0]
    
    # Initialize history list
    expectation_value_history = [[] for _ in range(N)]
 #save num op for each site
 # add hopping
    # Document initial values    
    for site_idx, op in enumerate(num_ops_per_site):
        expectation_i = np.trace((eps * op) @ initial_state)
        expectation_value_history[site_idx].append(expectation_i.real)

    print("Initial expectation value of number operator:", expectation_value_history[0])


    # Time evolution loop
    for step in range(1,nt+1):
        rho_t = U @ rho_t
        rho_matrix = rho_t.reshape(2**N ,2**N)
        rho_matrix /= np.trace(rho_matrix)

        for site_idx, op in enumerate(num_ops_per_site):
            exp_val = np.trace((eps * op) @ rho_matrix)
            expectation_value_history[site_idx].append(exp_val.real)
        time_points.append(step * dt)


    return expectation_value_history, time_points



def build_exact_diag_hamiltonian(N, j, eps):

    # H = np.zeros((2**N, 2**N), dtype=complex)


    H = np.zeros((2**N, 2**N), dtype=complex)

    # for i in range(N - 1):
    #     # Use the correlation matrix helper for a cleaner implementation
    #     term1 = Correlation_Matrix_i_Matrix_j(i, i + 1, N, Sigma_plus, Sigma_minus)
    #     term2 = Correlation_Matrix_i_Matrix_j(i, i + 1, N, Sigma_minus, Sigma_plus)
    #     H += j * (term1 + term2)
        

    for i in range(N - 1):
        term1 = Correlation_Matrix_i_Matrix_j(i, i + 1, N, Sigma_plus, Sigma_minus)
        term2 = Correlation_Matrix_i_Matrix_j(i, i + 1, N, Sigma_minus, Sigma_plus)
        H += j * (term1 + term2)
    
    # On-site energies (alternating signs)
    for i in range(N):
        sign = (-1)**i
        Z_i = Correlation_Matrix_i_Matrix_j(i, i, N, Sigma_plus, Sigma_minus) - \
              Correlation_Matrix_i_Matrix_j(i, i, N, Sigma_minus, Sigma_plus)
        H += sign * (eps / 2) * Z_i
     
    
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


