from imports import *
from defs import numberop, Sigma_minus, Sigma_plus, Sigma_x, Sigma_y, Sigma_z

def Liouvillian(H, Ls, hbar = 1):
    d = len(H) # dimension of the system
    superH = -1j/hbar * ( np.kron(np.eye(d),H)-np.kron(H.T,np.eye(d)) ) # Hamiltonian part
    superL = sum([np.kron(L.conjugate(),L) 
                  - 1/2 * ( np.kron(np.eye(d),L.conjugate().T.dot(L)) +
                            np.kron(L.T.dot(L.conjugate()),np.eye(d)) 
                          ) for L in Ls])
    return superH + superL

def verify_density_matrix(rho):
    # Check Hermitian
    hermitian = np.allclose(rho, rho.conj().T)
    print(f"Is Hermitian: {hermitian}")
    
    # Check trace is 1
    trace = np.trace(rho)
    print(f"Trace: {trace} (should be 1)")
    
    # Check positive semidefinite?
    eigenvalues = np.linalg.eigvalsh(rho)
    print(f"Eigenvalues: {eigenvalues}")
    print(f"All eigenvalues ≥ 0: {np.all(eigenvalues >= 0)}")
    
    # Check purity 
    purity = np.trace(rho @ rho)
    print(f"Purity (Tr(ρ²)): {purity} (should be 1 for pure state)")



def perform_exact_diag(eps, gamma, F, dt, nt, initial_state):
    H = eps*Sigma_minus@Sigma_plus

    #Define lindblad operators
    L_plus = np.sqrt(gamma*(1-F)) * Sigma_plus
    L_minus = np.sqrt(gamma*F) * Sigma_minus

    L_K = [L_minus, L_plus] 



    Superoperator = Liouvillian(H, L_K)
    null = null_space(Superoperator)
    NULL = null[:, 0]
    rho_ss = NULL.reshape(2, 2)
    rho_ss = rho_ss / np.trace(rho_ss)

    referenceN = np.trace(numberop @ rho_ss)
    print(f"Reference number operator expectation value: {referenceN}")

    # verify_density_matrix(rho_ss)
    verify_density_matrix(initial_state)

    # Create time evolution operator
    U = scipy.linalg.expm(Superoperator * dt)
    rho_t = initial_state.reshape(4,1)  # Vectorized  state

    expectation_value_history = [np.trace(numberop @ initial_state) / np.trace(initial_state)]
    print("Initial expectation value of number operator:", expectation_value_history[0])
    time_points = [0]

    # Time evolution loop
    for step in range(nt):
        rho_t = U @ rho_t
        rho_matrix = rho_t.reshape(2 ,2)
        rho_matrix = rho_matrix / np.trace(rho_matrix)
        expectation_value_history.append(np.trace(numberop @ rho_matrix))
        time_points.append((step + 1) * dt)

def build_exact_diag_hamiltonian():
    print("Gonna do it one day")
def output_results(vqte_results, exact_diag_results, time, nt, eps, mu, T):
    plt.figure(figsize=(10, 6))
    
    # Plot VQTE results
    plt.plot(np.linspace(0, time, nt), vqte_results, marker='o', linestyle='-', label='VQTE Result')
    plt.plot(np.linspace(0, time, nt), [1 / (1 + np.exp((eps - mu) / T))] * nt, label='Steady State Expectation Value', linestyle='--')

    # Plot Exact results
    plt.plot(np.linspace(0, time, nt), exact_diag_results, marker='', linestyle='--', color='red', label='Exact Result')
    
    plt.title("Comparison of VQTE and Exact Time Evolution 🎯")
    plt.xlabel("Time (t)")
    plt.ylabel("⟨n⟩ (Expectation Value)")
    plt.grid(True)
    plt.legend()
    
    plt.show()