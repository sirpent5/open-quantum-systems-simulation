from imports import *
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

def build_initial_states(ham_real):

    ansatz = EfficientSU2(ham_real.num_qubits, reps = 8)

    init_param_values = {}
    for i in range(len(ansatz.parameters)):
        init_param_values[ansatz.parameters[i]] = np.random.uniform(0, 2 * np.pi)
        #init_param_values[ansatz.parameters[i]] = 2*np.pi

    init_state = Statevector(ansatz.assign_parameters(init_param_values))
    
    psi_vector = init_state.data


    initial_state = np.outer(psi_vector, psi_vector.conj())
   
    return init_state, initial_state, ansatz, init_param_values




def output_results(vqte_results, exact_diag_results, time, nt, time_points):
    """
    Plots a comparison of VQTE and exact diagonalization results.
    """
    plt.figure(figsize=(10, 6))
    
    # This time axis likely has nt + 1 points (e.g., 31 points)
    time_axis = np.linspace(0, time, nt + 1)

    # Plot Exact Diagonalization Results
    # Ensure the time axis slice matches the length of the results data
    for site_idx in range(len(exact_diag_results)):
        num_points = len(exact_diag_results[site_idx])
        plt.plot(time_axis[:num_points], 
                 exact_diag_results[site_idx], 
                 label=f'Exact Diag Site {site_idx} Occupation', 
                 marker='', 
                 linestyle='dashed')

    # Plot VQTE Results
    # Ensure the time axis slice matches the length of the results data
    for site_idx in range(len(vqte_results)):
        num_points = len(vqte_results[site_idx])
        plt.plot(time_axis[:num_points], 
                 vqte_results[site_idx], 
                 label=f'VQTE Site {site_idx} Occupation', 
                 marker='', 
                 linestyle='solid')

    plt.title("Comparison of VQTE and Exact Time Evolution")
    plt.xlabel("Time (t)")
    plt.ylabel("⟨n⟩ (Expectation Value)")
    plt.grid(True)
    plt.legend()
    
    plt.show()