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

def build_initial_states(ham_real, N):
    ansatz = EfficientSU2(ham_real.num_qubits, reps=1)
    print(ham_real.num_qubits)
    init_param_values = {}
    for i in range(len(ansatz.parameters)):
        init_param_values[ansatz.parameters[i]] = (
        2*np.pi
    )  # initialize the parameters which also decide the initial state
    init_state = Statevector(ansatz.assign_parameters(init_param_values))
    
    psi_vector = init_state.data

    rho_matrix = psi_vector.reshape(2**N ,2**N, order='F')
    initial_state = np.matrix(rho_matrix)



    # #initial_state = np.outer(psi_vector, psi_vector.conj())
    # state = np.zeros(2**N, dtype=complex)
    # state[1] = 1.0  # Binary '00...01' (first qubit = |1⟩)

    # initial_state = np.outer(state, state.conj())  #

    #For exact diag
    return init_state, initial_state, ansatz, init_param_values



def output_results(vqte_results, exact_diag_results, time, nt, time_points):
    plt.figure(figsize=(10, 6))
    time_axis = np.linspace(0, time, nt + 1)
    vqte_results_for_plot = np.asarray(vqte_results).T
    #exact_results_for_plot = np.asarray(exact_diag_results).T
    # for site_idx, vqte_data in enumerate(vqte_results_for_plot):
    #     plt.plot(time_axis, vqte_data, label=f'VQTE Site {site_idx}', linestyle='dashed', marker='x')
    for site_idx in range(len(exact_diag_results)): # Iterate through each site's data
            plt.plot(time_points, exact_diag_results[site_idx], label=f'Exact Diag Site {site_idx} Occupation', marker='', linestyle='dashed')

    for site_idx in range(len(vqte_results)): # Iterate through each site's data
            plt.plot(time_points, vqte_results[site_idx], label=f'VQTE Site {site_idx} Occupation', marker='', linestyle='solid')

    plt.title("Comparison of VQTE and Exact Time Evolution")
    plt.xlabel("Time (t)")
    plt.ylabel("⟨n⟩ (Expectation Value)")
    plt.grid(True)
    plt.legend()
    
    plt.show()

