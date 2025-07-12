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
    # initial_state = np.outer(psi_vector, np.conj(psi_vector))
    rho_matrix = psi_vector.reshape(2**N ,2**N, order='F')


    initial_state = np.matrix(rho_matrix)
    #For exact diag

    return init_state, initial_state, ansatz, init_param_values
def output_vqte_results(vqte_results, time, nt):
    plt.figure(figsize=(10, 6))
    time_axis = np.linspace(0, time, nt + 1)

    vqte_results_for_plot = np.asarray(vqte_results).T
   
    plt.plot(time_axis, vqte_results_for_plot, linestyle='dashed', label='VQTE Result(s)')
    plt.title("Comparison of VQTE and Exact Time Evolution")
    plt.xlabel("Time (t)")
    plt.ylabel("⟨n⟩ (Expectation Value)")
    plt.grid(True)
    plt.legend()
    
    plt.show()

def output_results(vqte_results, exact_diag_results, time, nt, time_points):
    plt.figure(figsize=(10, 6))
    time_axis = np.linspace(0, time, nt + 1)

    #exact_results_for_plot = np.asarray(exact_diag_results).T
    vqte_results_for_plot = np.asarray(vqte_results).T
    plt.plot(time_axis, vqte_results_for_plot, linestyle='dashed', label='VQTE Result(s)')
    for site_idx in range(len(exact_diag_results)): # Iterate through each site's data
        if(site_idx % 2 == 0):
            plt.plot(time_points, exact_diag_results[site_idx], label=f'Site {site_idx} Occupation', marker='', linestyle='solid')
        else:
            plt.plot(time_points, exact_diag_results[site_idx], label=f'Site {site_idx} Occupation', marker='', linestyle='dashed')
    plt.title("Comparison of VQTE and Exact Time Evolution")
    plt.xlabel("Time (t)")
    plt.ylabel("⟨n⟩ (Expectation Value)")
    plt.grid(True)
    plt.legend()
    
    plt.show()



def update_live_plot(expectation_value_history, time_points, N):
    """
    Clears the current plot and redraws it with the updated data.
    """
    plt.clf()  # Clear the current figure to prepare for the new frame
    
    for site_idx in range(N):
        # Plot the history for each site up to the current time
        plt.plot(time_points, expectation_value_history[site_idx], label=f'Site {site_idx}', marker='o', markersize=3, linestyle='-')

    plt.title("Live Exact Diagonalization: Qubit Occupation")
    plt.xlabel("Time (t)")
    plt.ylabel("⟨n⟩ (Expectation Value)")
    plt.grid(True)
    plt.legend()
    
    # Set fixed y-axis limits for stability, e.g., from 0 to 1 for occupation
    plt.ylim(-0.05, 1.05) 
    
    plt.draw()
    plt.pause(0.01) 