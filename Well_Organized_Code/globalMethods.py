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
    ansatz = EfficientSU2(ham_real.num_qubits, reps=1)
    init_param_values = {}
    for i in range(len(ansatz.parameters)):
        init_param_values[ansatz.parameters[i]] = (
        2*np.pi
    )  # initialize the parameters which also decide the initial state
    init_state = Statevector(ansatz.assign_parameters(init_param_values))
    
    psi_vector = init_state.data
    rho_matrix = psi_vector.reshape(2 ,2, order='F')


    initial_state = np.matrix(rho_matrix)
    return init_state, initial_state, ansatz, init_param_values


def output_results(vqte_results, exact_diag_results, time, nt, eps, mu, T, time_points, trace_list):
    plt.figure(figsize=(10, 6))
    print("Je")
    # Plot VQTE resultsW
    plt.plot(np.linspace(0, time, nt), trace_list, label='Trace of Density Matrix')
    plt.plot(np.linspace(0, time, nt), [1 / (1 + np.exp((eps - mu) / T))] * nt, label='Steady State Expectation Value', linestyle='--')
    plt.plot(time_points, exact_diag_results, label='Expectation Value (Simulated)', marker='', linestyle='-')
    # Plot Exact results
    #plt.plot(np.linspace(0, time, nt), exact_diag_results, marker='', linestyle='--', color='red', label='Exact Result')
    plt.plot(np.linspace(0, time, nt)/2, vqte_results,marker='', linestyle='-', label='VQTE Result')
    plt.title("Comparison of VQTE and Exact Time Evolution 🎯")
    plt.xlabel("Time (t)")
    plt.ylabel("⟨n⟩ (Expectation Value)")
    plt.grid(True)
    plt.legend()
    
    plt.show()