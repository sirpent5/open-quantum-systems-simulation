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


    """
    Builds Initial States for Exact Diagonalization and VQTE simulations.
    
    This function creates:
    1. A parameterized quantum circuit (ansatz) for VQTE
    2. Initial parameter values for the ansatz
    3. The corresponding quantum statevector
    4. A density matrix representation for exact diagonalization
    
    Inputs:
        ham_real : The real Hamiltonian for the VQTE simulation (determines number of qubits)
    
    Returns:
        init_state : Statevector - Initial quantum state for VQTE
        initial_state : np.matrix - Density matrix representation for exact diagonalization
        ansatz : QuantumCircuit - Parameterized circuit used for VQTE
        init_param_values : dict - Dictionary of initial parameter values for the ansatz
    """
    # Create an ansatz circut with reps
    ansatz = EfficientSU2(ham_real.num_qubits, reps = 1)

    #Initialize param dictionary
    init_param_values = {}

    # Set all params to 2π initially
    for i in range(len(ansatz.parameters)):
        init_param_values[ansatz.parameters[i]] = (
        2*np.pi)
      
    # Assign params to the ansatz
    vqte_init_state = Statevector(ansatz.assign_parameters(init_param_values))
    

    # Copy initial state data to a vector
    psi_vector = vqte_init_state.data

    # Reshape to a matrix
    rho_matrix = psi_vector.reshape(2 ,2, order='F')
    exact_diag_initial_state = np.matrix(rho_matrix)

    return vqte_init_state, exact_diag_initial_state, ansatz, init_param_values

def output_results(vqte_results, exact_diag_results, time, nt,time_points):
    plt.figure(figsize=(10, 6))
    time_axis = np.linspace(0, time, nt+1)
 
    plt.plot(time_points, exact_diag_results, label='Exact Diag Result', marker='', linestyle='solid')
    plt.plot(time_axis, vqte_results,marker='', linestyle='dashed', label='VQTE Result', color='blue')


    plt.title("Comparison of VQTE and Exact Time Evolution for a Qubit coupled to two thermal reserviors")
    plt.xlabel("Time (t)")
    plt.ylabel("⟨n⟩ (Expectation Value)")
    plt.grid(True)
    plt.legend()
    
    plt.show()