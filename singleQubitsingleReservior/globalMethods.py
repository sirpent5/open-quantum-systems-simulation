# Imports
from imports import *

def verify_density_matrix(rho):

    """
    Verifies that a Density Matrix is Hermatian

    Inputs:
        rho: (density matrix) The density matrix to be checked.
    
    """
        
    # Check Hermitian
    hermitian = np.allclose(rho, rho.conj().T)
    print(f"Is Hermitian: {hermitian}")
    
    # Check trace is 1
    trace = np.trace(rho)
    print(f"Trace: {trace} (should be 1)")
    
    # Check positive semidefinite
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


def output_results(vqte_results, exact_diag_results, time, nt, eps, mu, T, time_points):

    """
    Plots and compares the results of VQTE simulation, exact diagonalization, and steady state value.
    
    Parameters:
        vqte_results (array-like): Array of expectation values from VQTE simulation
        exact_diag_results (array-like): Array of expectation values from exact diagonalization
        time (float): Total simulation time
        nt (int): Number of time steps
        eps (float): Energy level of the qubit
        mu (float): Chemical potential of the reservoir
        T (float): Temperature of the reservoir
        time_points (array-like): Array of time points for exact diagonalization results
        
    Returns:
        None (displays a matplotlib plot)
    """

    # Set up chart
    plt.figure(figsize=(10, 6))
    time_axis = np.linspace(0, time, nt+1)

    # Plot exact values
    plt.plot(time_axis, [1 / (1 + np.exp((eps - mu) / T))] * (nt+1), label='Steady State Expectation Value', linestyle='solid')
    plt.plot(time_points, exact_diag_results, label='Exact Classical Solution', marker='', linestyle='solid')

    # Plot VQTE values
    plt.plot(time_axis, vqte_results,marker='', linestyle='dashed', label='VQTE Simulation Result', color='blue')

    # Labels
    plt.title("Comparison of VQTE and Exact Time Evolution for a Single Qubit coupled to a Single Reservoir")
    plt.xlabel("Time (t)")
    plt.ylabel("Expectation Value")
    plt.grid(True)
    plt.legend()
    
    plt.show()
    
def compare_superoperator_to_vqte(superoperator, ham_real, ham_imag):
    """
    Compares the superoperator obtained from exact diagonalization
    with the Hamiltonian constructed from VQTE components.

    Inputs:
        superoperator : darray - The superoperator from exact diagonalization
        ham_real : SparsePauliOp - The real part of the Hamiltonian from VQTE
        ham_imag : SparsePauliOp - The imaginary part of the Hamiltonian from VQTE

    Returns:
        difference : array - The difference between the two matrices
    """
    vqte_hamiltonian = ham_real.to_matrix()- 1j * ham_imag.to_matrix()
    difference = superoperator - (-1j*vqte_hamiltonian)
    return difference