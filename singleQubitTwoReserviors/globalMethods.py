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




def build_initial_states(ham_real, reps):


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
    ansatz = EfficientSU2(ham_real.num_qubits, reps)

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


def output_results(
    vqte_results, 
    exact_diag_results, 
    time_points, 
    mu_left, temp_left, gamma_left, 
    mu_right, temp_right, gamma_right
):
    """
    Generates a presentation-quality plot with dynamic axes and improved legend/parameter positioning.
    """
    # Use a professional style for the plot
    plt.style.use('seaborn-v0_8-talk')
    fig, ax = plt.subplots(figsize=(12, 7))

    # Assuming VQTE results correspond to a linspace from 0 to the max time
    time_axis_vqte = np.linspace(0, np.max(time_points), len(vqte_results))


    ax.plot(time_points, exact_diag_results, label='Exact Diagonalization', linestyle='-', linewidth=4, color='#003594')
    ax.plot(time_axis_vqte, vqte_results, label='VQTE', linestyle='--', linewidth=2.5, color="#FFB81C")

    # --- Set titles and labels ---
    ax.set_title("Population of a Single Qubit Coupled to Two Reservoirs", fontsize=20, pad=20)
    ax.set_xlabel("Time (t)", fontsize=20)
    ax.set_ylabel("$\\langle n \\rangle$", fontsize=20)
    
    # --- Create text for reservoir parameters ---
    params_text = (
        f"$\\bf{{Left~Reservoir}}$\n"
        f"$\\mu_L = {mu_left:.2f}$ | $T_L = {temp_left:.2f}$ | $\\gamma_L = {gamma_left:.2f}$\n\n"
        f"$\\bf{{Right~Reservoir}}$\n"
        f"$\\mu_R = {mu_right:.2f}$ | $T_R = {temp_right:.2f}$ | $\\gamma_R = {gamma_right:.2f}$"
    )
    
    # Add parameters text box with light background
    param_box = ax.text(0.97, 0.12, params_text, transform=ax.transAxes, fontsize=20,
                       verticalalignment='bottom', horizontalalignment='right',
                       bbox=dict(facecolor='white', edgecolor = 'white', alpha=0.8,boxstyle='round,pad=0.5'))
    
    # Add legend just above the parameters box
    legend = ax.legend(fontsize=20, loc='lower right', 
                      bbox_to_anchor=(0.97, 0.4),  # Positioned above the params box
                      frameon=True, framealpha=0.8,
                      facecolor='white', edgecolor='lightgray')
    
    # --- Finalize plot ---
    # Dynamically set axis limits
    x_max = np.max(time_points)
    y_max = np.max(np.concatenate([vqte_results, exact_diag_results]))
    
    ax.set_xlim(0, x_max) 
    ax.set_ylim(bottom=0, top=y_max * 1.1)
    
    ax.grid(False)
    ax.tick_params(axis='both', which='major', labelsize=20)  # Adjust font size of axis numbers
    fig.tight_layout()
    plt.show()