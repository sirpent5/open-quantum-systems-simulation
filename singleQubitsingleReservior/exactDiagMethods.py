# Imports
from imports import *
from defs import numberop, Sigma_minus, Sigma_plus
from globalMethods import verify_density_matrix

def Liouvillian(H, Ls, hbar = 1):
    d = len(H) # dimension of the system
    superH = -1j/hbar * ( np.kron(np.eye(d),H)-np.kron(H.T,np.eye(d)) ) # Hamiltonian part
    superL = sum([np.kron(L.conjugate(),L) 
                  - 1/2 * ( np.kron(np.eye(d),L.conjugate().T.dot(L)) +
                            np.kron(L.T.dot(L.conjugate()),np.eye(d)) 
                          ) for L in Ls])
    return superH + superL


def perform_exact_diag(gamma, F, dt, nt, initial_state, H):

    """
    Performs exact diagonalization of a single qubit system with Lindblad dynamics.
    
    Parameters:
        gamma (float): Decay rate/damping coefficient
        F (float): Fermi-Dirac distribution value (occupation probability)
        dt (float): Time step size
        nt (int): Number of time steps
        initial_state (np.array): Initial density matrix (2x2)
        H (np.array): System Hamiltonian (2x2)
        
    Returns:
        tuple: (expectation_value_history, time_points)
            - expectation_value_history: List of expectation values at each time step
            - time_points: List of corresponding time points
    """
    

    #Define lindblad operators
    L_plus = np.sqrt(gamma*(1-F)) * Sigma_plus
    L_minus = np.sqrt(gamma*F) * Sigma_minus
    L_K = [L_minus, L_plus] 

    # Create Superoperator
    Superoperator = Liouvillian(H, L_K)

    # verify_density_matrix(rho_ss)
    verify_density_matrix(initial_state)

    # Create time evolution operator
    U = scipy.linalg.expm(Superoperator * dt)
    rho_t = initial_state.reshape(4,1)  # Vectorized  state

    # Get inital values
    expectation_value_history = [np.trace(numberop @ initial_state) / np.trace(initial_state)]
    print("Initial expectation value of number operator:", expectation_value_history[0])
    time_points = [0]

    # Time evolution loop
    for step in range(1,nt+1):
        rho_t = U @ rho_t

        # Reshape into a density matrix and normalize
        rho_matrix = rho_t.reshape(2 ,2)
        rho_matrix = rho_matrix / np.trace(rho_matrix)

        # Calulate and store new values
        expectation_value_history.append(np.trace(numberop @ rho_matrix))
        time_points.append(step * dt)
    return expectation_value_history, time_points, Superoperator



def build_exact_diag_hamiltonian(eps):

    """
    Constructs the Hamiltonian for exact diagonalization of a two-level system (qubit).
    
    The Hamiltonian represents the energy of the excited state, with:
    H = ε|1⟩⟨1| = εσ₊σ₋
    where |1⟩ is the excited state and ε is its energy.

    Parameters:
        eps (float): The energy splitting/level spacing between ground |0⟩ and excited |1⟩ states
        
    Returns:
        numpy.ndarray: The 2×2 Hamiltonian matrix for the qubit system
    """

    H = eps*Sigma_minus@Sigma_plus
    return H
