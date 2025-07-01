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


def perform_exact_diag(gamma, F, dt, nt, initial_state, H):

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
    for step in range(1,nt+1):
        rho_t = U @ rho_t
        rho_matrix = rho_t.reshape(2 ,2)
        rho_matrix = rho_matrix / np.trace(rho_matrix)
        expectation_value_history.append(np.trace(numberop @ rho_matrix))
        time_points.append(step * dt)
    return expectation_value_history, time_points

def build_exact_diag_hamiltonian(eps):
    H = eps*Sigma_minus@Sigma_plus
    return H
