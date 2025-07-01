
from imports import *
def hamiltonian_generation(eps, gamma, mu, T):
    """
    Generates the Hamiltonian for the system of a single qubit coupled to a reservoir.

    Inputs:
        eps: (float) Coupling strength of the qubit to the reservoir.
        gamma: (float) Coupling strength of the qubit to the environment.
        mu: (float) Chemical potential of the reservoir.
        T: (float) Temperature of the reservoir.
    Returns:
        hamiltonian_re: SparsePauliOp representing the real part of the Hamiltonian of the system.
        hamiltonian_im: SparsePauliOp representing the imaginary part of the Hamiltonian of the system.
    """
    
    F = 1 / (1 + np.exp((eps - mu) / T))
    hamiltonian_re = SparsePauliOp(["IZ", "ZI", "XY", "YX"], coeffs=[-eps / 2, eps / 2, -(gamma * (1 - 2*F)) / 4, -(gamma * (1 - 2*F)) / 4])
    hamiltonian_im = -1 * SparsePauliOp(["XX", "YY", "II", "IZ", "ZI"], coeffs=[gamma / 4, -gamma / 4, -gamma / 2, (gamma * (1 - 2*F)) / 4, (gamma * (1 - 2*F)) / 4])
    
    return hamiltonian_re, hamiltonian_im

def hamiltonian_generation_simple():
    """
    Generates a simple Hamiltonian for a single qubit system.

    Returns:
        hamiltonian_re: SparsePauliOp representing the Hamiltonian.
    """
    return SparsePauliOp(["IX", "XI"], coeffs=[1, -1])  # Example coefficients

def statevector_to_densitymatrix(v):
    """
    Converts a Statevector to a density matrix.

    Inputs:
        v: (numpy.ndarray) The state vector to be converted.
    
    Returns:
        (numpy.ndarray) The corresponding density matrix.
    """
    m = int(np.sqrt(len(v)))
    return np.reshape(v, (m, m), order='F')