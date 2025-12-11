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
    difference = superoperator + (1j*vqte_hamiltonian)
    return difference