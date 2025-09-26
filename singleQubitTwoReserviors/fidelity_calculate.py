
def calculateFidelity(vqte_results, exact_results):
    """
    Calculate the fidelity between VQTE results and exact results.

    Parameters:
    vqte_results (list of np.ndarray): List of density matrices from VQTE simulation.
    exact_results (list of np.ndarray): List of density matrices from exact simulation.

    Returns:
    list of float: Fidelity values at each time step.
    """
    fidelities = []
    for rho_vqte, rho_exact in zip(vqte_results, exact_results):
        sqrt_rho_vqte = scipy.linalg.sqrtm(rho_vqte)
        fidelity = np.real(np.trace(scipy.linalg.sqrtm(sqrt_rho_vqte @ rho_exact @ sqrt_rho_vqte)))**2
        fidelities.append(fidelity)
    return fidelities