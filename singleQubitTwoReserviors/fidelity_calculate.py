import numpy as np
import matplotlib.pyplot as plt
from imports import*


def calculate_fidelity(vqte_results, exact_results):
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

def extract_density_matrix_components(results_list):
    """
    Extracts the four components (a, b, c, d) of 2x2 density matrices
    for every time step.

    Parameters:
    results_list (list of np.ndarray): List of 2x2 density matrices (rho).

    Returns:
    dict: Dictionary containing lists for 'a', 'b', 'c', 'd' components.
    """
    a_list, b_list, c_list, d_list = [], [], [], []
    

    for rho in results_list:
        a_list.append(np.real(rho[0, 0]))
        b_list.append(np.real(rho[0, 1])) # Note: b and c should generally be complex conjugates
        c_list.append(np.real(rho[1, 0]))
        d_list.append(np.real(rho[1, 1]))
        
    return {'a': a_list, 'b': b_list, 'c': c_list, 'd': d_list}

def plot_multiple_fidelity_vs_layers(all_results, common_params):
    """
    Plots multiple fidelity vs layers results on the same axes.

    Args:
        all_results (list): A list of dictionaries, where each dict has:
                            'label': str (for the legend)
                            'layers_list': list[int]
                            'fidelity_results': list[float]
        common_params (dict): Dictionary of common simulation parameters 
                              to be displayed on the plot.
    """
    plt.style.use('seaborn-v0_8-talk')
    fig, ax = plt.subplots(figsize=(10, 6))

    # Define a color palette (you can extend this for more lines)
    colors = ['#003594', '#940000', '#009435', '#946C00', '#6C0094']
    
    # 1. Loop through all results and plot each line
    for i, result in enumerate(all_results):
        label = result['label']
        layers_list = result['layers_list']
        fidelity_results = result['fidelity_results']
        color = colors[i % len(colors)] # Cycle through colors
        
        ax.plot(
            layers_list, 
            fidelity_results, 
            'o-', 
            linewidth=2, 
            markersize=8, 
            color=color,
            label=label # Set label for the legend
        )
    
    # 2. Add titles, labels, and legend
    ax.set_title("Fidelity vs Ansatz Layers for Different Scenarios", fontsize=20, pad=20)
    ax.set_xlabel("Number of Ansatz Layers", fontsize=16)
    ax.set_ylabel("Fidelity", fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    ax.legend(loc='lower right', fontsize=12, title='Scenario')

    # 3. Add common parameter info
    params_text = (
        f"$\\gamma_L = {common_params['gamma_L']:.1f}$, $\\gamma_R = {common_params['gamma_R']:.1f}$\n"
        f"$\\mu_L = {common_params['mu_L']:.1f}$, $\\mu_R = {common_params['mu_R']:.1f}$\n"
        f"$T_L = {common_params['T_L']:.1f}$, $T_R = {common_params['T_R']:.1f}$\n"
        f"$\\Delta t = {common_params['dt']:.2f}$, Time = {common_params['time']:.1f}$"
    )
    
    ax.text(0.02, 0.98, params_text, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', horizontalalignment='left',
            bbox=dict(facecolor='white', edgecolor='lightgray', alpha=0.8, boxstyle='round,pad=0.5'))
    
    # 4. Final aesthetics
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.05)
    plt.tight_layout()
    plt.show()


