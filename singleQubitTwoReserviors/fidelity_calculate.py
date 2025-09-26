import numpy as np
import matplotlib.pyplot as plt
from imports import*


def calculate_fidelity(vqte_results, exact_results):
    """
    Simple fidelity calculation with error handling.
    """
    fidelities = []
    
    for i in range(min(len(vqte_results), len(exact_results))):
        try:
            rho_vqte = vqte_results[i]
            rho_exact = exact_results[i]
            
            if (isinstance(rho_vqte, np.ndarray) and isinstance(rho_exact, np.ndarray) and
                rho_vqte.shape == (2, 2) and rho_exact.shape == (2, 2)):
                
                sqrt_rho_vqte = scipy.linalg.sqrtm(rho_vqte)
                product = sqrt_rho_vqte @ rho_exact @ sqrt_rho_vqte
                sqrt_product = scipy.linalg.sqrtm(product)
                fidelity = np.real(np.trace(sqrt_product))**2
                fidelities.append(fidelity)
            else:
                fidelities.append(0.0)
        except:
            fidelities.append(0.0)
    
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
def plot_multiple_fidelity_vs_layers(all_results):
    """
    Plots multiple fidelity vs layers results on the same axes.

    Args:
        all_results (list): A list of dictionaries, where each dict has:
                            'label': str (for the legend)
                            'layers_list': list[int]
                            'fidelity_results': list[float]
    """
    plt.style.use('seaborn-v0_8-talk')
    fig, ax = plt.subplots(figsize=(10, 6))

    # Define a color palette
    colors = ['#003594', '#940000', '#009435', '#946C00', '#6C0094']
    
    # Check if all_results is a list of dictionaries (multiple scenarios)
    # or if it's just a single list of fidelity results (backward compatibility)
    if not isinstance(all_results, list) or (all_results and not isinstance(all_results[0], dict)):
        # Convert single list to the expected format
        all_results = [{
            'label': 'Single Scenario',
            'layers_list': list(range(1, len(all_results) + 1)),
            'fidelity_results': all_results
        }]
    
    # Loop through all results and plot each line
    for i, result in enumerate(all_results):
        label = result['label']
        layers_list = result['layers_list']
        fidelity_results = result['fidelity_results']
        color = colors[i % len(colors)]
        
        ax.plot(
            layers_list, 
            fidelity_results, 
            'o-', 
            linewidth=2, 
            markersize=8, 
            color=color,
            label=label
        )
    
    # Add titles, labels, and legend
    ax.set_title("Fidelity vs Ansatz Layers for Different Scenarios", fontsize=20, pad=20)
    ax.set_xlabel("Number of Ansatz Layers", fontsize=16)
    ax.set_ylabel("Fidelity", fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    
    # Only show legend if there are multiple scenarios
    if len(all_results) > 1:
        ax.legend(loc='lower right', fontsize=12, title='Scenario')

    # Final aesthetics
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.05)
    ax.set_xlim(0.5, max(layers_list) + 0.5)  # Adjust x-axis limits
    plt.tight_layout()
    plt.show()