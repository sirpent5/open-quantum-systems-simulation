import numpy as np
import matplotlib.pyplot as plt
from imports import*
import cirq


def calculate_fidelity(vqte_results, exact_results):
    """
    Simple fidelity calculation with error handling.
    """
    fidelities = []
    
    for i in range(min(len(vqte_results), len(exact_results))):
        rho_vqte_dm = vqte_results[i], qid_shape=(2,)
        rho_exact_dm = exact_results[i], qid_shape=(2,)
        fidelity_value = cirq.fidelity(rho_vqte_dm, rho_exact_dm)

        # rho_vqte = vqte_results[i]
        # rho_exact = exact_results[i]
            
      
                        
        # sqrt_rho_vqte = scipy.linalg.sqrtm(rho_vqte)
        # product = sqrt_rho_vqte @ rho_exact @ sqrt_rho_vqte
        # sqrt_product = scipy.linalg.sqrtm(product)
        # fidelity = np.real(np.trace(sqrt_product))**2
        fidelities.append(fidelity_value)

    return fidelities

def extract_density_matrix_components(vqte_list, exact_list):
    """
    Extracts the four components (a, b, c, d) of 2x2 density matrices
    for every time step.

    Parameters:
    results_list (list of np.ndarray): List of 2x2 density matrices (rho).

    Returns:
    dict: Dictionary containing lists for 'a', 'b', 'c', 'd' components.
    """
    va_list, vb_list, vc_list, vd_list = [], [], [], []
    ea_list, eb_list, ec_list, ed_list = [], [], [], []

    for rho in vqte_list:
        va_list.append(np.real(rho[0, 0]))
        vb_list.append(np.real(rho[0, 1])) # Note: b and c should generally be complex conjugates
        vc_list.append(np.real(rho[1, 0]))
        vd_list.append(np.real(rho[1, 1]))
        
        va_list, vb_list, vc_list, vd_list = [], [], [], []
    

    for rho in exact_list:
        ea_list.append(np.real(rho[0, 0]))
        eb_list.append(np.real(rho[0, 1])) # Note: b and c should generally be complex conjugates
        ec_list.append(np.real(rho[1, 0]))
        ed_list.append(np.real(rho[1, 1]))

    return {'va': va_list, 'vb': vb_list, 'vc': vc_list, 'vd': vd_list, 'ea': ea_list, 'eb': eb_list, 'ec': ec_list, 'ed': ed_list}

# def plot_matrix_components(all_matrix_components_by_layer, time, nt, layers):
#     """
#     Plots the real parts of the four components of the density matrix 
#     (rho_00, rho_01, rho_10, rho_11) over time, comparing VQTE and Exact results
#     for different ansatz layers.

#     Args:
#         all_matrix_components_by_layer (list): List of dictionaries, where each dict 
#                                                contains {'layers': int, 'components': dict}.
#         time (float): Total simulation time.
#         nt (int): Number of time steps.
#     """
#     plt.style.use('seaborn-v0_8-talk')
#     fig, axes = plt.subplots(2, 2, figsize=(16, 12))
#     axes = axes.flatten()
#     time_axis = np.linspace(0, time, nt + 1)
    
#     # Titles and component keys (VQTE and Exact have same structure)
#     titles = [
#         r'Diagonal Element $\rho_{00}$ (Occupation of $\ket{0}$)',
#         r'Off-Diagonal Element $\text{Re}(\rho_{01})$',
#         r'Off-Diagonal Element $\text{Re}(\rho_{10})$', 
#         r'Diagonal Element $\rho_{11}$ (Occupation of $\ket{1}$)'
#     ]
#     # Keys for extraction: VQTE (v) and Exact (e)
#     component_v_keys = ['va', 'vb', 'vc', 'vd']
#     component_e_keys = ['ea', 'eb', 'ec', 'ed'] 

#     # Define a color palette for different layers (e.g., Viridis for VQTE lines)
#     num_layers = len(all_matrix_components_by_layer)
#     vqte_colors = plt.cm.plasma(np.linspace(0.2, 0.9, num_layers))
    
#     # Plotting loop
#     for i, data in enumerate(all_matrix_components_by_layer):

#         components = data['components']
        
#         # Plot VQTE lines (solid, colored based on layer)
#         for j in range(4):
#             ax = axes[j]
#             v_key = component_v_keys[j]
#             component_data = components[v_key]
            
#             num_points = len(component_data)
#             ax.plot(time_axis[:num_points], 
#                     component_data, 
#                     label=f'VQTE, {layers} Layers', 
#                     color=vqte_colors[i],
#                     linestyle='-',
#                     linewidth=2,
#                     alpha=0.7)

#         # Plot Exact lines (dashed, black/dark gray, only once per subplot)
#         if i == 0:
#             for j in range(4):
#                 ax = axes[j]
#                 e_key = component_e_keys[j]
#                 component_data = components[e_key]
                
#                 num_points = len(component_data)
#                 ax.plot(time_axis[:num_points], 
#                         component_data, 
#                         label=f'Exact', 
#                         color='black',
#                         linestyle='--',
#                         linewidth=3,
#                         alpha=0.9)


#     # Finalize plots
#     for ax, title in zip(axes, titles):
#         ax.set_title(title, fontsize=18, pad=10)
#         ax.set_xlabel('Time (t)', fontsize=16)
#         ax.set_ylabel('Real Value', fontsize=16)
#         ax.grid(True, alpha=0.3)
#         ax.legend(title='Method/Layers', fontsize=12, loc='best')

#     plt.tight_layout()
#     plt.show()


def plot_matrix_components(all_matrix_components_by_layer, time, nt, layers):

    plt.style.use('seaborn-v0_8-talk')
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()
    time_axis = np.linspace(0, time, nt + 1)
    

    titles = [
        r'Diagonal Element $\rho_{00}$ (Occupation of $\ket{0}$)',
        r'Off-Diagonal Element $\text{Re}(\rho_{01})$',
        r'Off-Diagonal Element $\text{Re}(\rho_{10})$', 
        r'Diagonal Element $\rho_{11}$ (Occupation of $\ket{1}$)'
    ]
    # Keys for extraction: VQTE (v) and Exact (e)
    component_v_keys = ['va', 'vb', 'vc', 'vd']
    component_e_keys = ['ea', 'eb', 'ec', 'ed'] 


    num_layers = len(all_matrix_components_by_layer)
    vqte_colors = plt.cm.plasma(np.linspace(0.2, 0.9, num_layers))
    
    # Plotting loop
    for i, data in enumerate(all_matrix_components_by_layer):
  

        
        # Plot VQTE lines (solid, colored based on layer)
        for j in range(4):
            ax = axes[j]
            v_key = component_v_keys[j]
            # --- USE .get() with a fallback to avoid KeyError ---
            component_data = data.get(v_key, [])
            

            num_points = len(component_data)
            ax.plot(time_axis[:num_points], 
                    component_data, 
                    label=f'VQTE, {layers} Layers', 
                    color=vqte_colors[i],
                    linestyle='-',
                    linewidth=2,
                    alpha=0.7)

        # Plot Exact lines (dashed, black/dark gray, only once per subplot)
        if i == 0:
            for j in range(4):
                ax = axes[j]
                e_key = component_e_keys[j]
                component_data = data.get(e_key, []) 
                
                if not component_data:
                     print(f"Warning: No data found for Exact key '{e_key}'.")
                     continue

                num_points = len(component_data)
                ax.plot(time_axis[:num_points], 
                        component_data, 
                        label=f'Exact', 
                        color='black',
                        linestyle='--',
                        linewidth=3,
                        alpha=0.9)


    # Finalize plots
    for ax, title in zip(axes, titles):
        ax.set_title(title, fontsize=18, pad=10)
        ax.set_xlabel('Time (t)', fontsize=16)
        ax.set_ylabel('Real Value', fontsize=16)
        ax.grid(True, alpha=0.3)
        ax.legend(title='Method/Layers', fontsize=12, loc='best')

    plt.tight_layout()
    plt.show()


def plot_multiple_fidelity_vs_layers(results, time, nt):

    plt.style.use('seaborn-v0_8-talk')
    fig, ax = plt.subplots(figsize=(10, 6))
    time_axis = np.linspace(0, time, nt + 1)
    # Define a color palette
    colors = ['#003594', '#940000', '#009435', '#946C00', '#6C0094']
    
    # Check if all_results is a list of dictionaries (multiple scenarios)
    # or if it's just a single list of fidelity results (backward compatibility)
    for layer in range(len(results)):
            num_points = len(results[layer])
            plt.plot(time_axis[:num_points], 
                    results[layer], 
                    label=f'Exact Diag Site {layer+1} Occupation', 
                    marker='', 
                    linestyle='dashed')
    
    # Loop through all results and plot each line
    # for i, result in enumerate(all_results):
    #     label = result['label']
    #     layers_list = result['layers_list']
    #     fidelity_results = result['fidelity_results']
    #     color = colors[i % len(colors)]
        
    #     ax.plot(
    #         layers_list, 
    #         fidelity_results, 
    #         'o-', 
    #         linewidth=2, 
    #         markersize=8, 
    #         color=color,
    #         label=label
    #     )
    
    # Add titles, labels, and legend
    ax.set_title("Fidelity vs Ansatz Layers for Different Scenarios", fontsize=20, pad=20)
    ax.set_xlabel("Number of Ansatz Layers", fontsize=16)
    ax.set_ylabel("Fidelity", fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    
    # Only show legend if there are multiple scenarios

    # Final aesthetics
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.05)
    # ax.set_xlim(0.5, max(len(results)) )  # Adjust x-axis limits
    plt.tight_layout()
    plt.show()