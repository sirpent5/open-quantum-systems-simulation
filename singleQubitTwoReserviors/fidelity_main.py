import os
from imports import *
from exactDiagMethods import perform_exact_diag, build_exact_diag_hamiltonian
from globalMethods import build_initial_states
from vqteMethods import  hamiltonian_generation, perform_vqte
from fidelity_calculate import calculate_fidelity, plot_multiple_fidelity_vs_layers, plot_matrix_components,extract_density_matrix_components


def run_multiple_layers(maxLayers):
    params = {
        'gamma_L': 2,
        'gamma_R': 3,
        'eps': 1.0,
        'mu_L': 1,
        'mu_R': 2,
        'T_L': 1,
        'T_R': 1,
        'time': 4,
        'dt': 0.1,
    }
    nt = int(params['time']/params['dt'])
    
    
    layers_list = list(range(1, maxLayers + 1))
    fidelity_results = []
    all_fidelities_over_time = []
    all_matrix_components = []

    
    os.makedirs('fidelity_results', exist_ok=True)
    
    for layers in layers_list:
        print(f"Running simulation with {layers} layers...")
        
        eps = params['eps']
        nt = int(params['time'] / params['dt'])
        
        beta_L = 1 / params['T_L']
        beta_R = 1 / params['T_R']
        F_L = 1 / (np.exp(beta_L * (eps - params['mu_L'])) + 1)
        F_R = 1 / (np.exp(beta_R * (eps - params['mu_R'])) + 1)
        
        exact_diag_ham = build_exact_diag_hamiltonian(eps)
        ham_real, ham_imag = hamiltonian_generation(eps, params['gamma_L'], params['gamma_R'], F_R, F_L)
        
        ansatz = EfficientSU2(ham_real.num_qubits, reps = layers)
        print('Number of layers', layers)
        # Build initial states with specified number of layers
        vqte_init_state, exact_diag_init_state, ansatz, init_param_values = build_initial_states(ham_real, ansatz)
        
        # Run simulations
        exact_diag_results, time_points, exact_fidelity = perform_exact_diag(
            params['gamma_L'], F_L, params['gamma_R'], F_R, 
            params['dt'], nt, exact_diag_init_state, exact_diag_ham
        )
        
        vqte_results, vqte_fidelity = perform_vqte(
            ham_real, ham_imag, vqte_init_state, 
            params['dt'], nt, ansatz, init_param_values
        )
        
        components_over_time = extract_density_matrix_components(vqte_fidelity, exact_fidelity)
        all_matrix_components.append(components_over_time)
        
        
    # Calculate fidelity over time between VQTE and exact results
        #fidelity_over_time = calculate_fidelity(vqte_fidelity, exact_fidelity)
        #all_fidelities_over_time.append(fidelity_over_time)
        

    # Store the final fidelity value for plotting vs layers

    plot_matrix_components(all_matrix_components, params['time'],nt, maxLayers)
    # plot_data = [{
    # 'label': 'Default Scenario',
    # 'layers_list': layers_list,
    # 'fidelity_results': all_fidelities_over_time
    # # 'time' : time_points
    # }]
    #plot_multiple_fidelity_vs_layers(all_fidelities_over_time, params['time'], nt)
    #plot_matrix_components(all_matrix_components, params['time'],nt)
    
    #print(all_fidelities_over_time[-1][-1])
    
    return layers_list, fidelity_results

# Run the simulation
# layers, fidelities = run_multiple_layers()