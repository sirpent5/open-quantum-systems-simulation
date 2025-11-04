import os
from imports import *
from exactDiagMethods import *
from globalMethods import *
from vqteMethods import  *
from fidelity_calculate import *
    
# layers_fidelity = []

    # build hamiltonian(layer)
    
    # perform vqte
        # for timestep
            #store density matrix ---> vqte_matrix
        # return vqte_matrix
    # perform exact
        #for timestep
            #store matrix ---> exact matrix
    
    #calculate fidelity (exact_matrix, vqte_matrix)
        #for timestep
            #calculate fidelity ----> fidelity_list
        #return fidelity list
    #layers_fidelity.append(fidelity_list)
#plot layers_fidelity


def run_multiple_layers(maxLayers,time, dt):
    params = {
        'gamma_L': 2.0,
        'gamma_R': 2.0,
        'eps': 1.0,
        'mu_L': 1.5,
        'mu_R': 3,
        
        'T_L': 1,
        'T_R': 1,
        'time': time,
        'dt': dt,
    }
    nt = int(time/dt)
    eps = params['eps']
    nt = int(params['time'] / params['dt'])
        
    beta_L = 1 / params['T_L']
    beta_R = 1 / params['T_R']
    F_L = 1 / (np.exp(beta_L * (eps - params['mu_L'])) + 1)
    F_R = 1 / (np.exp(beta_R * (eps - params['mu_R'])) + 1)
        

    layers_list = list(range(1, maxLayers + 1))
    layers_fidelity = []
    all_matrix_components = []
 

    for layers in layers_list:
        print(f"Running simulation with {layers} layers...")
        

        exact_diag_ham = build_exact_diag_hamiltonian(eps)
        ham_real, ham_imag = hamiltonian_generation(eps, params['gamma_L'], params['gamma_R'], F_R, F_L)
        

        print('Number of layers', layers)
        # Build initial states with specified number of layers
        vqte_init_state, exact_diag_init_state, init_param_values, ansatz = build_initial_states(ham_real, layers)
        
        # Run simulations
        exact_results, time_points,exact_fidelity = perform_exact_diag(
            params['gamma_L'], F_L, params['gamma_R'], F_R, 
            params['dt'], nt, exact_diag_init_state, exact_diag_ham
        )
        
        vqte_results, vqte_fidelity = perform_vqte(
            ham_real, ham_imag, vqte_init_state, 
            params['dt'], nt, ansatz, init_param_values
        )
        
        
        fidelities = calculate_fidelity(vqte_fidelity, exact_fidelity)
        
        layers_fidelity.append(fidelities)
        
        components_over_time = extract_density_matrix_components(vqte_fidelity, exact_fidelity)
        all_matrix_components.append(components_over_time)

    plot_matrix_components(all_matrix_components, time, nt, layers_list)

    plot_multiple_fidelity_vs_layers(layers_fidelity, time, nt)


    # return layers_list, fidelity_results




# def run_multiple_layers(maxLayers,time, dt):
#     params = {
#         'gamma_L': 2.0,
#         'gamma_R': 2.0,
#         'eps': 1.0,
#         'mu_L': 1.5,
#         'mu_R': 3,
        
#         'T_L': 1,
#         'T_R': 1,
#         'time': time,
#         'dt': dt,
#     }
#     nt = int(time/dt)
    

#     layers_list = list(range(1, maxLayers + 1))
#     fidelity_results = []
#     all_fidelities_over_time = []
#     all_matrix_components = []

    
#     # os.makedirs('fidelity_results', exist_ok=True)
    
#     for layers in layers_list:
#         print(f"Running simulation with {layers} layers...")
        
#         eps = params['eps']
#         nt = int(params['time'] / params['dt'])
        
#         beta_L = 1 / params['T_L']
#         beta_R = 1 / params['T_R']
#         F_L = 1 / (np.exp(beta_L * (eps - params['mu_L'])) + 1)
#         F_R = 1 / (np.exp(beta_R * (eps - params['mu_R'])) + 1)
        
#         exact_diag_ham = build_exact_diag_hamiltonian(eps)
#         ham_real, ham_imag = hamiltonian_generation(eps, params['gamma_L'], params['gamma_R'], F_R, F_L)
        
#         ansatz = EfficientSU2(ham_real.num_qubits, reps = layers)
#         print('Number of layers', layers)
#         # Build initial states with specified number of layers
#         vqte_init_state, exact_diag_init_state, ansatz, init_param_values = build_initial_states(ham_real, ansatz)
        
#         # Run simulations
#         exact_results, time_points,exact_fidelity = perform_exact_diag(
#             params['gamma_L'], F_L, params['gamma_R'], F_R, 
#             params['dt'], nt, exact_diag_init_state, exact_diag_ham
#         )
        
#         vqte_results, vqte_fidelity = perform_vqte(
#             ham_real, ham_imag, vqte_init_state, 
#             params['dt'], nt, ansatz, init_param_values
#         )
        
        
#         fidelities = calculate_fidelity(vqte_fidelity, exact_fidelity)
        
#         all_fidelities_over_time.append(fidelities)
        
#         print(f"Fidelities for {layers} layers:", fidelities)
        
#         components_over_time = extract_density_matrix_components(vqte_fidelity, exact_fidelity)
#         all_matrix_components.append(components_over_time)

#     plot_matrix_components(all_matrix_components, time, nt, layers_list)
#     plot_multiple_fidelity_vs_layers(all_fidelities_over_time, time, nt)

#     output_results(vqte_results, 
#     exact_results, 
#     time_points, 
#     params['mu_L'], params['T_L'], params['gamma_L'], 
#     params['mu_R'], params['T_R'], params['gamma_R']
#     )
#     return layers_list, fidelity_results
