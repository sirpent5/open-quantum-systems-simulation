
# from qiskit.quantum_info import SparsePauliOp
# import numpy as np
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


def perform_vqte(ham_real, ham_imag, init_state, mu, T, dt, nt, ansatz, init_param_values):
    real_var_principle = RealMcLachlanPrinciple(qgt=ReverseQGT(), gradient=ReverseEstimatorGradient(derivative_type=DerivativeType.IMAG))
    imag_var_principle = ImaginaryMcLachlanPrinciple(qgt=ReverseQGT(), gradient=ReverseEstimatorGradient())


# Initialize lists to store results
#second is the is expectation value of the number operator
    trace_list = [1.0]
    num_op_list = [np.trace(statevector_to_densitymatrix(init_state.data) @ np.array([[0, 0], [0, 1]])) / np.trace(statevector_to_densitymatrix(init_state.data))]
    print("Initial expectation value of number operator using VQE:", num_op_list[0])



    # Perform time evolution
    for t in range(nt):
        # Real evolution
        evolution_problem = TimeEvolutionProblem(ham_real, dt/2)
        var_qrte = VarQRTE(ansatz, init_param_values, real_var_principle, num_timesteps=1)
        evolution_result_re = var_qrte.evolve(evolution_problem)
        init_param_values = evolution_result_re.parameter_values[-1]
        
        # Imaginary evolution
        evolution_problem = TimeEvolutionProblem(ham_imag, dt/2)
        var_qite = VarQITE(ansatz, init_param_values, imag_var_principle, num_timesteps=1)
        evolution_result_im = var_qite.evolve(evolution_problem)
        init_param_values = evolution_result_im.parameter_values[-1]
        
        # Calculate the trace and expectation value of the number operator
        trace = np.trace(statevector_to_densitymatrix(Statevector(ansatz.assign_parameters(init_param_values)).data))
        trace_list.append(1.0) # Normalized so the trace is always 1
        num_op_list.append(np.trace(statevector_to_densitymatrix(Statevector(ansatz.assign_parameters(init_param_values)).data) @ np.array([[0, 0], [0, 1]])) / trace)
        #Stop
    return num_op_list, trace_list