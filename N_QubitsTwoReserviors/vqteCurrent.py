import numpy as np
import matplotlib.pyplot as plt
from qiskit.quantum_info import SparsePauliOp, Statevector
from qiskit_algorithms import TimeEvolutionProblem, VarQRTE, VarQITE
from qiskit_algorithms.time_evolvers.variational import RealMcLachlanPrinciple, ImaginaryMcLachlanPrinciple
from qiskit.circuit.library import EfficientSU2
from qiskit_algorithms.gradients import ReverseQGT, ReverseEstimatorGradient
from qiskit.primitives import Estimator, Sampler

# --- Helper Functions ---

def fermi_dirac(energy, mu, T, k_B=1.0):
    """Calculates the Fermi-Dirac distribution."""
    if T < 1e-9: # Avoid division by zero at T=0
        return 1.0 if energy < mu else 0.0
    arg = (energy - mu) / (k_B * T)
    # Avoid overflow for large arguments in exp
    if arg > 500:
        return 0.0
    return 1.0 / (np.exp(arg) + 1.0)

def statevector_to_densitymatrix(v):
    """Converts a vectorized density matrix back to its matrix form."""
    m = int(np.sqrt(len(v)))
    return np.reshape(v, (m, m), order='F')

def liouvillian_generation(eps, gamma, mu, T):
    """
    Generates the Liouvillian for a single qubit coupled to a fermionic bath.
    """
    # System Hamiltonian
    H = 0.5 * eps * SparsePauliOp("Z")
    
    # Define jump operators
    sm = SparsePauliOp("X", coeffs=[0.5]) + SparsePauliOp("Y", coeffs=[0.5j])
    sp = SparsePauliOp("X", coeffs=[0.5]) - SparsePauliOp("Y", coeffs=[0.5j])

    # Thermal occupation from Fermi-Dirac distribution
    n_th = fermi_dirac(eps, mu, T)

    # Collapse operators
    c_ops = []
    # Decay
    if gamma * (1 - n_th) > 0:
        c_ops.append(np.sqrt(gamma * (1 - n_th)) * sm)
    # Excitation
    if gamma * n_th > 0:
        c_ops.append(np.sqrt(gamma * n_th) * sp)

    # Construct the Liouvillian using the qiskit.quantum_info.Lindblad class
    # The Liouvillian L is such that d(rho)/dt = L(rho)
    from qiskit.quantum_info.operators.channel import Lindblad
    lindblad = Lindblad(H, c_ops)
    L = lindblad.to_operator() # This is the Liouvillian superoperator

    # The VQTE simulates d(psi)/dt = -i * H_eff * psi, where psi = vec(rho)
    # So, we need H_eff = i * L
    H_eff = 1j * L

    # Decompose H_eff into real and imaginary parts (both must be Hermitian)
    H_re = 0.5 * (H_eff + H_eff.adjoint())
    H_im = -0.5j * (H_eff - H_eff.adjoint())

    return H_re, H_im


def perform_vqte(ham_real, ham_imag, ansatz, init_param_values, convergence_threshold=1e-4, max_steps=100):
    """
    Performs the VQTE simulation to find the steady state.
    """
    dt = 0.1 # Time step
    
    # Define variational principles
    real_var_principle = RealMcLachlanPrinciple(qgt=ReverseQGT(), gradient=ReverseEstimatorGradient())
    imag_var_principle = ImaginaryMcLachlanPrinciple(qgt=ReverseQGT(), gradient=ReverseEstimatorGradient())

    # The observable is the population of the excited state |1>, which is rho_11
    # In the vectorized picture, this corresponds to the operator |1><1| (x) I
    # For a 2-qubit system representing the density matrix, this is (I-Z)/2 (x) I
    num_op = SparsePauliOp("ZI", coeffs=[-0.5]) + SparsePauliOp("II", coeffs=[0.5])

    param_values = init_param_values
    last_exp_val = -1

    for t in range(max_steps):
        # Imaginary time evolution to find the ground state of H_re (part of finding the null space)
        evolution_problem_im = TimeEvolutionProblem(ham_imag, dt)
        var_qite = VarQITE(ansatz, param_values, imag_var_principle)
        evolution_result_im = var_qite.evolve(evolution_problem_im)
        param_values = evolution_result_im.parameter_values[-1]

        # Real time evolution
        evolution_problem_re = TimeEvolutionProblem(ham_real, dt)
        var_qrte = VarQRTE(ansatz, param_values, real_var_principle)
        evolution_result_re = var_qrte.evolve(evolution_problem_re)
        param_values = evolution_result_re.parameter_values[-1]
        
        # Calculate expectation value
        psi = Statevector(ansatz.assign_parameters(param_values))
        exp_val = psi.expectation_value(num_op).real
        
        # Check for convergence
        if abs(exp_val - last_exp_val) < convergence_threshold:
            print(f"Converged at step {t} with expectation value {exp_val:.4f}")
            return exp_val
        last_exp_val = exp_val

    print("Warning: VQTE did not converge within max steps.")
    return last_exp_val

# --- Simulation Parameters ---
eps = 1.0  # Energy splitting
num_qubits = 2 # 2 qubits to represent the 2x2 density matrix
ansatz = EfficientSU2(num_qubits, reps=1)
initial_parameters = np.full(ansatz.num_parameters, 2 * np.pi)
# --- 1. Plot vs. Gamma ---
gamma_list = np.linspace(0.1, 2.0, 10)
mu_fixed = 0.5
T_fixed = 0.1
pop_vs_gamma = []
for g in gamma_list:
    print(f"\n--- Running for Gamma = {g:.2f} ---")
    H_re, H_im = liouvillian_generation(eps, g, mu_fixed, T_fixed)
    ss_pop = perform_vqte(H_re, H_im, ansatz, initial_parameters)
    pop_vs_gamma.append(ss_pop)

# --- 2. Plot vs. Mu ---
mu_list = np.linspace(-2.0, 4.0, 10)
gamma_fixed = 0.5
T_fixed = 0.1
pop_vs_mu = []
for m in mu_list:
    print(f"\n--- Running for Mu = {m:.2f} ---")
    H_re, H_im = liouvillian_generation(eps, gamma_fixed, m, T_fixed)
    ss_pop = perform_vqte(H_re, H_im, ansatz, initial_parameters)
    pop_vs_mu.append(ss_pop)

# --- 3. Plot vs. Temperature ---
temp_list = np.linspace(0.1, 2.0, 10)
gamma_fixed = 0.5
mu_fixed = 0.5
pop_vs_temp = []
for T in temp_list:
    print(f"\n--- Running for Temp = {T:.2f} ---")
    H_re, H_im = liouvillian_generation(eps, gamma_fixed, mu_fixed, T)
    ss_pop = perform_vqte(H_re, H_im, ansatz, initial_parameters)
    pop_vs_temp.append(ss_pop)

# --- Plotting Results ---
plt.figure(figsize=(18, 5))
plt.subplot(1, 3, 1)
plt.plot(gamma_list, pop_vs_gamma, 'o-', label='VQTE')
plt.title(f"vs. Gamma (μ={mu_fixed}, T={T_fixed})")
plt.xlabel("Coupling Strength (γ)")
plt.ylabel("Excited State Population")
plt.grid(True)

plt.subplot(1, 3, 2)
plt.plot(mu_list, pop_vs_mu, 'o-', label='VQTE')
plt.title(f"vs. Chemical Potential (γ={gamma_fixed}, T={T_fixed})")
plt.xlabel("Chemical Potential (μ)")
plt.grid(True)

plt.subplot(1, 3, 3)
plt.plot(temp_list, pop_vs_temp, 'o-', label='VQTE')
plt.title(f"vs. Temperature (γ={gamma_fixed}, μ={mu_fixed})")
plt.xlabel("Temperature (T)")
plt.grid(True)

plt.tight_layout()
plt.legend()
plt.show()