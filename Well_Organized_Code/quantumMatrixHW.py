import numpy as np
import matplotlib.pyplot as plt

# Define identity and annihilation operator
I = np.eye(2)  # Identity matrix
a = np.array([[0, 1], [0, 0]])  # Annihilation operator
a_dag = a.T  # Creation operator

# Define annihilation and creation operators for two qubits
a_1 = np.kron(a, I)        # Acts on qubit 1
a_dag_1 = np.kron(a_dag, I)
a_2 = np.kron(I, a)        # Acts on qubit 2
a_dag_2 = np.kron(I, a_dag)

# Define hopping parameter
J = 1.0
# Construct Hamiltonian
H = J * (np.dot(a_dag_1, a_2) + np.dot(a_1, a_dag_2))

# Define basis states
basis_states = {
    "|00>": np.array([1, 0, 0, 0]),
    "|01>": np.array([0, 1, 0, 0]),
    "|10>": np.array([0, 0, 1, 0]),
    "|11>": np.array([0, 0, 0, 1])
}

# Time evolution parameters
# times = np.linspace(0, 2*np.pi, 100)  # Full oscillation period
times = np.linspace(0, 1, 100)  # Full oscillation period

def time_evolution(initial_state, H, times):
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    return [eigenvectors @ np.diag(np.exp(-1j * eigenvalues * t)) @ 
            np.linalg.inv(eigenvectors) @ initial_state for t in times]

# Evolve |01> state
evolved_states = time_evolution(basis_states["|01>"], H, times)

# Calculate probabilities
p_01 = [np.abs(np.vdot(basis_states["|01>"], state))**2 for state in evolved_states]
p_10 = [np.abs(np.vdot(basis_states["|10>"], state))**2 for state in evolved_states]

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(times, p_01, label='P(|01>)', linestyle='--')
plt.plot(times, p_10, label='P(|10>)', linestyle='-')
plt.xlabel('Time')
plt.ylabel('Probability')
plt.title('Hopping Dynamics')
plt.legend()
plt.grid(True)
plt.show()
