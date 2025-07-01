##Define Pauli Matrices, Anihilation/ Creation and Number Operator
import numpy as np
Sigma_x = np.matrix([[0, 1], [1, 0]])
Sigma_y = np.matrix([[0, -1j], [1j, 0]])
Sigma_z = np.matrix([[1, 0], [0, -1]])
Sigma_plus = (Sigma_x +1j*Sigma_y)/2
Sigma_minus = Sigma_plus.getH()
numberop = Sigma_minus@Sigma_plus

##Parameters