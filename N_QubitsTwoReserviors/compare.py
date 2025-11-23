import numpy as np

def verify_algebra(superoperator, ham_real, ham_imag):

    print("\n Algebra Check")

    h_real_mat = ham_real.to_matrix()
    h_imag_mat = ham_imag.to_matrix()

    h_end = h_real_mat + 1j * h_imag_mat

    
    if np.allclose(h_end, superoperator):
        print("h_end matches Superoperator")
    else:
        print(f"h_end mismatch Norm diff: {np.linalg.norm(h_end - superoperator)}")
        all_passed = False

    if np.allclose(h_real_mat, h_real_mat.conj().T):
        print("h_real is Hermitian")
    else:
        print("h_real is NOT Hermitian")
 

    lhs = 1j * h_imag_mat
    rhs = -(1j * h_imag_mat).conj().T
    if np.allclose(lhs, rhs):
        print("i*h_imag is Anti-Hermitian")
    else:
        print(" i*h_imag is not Anti-Hermitian")
  

    if np.allclose(h_real_mat, 0.5 * (h_end + h_end.conj().T)):
        print("h_real extraction correct")
    else:
        print("h_real extraction failed")
    

    if np.allclose(h_imag_mat, 0.5 * (h_end - h_end.conj().T) / 1j):
        print(" h_imag extraction correct")
    else:
        print("h_imag extraction failed")
  

    return 

def diagnose_mismatch(superoperator, ham_real, ham_imag):
    
    H_end = ham_real.to_matrix() + 1j * ham_imag.to_matrix()
    diff = H_end - superoperator
    
    print(f"Max difference: {np.max(np.abs(diff))}")
    
    # Check diag
    diag_diff = np.diag(diff)
    print(f"Diagonal Difference Norm: {np.linalg.norm(diag_diff)}")
    
    # Check the coherant
    off_diag_diff = diff - np.diag(diag_diff)
    print(f"Off-Diagonal Difference Norm: {np.linalg.norm(off_diag_diff)}")
    
    if np.linalg.norm(diag_diff) < 1e-5:
        print(">> Diagnosis: Your DISSIPATOR (decay) terms match, but your HAMILTONIAN (energy) terms are wrong.")
    elif np.linalg.norm(off_diag_diff) < 1e-5:
        print(">> Diagnosis: Your HAMILTONIAN terms match, but your DISSIPATOR terms are wrong.")
    else:
        print(">> Diagnosis: Both parts are different. Likely a Basis/Vectorization mismatch (e.g., A⊗B vs B⊗A).")