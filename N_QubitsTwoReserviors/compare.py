import numpy as np

def verify_algebra(superoperator, ham_real, ham_imag):
    """
    Verifies the algebraic properties of the split Hamiltonians against the 
    exact diagonalization superoperator.
    
    """
    print("\n Algebra Check")

    # 1. Convert inputs to numpy matrices

    h_real_mat = ham_real.to_matrix()
    h_imag_mat = ham_imag.to_matrix()
    
    # 2. Construct h_end
    h_end = h_real_mat + 1j * h_imag_mat

    # 3. Run Assertions
    all_passed = True

    # Check 1: h_end == superoperator
    # Corresponds to: @assert h_end == superoperator
    if np.allclose(h_end, superoperator):
        print("✓ h_end matches Superoperator (superoperator)")
    else:
        print(f"X h_end mismatch! Norm diff: {np.linalg.norm(h_end - superoperator):.2e}")
        all_passed = False

    # Check 2: h_real is Hermitian
    # Corresponds to: @assert h_real == adjoint(h_real)
    if np.allclose(h_real_mat, h_real_mat.conj().T):
        print("✓ h_real is Hermitian")
    else:
        print("X h_real is NOT Hermitian")
        all_passed = False

    # Check 3: i*h_imag is Anti-Hermitian
    # Corresponds to: @assert 1im*h_imag == -adjoint(1im*h_imag)
    lhs = 1j * h_imag_mat
    rhs = -(1j * h_imag_mat).conj().T
    if np.allclose(lhs, rhs):
        print("✓ i*h_imag is Anti-Hermitian")
    else:
        print("X i*h_imag is NOT Anti-Hermitian")
        all_passed = False

    # Check 4: Real extraction logic
    # Corresponds to: @assert h_real == 0.5*(h_end + adjoint(h_end))
    if np.allclose(h_real_mat, 0.5 * (h_end + h_end.conj().T)):
        print("✓ h_real extraction correct")
    else:
        print("X h_real extraction failed")
        all_passed = False

    # Check 5: Imag extraction logic
    # Corresponds to: @assert h_imag == 0.5*(h_end - adjoint(h_end))/1im
    if np.allclose(h_imag_mat, 0.5 * (h_end - h_end.conj().T) / 1j):
        print("✓ h_imag extraction correct")
    else:
        print("X h_imag extraction failed")
        all_passed = False

    print("------------------------------\n")
    return all_passed

def diagnose_mismatch(superoperator, ham_real, ham_imag):
    H_end = ham_real.to_matrix() + 1j * ham_imag.to_matrix()
    diff = H_end - superoperator
    
    print("\n--- Diagnostic ---")
    print(f"Max difference: {np.max(np.abs(diff))}")
    
    # Check the Diagonal (Population decay terms)
    diag_diff = np.diag(diff)
    print(f"Diagonal Difference Norm: {np.linalg.norm(diag_diff)}")
    
    # Check the Off-Diagonal (Coherence/Hamiltonian terms)
    off_diag_diff = diff - np.diag(diag_diff)
    print(f"Off-Diagonal Difference Norm: {np.linalg.norm(off_diag_diff)}")
    
    if np.linalg.norm(diag_diff) < 1e-5:
        print(">> Diagnosis: Your DISSIPATOR (decay) terms match, but your HAMILTONIAN (energy) terms are wrong.")
    elif np.linalg.norm(off_diag_diff) < 1e-5:
        print(">> Diagnosis: Your HAMILTONIAN terms match, but your DISSIPATOR terms are wrong.")
    else:
        print(">> Diagnosis: Both parts are different. Likely a Basis/Vectorization mismatch (e.g., A⊗B vs B⊗A).")