function norm2(pauli_algebra::PauliAlgebra)
    coeffs = values(pauli_algebra.terms)
    return sqrt(sum(abs2, coeffs))
end

