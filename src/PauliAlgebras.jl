module PauliAlgebras

# using Yao
import Yao: mat, kron, I2, X, Y, Z, expect, sandwich, density_matrix, state, subblocks, nqubits
import Yao: DitStr
using LinearAlgebra: tr
using Base.Iterators: product, repeated


export 
# Pauli operations
mat, bitexpect, expect_pauli_algebra, expect_pauli_algebra_trace,
# Pauli and its operations
Pauli, is_commute, to_yao, 
# PauliAlgebra 
PauliAlgebra, rand_pauli_algebra, tomography,
# PauliGroup and its operations
find_basis, column_echelon,
# Pauli Linear Algebra
norm2,
# Jordan Wigner Transformation
JordanWigner,FermionicLadder,
# PauliAlgebra and its operations
get_depth

include("ExtendYao.jl")
include("PauliModule/Pauli.jl")
include("PauliModule/PauliAlg.jl")
include("PauliModule/Paulimeasure.jl")
include("PauliModule/PauliArith.jl")
include("PauliModule/PauliLA.jl")
include("PauliGroup/FindBasis.jl")
include("Transformation/JordanWigner.jl")

end # module PauliAlgebras

