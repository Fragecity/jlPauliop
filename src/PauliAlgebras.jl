module PauliAlgebras

# using Yao
import Yao: mat, kron, I2, X, Y, Z, expect
import Yao: DitStr
using LinearAlgebra: tr
using Base.Iterators: product, repeated


export 
# Pauli operations
mat, bitexpect, expect_pauli_algebra,
# Pauli and its operations
Pauli, is_commute, to_yao, 
# PauliAlgebra 
PauliAlgebra, rand_pauli_algebra, tomography,
# PauliGroup and its operations
find_basis, column_echelon,
# Pauli Linear Algebra
norm2,
# Jordan Wigner Transformation
JordanWigner,FermionicLadder

include("PauliModule/Pauli.jl")
include("PauliModule/PauliAlg.jl")
include("PauliModule/Paulimeasure.jl")
include("PauliModule/PauliArith.jl")
include("PauliModule/PauliLA.jl")
include("PauliGroup/FindBasis.jl")
include("Transformation/JordanWigner.jl")

end # module PauliAlgebras

