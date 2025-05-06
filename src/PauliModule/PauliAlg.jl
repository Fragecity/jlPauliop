import Base: show

struct PauliAlgebra
    nqubits::Int
    terms::Dict{String, ComplexF32}
    bvecs_dict::Dict{String, BitVector}
end

function Base.show(io::IO, p::PauliAlgebra)
    println(io, "\nPauliAlgebra($(p.nqubits) qubits):")
    terms_array = [(pauli, coeff) for (pauli, coeff) in p.terms if !(real(coeff) == 0 && imag(coeff) == 0)]
    for i in 1:5:length(terms_array)
        chunk = terms_array[i:min(i+4, length(terms_array))]
        for (pauli, coeff) in chunk   
            # 使用内置函数格式化复数
            formatted_coeff = if iszero(imag(coeff))
                string(round(real(coeff), digits=3))
            elseif iszero(real(coeff))
                string(round(imag(coeff), digits=3), "im")
            else
                string(round(real(coeff), digits=3), " + ", round(imag(coeff), digits=3), "im")
            end
            # 替换掉可能出现的"-0.0"为"0.0"
            formatted_coeff = replace(formatted_coeff, "-0.0" => "0.0")
            # 使用格式化后的系数进行输出
            print(io, "\t$(formatted_coeff) * $pauli")
        end
        println(io)
    end
end


function PauliAlgebra(paulis::Vector{Pauli}, coeffs::Vector{Number})
    nqubits = length(paulis[1].name)
    coeffs = ComplexF32.(coeffs)
    terms = Dict{String, ComplexF32}(zip((pauli.name for pauli in paulis), coeffs))
    bvecs_dict = Dict{String, BitVector}(zip((pauli.name for pauli in paulis), (pauli.bvec for pauli in paulis)))
    return PauliAlgebra(nqubits, terms, bvecs_dict)
end

function PauliAlgebra(terms::Dict{String, T}) where T<:Number
    nqubits = length(first(keys(terms)))
    terms_names = keys(terms)
    terms_coeffs = ComplexF32.(values(terms))
    terms_bvecs = name2bvec.(terms_names)
    terms = Dict{String, ComplexF32}(zip(terms_names, terms_coeffs))
    bvecs_dict = Dict{String, BitVector}(zip(terms_names, terms_bvecs))
    return PauliAlgebra(nqubits, terms, bvecs_dict)
end


function tomography(M::AbstractMatrix{<:Number})
    @assert size(M, 1) == size(M, 2) "矩阵必须是方阵"
    nqubits = size(M, 1) |> log2 |> Int
    single_paulis = (I2, X, Y, Z)
    pauli_dict = Dict(I2 => 'I', X => 'X', Y => 'Y', Z => 'Z')

    terms = Dict{String, Float32}()
    for combo in Iterators.product(Iterators.repeated(single_paulis, nqubits)...)
        pauli = kron(combo...)
        pstr = join(pauli_dict[pauli] for pauli in combo)
        coeff = tr(M * mat(pauli)).re / 2^nqubits |> ComplexF32
        if !isapprox(coeff.re, 0f0, atol=1e-9)
            terms[pstr] = coeff
        end
    end

    return PauliAlgebra(terms)
end


# function rand_pauli_algebra(nqubits::Int, real = false)
#     type = real ? Float32 : ComplexF32
#     coeff = rand(type, 4^nqubits-1)
#     terms_itr = product(repeated(['I', 'X', 'Y', 'Z'], nqubits)...) 
#     terms = (join(t) for t in terms_itr if join(t) != repeat('I', nqubits))
#     return PauliAlgebra(Dict{String, ComplexF32}(zip(terms, coeff)))
# end



function mat(pauli_algebra::PauliAlgebra)
    return sum(to_yao(pauli) * coeff for (pauli, coeff) in pauli_algebra.terms) |> mat
end


function shrink(palg::PauliAlgebra)
    nqubits = palg.nqubits
    terms = palg.terms
    
    for (name, val) in terms
        if isapprox(val, 0f0, atol=5e-10)
            delete!(terms, name)
        end
    end
    bvecs_dict = Dict{String, BitVector}(term => bvec for (term, bvec) in palg.bvecs_dict)
    
    return PauliAlgebra(nqubits, terms, bvecs_dict)
end
