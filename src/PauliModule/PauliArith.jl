import Base: *, +

###############################################################
# Pauli乘法
###############################################################
_phase_dict = Dict(0 => 1.0f0, 1 => 1.0f0im, 2 => -1.0f0, 3 => -1.0f0im)

function _mul_util(bvec1, bvec2)
    p1x = @view bvec1[1:end÷2]
    p1z = @view bvec1[end÷2+1:end]
    p2x = @view bvec2[1:end÷2]
    p2z = @view bvec2[end÷2+1:end]
    bvecz = p1z .⊻ p2z
    bvecx = p1x .⊻ p2x
    bvec = vcat(bvecx, bvecz)

    phase = (p1z'*p2x) % 2 == 0 ? 1.0f0 : -1.0f0
    phase1 = _phase_dict[(p1z'*p1x) % 4]
    phase2 = _phase_dict[(p2z'*p2x) % 4]
    phase3 = _phase_dict[(bvecz'*bvecx)*3 % 4]
    return bvec, phase * phase1 * phase2 * phase3
end

function Base.:*(pauli1::Pauli, pauli2::Pauli)
    bvec, phase = _mul_util(pauli1.bvec, pauli2.bvec)
    name = bvec2name(bvec)
    nqubits = length(pauli1.name)
	return PauliAlgebra(nqubits, Dict(name => phase), Dict(name => bvec))
end

function Base.:*(coeff::T, pauli::Pauli) where T <: Number
	name = pauli.name
	bvec = pauli.bvec
	nqubits = length(pauli.name)
	return PauliAlgebra(nqubits, Dict(name => coeff), Dict(name => bvec))
end


function Base.:*(pauli::Pauli, coeff::T) where T <: Number
	name = pauli.name
	bvec = pauli.bvec
	nqubits = length(pauli.name)
	return PauliAlgebra(nqubits, Dict(name => coeff), Dict(name => bvec))
end

function Base.:*(coeff::T, pauli::PauliAlgebra) where T <: Number
    nqubits = pauli.nqubits
    terms = Dict(name => coeff * val for (name, val) in pauli.terms)
    bvecs_dict = pauli.bvecs_dict

	return PauliAlgebra(nqubits, terms, bvecs_dict)
end


function Base.:*(pauli::PauliAlgebra, coeff::T) where T <: Number
    nqubits = pauli.nqubits
    terms = Dict(name => coeff * val for (name, val) in pauli.terms)
    bvecs_dict = pauli.bvecs_dict
    
	return PauliAlgebra(nqubits, terms, bvecs_dict)
end

function Base.:*(pauli1::PauliAlgebra, pauli2::Pauli)
    terms = Dict{String, ComplexF32}()
    bvecs = Dict{String, BitVector}()
    for (name, val) in pauli1.terms
        bvec, phase = _mul_util(pauli1.bvecs_dict[name], pauli2.bvec)
        phase = phase * pauli1.terms[name]
        
        name = bvec2name(bvec)
        terms[name] = val * phase
        bvecs[name] = bvec
    end
    return PauliAlgebra(pauli1.nqubits, terms, bvecs)
end

function Base.:*(pauli1::PauliAlgebra, pauli2::PauliAlgebra)
    terms = Dict{String, ComplexF32}()
    bvecs = Dict{String, BitVector}()
    for (name1, val1) in pauli1.terms, (name2, val2) in pauli2.terms
        bvec, phase = _mul_util(pauli1.bvecs_dict[name1], pauli2.bvecs_dict[name2])
        phase = phase * val1 * val2
        
            name = bvec2name(bvec)
            if haskey(terms, name)
                terms[name] += phase
            else
                terms[name] = phase
                bvecs[name] = bvec
            end
        
    end
    return PauliAlgebra(pauli1.nqubits, terms, bvecs)
end


###############################################################
# Pauli加法
###############################################################

function Base.:+(pauli1::PauliAlgebra, pauli2::PauliAlgebra)
    nqubits = pauli1.nqubits
	@assert nqubits == pauli2.nqubits "两个Pauli的qubit数量不一致"

    terms = pauli1.terms
    bvecs_dict = pauli1.bvecs_dict
    mergewith!(+, terms, pauli2.terms)
    merge!(bvecs_dict, pauli2.bvecs_dict)

	return PauliAlgebra(nqubits, terms, bvecs_dict)
end


function Base.:+(pauli1::Pauli, pauli2::PauliAlgebra)
	nqubits = length(pauli1.name)
	@assert nqubits == pauli2.nqubits "两个Pauli的qubit数量不一致"
	
	terms = copy(pauli2.terms)
	bvecs_dict = copy(pauli2.bvecs_dict)
	
	# 将Pauli添加到PauliAlgebra中
	if haskey(terms, pauli1.name)
		terms[pauli1.name] += ComplexF32(1.0)
	else
		terms[pauli1.name] = ComplexF32(1.0)
		bvecs_dict[pauli1.name] = pauli1.bvec
	end

	return PauliAlgebra(nqubits, terms, bvecs_dict)
end

function Base.:+(pauli1::PauliAlgebra, pauli2::Pauli)
	nqubits = pauli1.nqubits
	@assert nqubits == length(pauli2.name) "两个Pauli的qubit数量不一致"
	
	terms = copy(pauli1.terms)
	bvecs_dict = copy(pauli1.bvecs_dict)
	
	# 将Pauli添加到PauliAlgebra中
	if haskey(terms, pauli2.name)
		terms[pauli2.name] += ComplexF32(1.0)
	else
		terms[pauli2.name] = ComplexF32(1.0)
		bvecs_dict[pauli2.name] = pauli2.bvec
	end

	return PauliAlgebra(nqubits, terms, bvecs_dict)
end

function Base.:+(pauli1::Pauli, pauli2::Pauli)
    nqubits = length(pauli1.name)
    @assert nqubits == length(pauli2.name) "两个Pauli的qubit数量不一致"

	terms = Dict{String, ComplexF32}()
	bvecs_dict = Dict{String, BitVector}()
	
	# 将第一个Pauli添加到字典中
	terms[pauli1.name] = ComplexF32(1.0)
	bvecs_dict[pauli1.name] = pauli1.bvec
	
	# 将第二个Pauli添加到字典中
	if haskey(terms, pauli2.name)
		terms[pauli2.name] += ComplexF32(1.0)
	else
		terms[pauli2.name] = ComplexF32(1.0)
		bvecs_dict[pauli2.name] = pauli2.bvec
	end
	return PauliAlgebra(nqubits, terms, bvecs_dict)
end
