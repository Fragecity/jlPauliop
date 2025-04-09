pauli_dict = Dict('I' => [false, false], 'X' => [true, false], 'Y' => [true, true], 'Z' => [false, true])
pauli_dict_rev = Dict([false, false] => 'I', [true, false] => 'X', [true, true] => 'Y', [false, true] => 'Z')
pauli_yao_dict = Dict('I' => I2, 'X' => X, 'Y' => Y, 'Z' => Z)


struct Pauli
	name::String
	bvec::BitVector
end

function Base.show(io::IO, p::Pauli)
    println(io, p.name)
end


function name2bvec(name::String)
	bvec_x = BitVector([pauli_dict[char][1] for char in name])
	bvec_z = BitVector([pauli_dict[char][2] for char in name])
	return vcat(bvec_x, bvec_z)
end

function bvec2name(bvec::BitVector)
	n = length(bvec) ÷ 2
	bvec_x = bvec[1:n]
	bvec_z = bvec[n+1:end]
	return join([pauli_dict_rev[[bvec_x[i], bvec_z[i]]] for i in 1:n])
end


function Pauli(name::String)
	Pauli(name, name2bvec(name))
end

function Pauli(bvec::AbstractArray)
	@assert all(x -> x == 0 || x == 1, bvec) "bvec中的元素必须为0或1"
	@assert length(bvec) % 2 == 0 "bvec的长度必须为偶数"
	bvec_bit = BitVector(bvec .== 1)
	name = bvec2name(bvec_bit)
	Pauli(name, bvec_bit)
end

function Pauli(bvec::BitVector)
	Pauli(bvec2name(bvec), bvec)
end

function is_commute(pauli1::Pauli, pauli2::Pauli)
	p1x, p1z = pauli1.bvec[1:end÷2], pauli1.bvec[end÷2+1:end]
	p2x, p2z = pauli2.bvec[1:end÷2], pauli2.bvec[end÷2+1:end]
	return !(reduce(⊻, p1x .* p2z) ⊻ reduce(⊻, p1z .* p2x))
end

function is_commute(pauli1::BitVector, pauli2::BitVector)
	p1x, p1z = pauli1[1:end÷2], pauli1[end÷2+1:end]
	p2x, p2z = pauli2[1:end÷2], pauli2[end÷2+1:end]
	return !(reduce(⊻, p1x .* p2z) ⊻ reduce(⊻, p1z .* p2x))
end

function is_commute(pauli1::String, pauli2::String)
	pauli1 = Pauli(pauli1)
	pauli2 = Pauli(pauli2)
	return is_commute(pauli1, pauli2)
end

to_yao(pauli::Pauli) = kron((pauli_yao_dict[p.name] for p in pauli.name)...)
to_yao(pauli::String) = kron((pauli_yao_dict[p] for p in pauli)...)
