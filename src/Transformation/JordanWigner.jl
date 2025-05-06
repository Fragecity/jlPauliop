struct FermionicLadder
    annihilate::Vector{Int16}
    create::Vector{Int16}
    nqubits::Int16
end

function FermionicLadder(annihilate::Vector{Int64}, create::Vector{Int64}, nq::Int64) 
    annihilate = convert(Vector{Int16}, annihilate)
    create = convert(Vector{Int16}, create)
    nq = convert(Int16, nq) 
    FermionicLadder(annihilate, create, nq)
end

function _JWT(i::T, is_dagger::Bool, nq::T) where T <: Integer
    name = 'Z'^(i-1)
    rest = 'I'^(nq - i)
    if is_dagger
        return 1/2 * Pauli(name * 'X'*rest) + 1im/2 * Pauli(name * 'Y'*rest)
    else
        return 1/2 * Pauli(name * 'X'*rest) + (- 1im/2) * Pauli(name * 'Y'*rest)
    end
end

function JordanWigner(ladder::FermionicLadder)
    a = ladder.annihilate
    a_dagger = ladder.create
    nq = ladder.nqubits
    if !isempty(a)
        a_paulis = _JWT.(a, false, nq)
    else
        a_paulis = 1.0
    end
    if !isempty(a_dagger)
        a_dagger_paulis = _JWT.(a_dagger, true, nq)
    else
        a_dagger_paulis = 1.0
    end
    return prod(a_paulis) * prod(a_dagger_paulis)
end