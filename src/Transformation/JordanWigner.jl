struct FermionicLadder
    annihilate::Vector{Int16}
    create::Vector{Int16}
end

function FermionicLadder(annihilate::Vector{T}, create::Vector{T}) where T <: Integer
    annihilate = convert(Vector{Int16}, annihilate)
    create = convert(Vector{Int16}, create)
    FermionicLadder(annihilate, create)
end

function _JWT(i::T, is_dagger::Bool) where T <: Integer
    name = 'Z'^(i-1)
    if is_dagger
        return 1/2 * Pauli(name * 'X') + 1im/2 * Pauli(name * 'Y')
    else
        return 1/2 * Pauli(name * 'X') + (- 1im/2) * Pauli(name * 'Y')
    end
end

function JordanWigner(ladder::FermionicLadder)
    a = ladder.annihilate
    a_dagger = ladder.create
    if !isempty(a)
        a_paulis = _JWT.(a, false)
    else
        a_paulis = 1.0
    end
    if !isempty(a_dagger)
        a_dagger_paulis = _JWT.(a_dagger, true)
    else
        a_dagger_paulis = 1.0
    end

    return prod(a_paulis) * prod(a_dagger_paulis)

end