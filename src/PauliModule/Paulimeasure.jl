function bitexpect(bits::T, pauli::String)::Int8 where T<:DitStr

    if any(c -> c != 'I' && c != 'Z', pauli)
        @warn "pauli串中包含非I、Z的元素，返回0"
        return 0
    end
    @assert length(pauli) == length(bits) "pauli串和bits长度不一致"

    parity = false
    @inbounds for i in eachindex(pauli)
        if pauli[i] == 'Z' && bits[i] == 1
            parity = !parity
        end
    end

    return ifelse(parity, Int8(-1), Int8(1))
end

bitexpect(bits::T, pauli::Pauli) where T<:DitStr = bitexpect(bits, pauli.name)

function bitexpect(bits::T, paulialg::PauliAlgebra) where T<:DitStr
    terms = paulialg.terms
    return sum(bitexpect(bits, pauli) * coeff for (pauli, coeff) in terms)
end


function expect_pauli_algebra(Hamiltonian::PauliAlgebra, state)
    hamil_yao = sum(Hamiltonian.terms) do (pauli, coeff)
        if abs(coeff.im) > 5e-7
            @warn "Aborting imaginary part of coefficient $coeff"
        end
        coeff = coeff.re
        to_yao(pauli) * coeff
    end
    return expect(hamil_yao, state)
end


