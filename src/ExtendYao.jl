function _update_gatecount!(gates_count, gate_range)
    for i in gate_range
        gates_count[i] += 1
    end
end

function _update_gatecount_byblock!(gates_count, sb)
    if hasfield(typeof(sb), :ctrl_locs)
        _update_gatecount!(gates_count, sb.ctrl_locs)
    end
    if hasfield(typeof(sb), :locs)
        _update_gatecount!(gates_count, sb.locs)
    end
    if hasfield(typeof(sb), :blocks)
        for ssb in sb.blocks
            _update_gatecount_byblock!(gates_count, ssb)
        end
    end
end

function get_depth(circuit)
    n = nqubits(circuit)
    gates_count = zeros(Int, n)
    for sb in subblocks(circuit)
        _update_gatecount_byblock!(gates_count, sb)
    end
    return maximum(gates_count)
end
