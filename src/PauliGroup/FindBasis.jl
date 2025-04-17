using LinearAlgebra
function mat(paulis::Vector{Pauli})
    return hcat([pauli.bvec for pauli in paulis]...)
end

"""
    column_echelon(M::Matrix{Bool}) -> Matrix{Bool}

对布尔矩阵进行列消元操作，返回简化后的矩阵。

# 参数
- `M`: 输入的布尔矩阵

# 返回
- 经过列消元后的布尔矩阵，其中非零列被保留

# 算法
使用高斯消元法，但使用异或(⊻)操作代替加减操作，因为是在GF(2)域上运算。
"""
function column_echelon(M::BitMatrix)::BitMatrix
    M_copy = copy(M)
    m, n = size(M_copy)
    pivot_cols = zeros(Int, m)  
    # print
    for col in 1:n
        pivot_row = findfirst(x -> x, @view M_copy[:, col])
        isnothing(pivot_row) && continue  
        if pivot_cols[pivot_row] != 0
            
            pivot_col = pivot_cols[pivot_row]
            M_copy[:, col] .⊻= @view M_copy[:, pivot_col]
        else
            pivot_cols[pivot_row] = col
            for next_col in (col+1):n
                M_copy[pivot_row, next_col] && (M_copy[:, next_col] .⊻= @view M_copy[:, col])
            end
        end
        
    end
    
    non_zero_cols = findall(j -> any(@view M_copy[:, j]), 1:n)
    return M_copy[:, non_zero_cols]
end



"""
    powerset(v::Vector{T}) where T -> Vector{Vector{T}}

生成一个向量的所有可能子集。

# 参数
- `v`: 输入向量

# 返回
- 包含所有可能子集的向量

# 示例
```julia
powerset([1, 2, 3])  # 返回 [[],  [1], [1,2], [1,2,3], [1,3], [2], [2,3], [3]]
```
"""
function powerset(v::Vector{T}) where T
    n = length(v)
    result = Vector{Vector{T}}()
    for i in 0:(2^n-1)
        subset = Vector{T}()
        for j in 1:n
            if (i & (1 << (j-1))) != 0
                push!(subset, v[j])
            end
        end
        push!(result, subset)
    end
    return result
end


function find_basis(M::BitMatrix)::BitMatrix
    M_copy = copy(M)
    m, n = size(M_copy)
    pivot_cols = zeros(Int, m)  
    
    # 列消元部分保持不变
    for col in 1:n
        pivot_row = findfirst(x -> x, @view M_copy[:, col])
        isnothing(pivot_row) && continue  
        if pivot_cols[pivot_row] != 0
            pivot_col = pivot_cols[pivot_row]
            M_copy[:, col] .⊻= @view M_copy[:, pivot_col]
        else
            pivot_cols[pivot_row] = col
            for next_col in (col+1):n
                M_copy[pivot_row, next_col] && (M_copy[:, next_col] .⊻= @view M_copy[:, col])
            end
        end
    end
    

    non_zero_cols = findall(j -> any(@view M_copy[:, j]), 1:n)
    M_copy = M_copy[:, non_zero_cols]
    
    # 获取矩阵的新维度
    n_rows, n_cols = size(M_copy)
    n = n_rows ÷ 2
    k = n - n_cols
    k <= 0 && return M_copy

    r = n_cols
    while k > 0
        # 找到所有pivot_cols为0的索引
        zero_indices = findall(x -> x == 0, pivot_cols)
        isempty(zero_indices) && break  # 如果没有零索引，退出循环
        first_zero_idx = zero_indices[1]
        # println(zero_indices[2:end])
        
        indd = powerset(zero_indices)
        sort!(indd, by = x -> length(x))
        indd = indd[2:end]
        # 生成所有可能的组合
        for combination in indd
            # 创建基础向量，长度为n_rows
            result_vector = falses(n_rows)
            # last_zero_idx必须为true
            # result_vector[first_zero_idx] = true


            # 根据当前组合设置其他零索引位置
            for idx in combination
                result_vector[idx] = true
            end

            if rank(column_echelon(hcat(M_copy, result_vector))) == r
                continue
            end
            # 检查result_vector是否与M_copy的所有列都commute
            is_all_commute = true
            for col in 1:size(M_copy, 2)
                if !is_commute(result_vector, M_copy[:, col])
                    is_all_commute = false
                    break
                end
                
            end
            
            # 如果都commute，将result_vector添加到M_copy中并终止循环
            if is_all_commute
                # println(result_vector)
                M_copy = hcat(M_copy, result_vector)
                # display(M_copy)
                # pivot_cols[first_zero_idx] = n-k+1
                break
            end
        end
        
        k -= 1
        r += 1
    end
    
    return M_copy
end

function find_basis(paulis::Vector{Pauli})
    M = find_basis(mat(paulis))
    return [Pauli(M[:, i]) for i in 1:size(M, 2)]
end