"""
    tens4checksymmetry(C4th) 

If the fourth-order tensor of material elasticity has the full set of
symmetries, return true; otherwise false.
"""
function tens4checksymmetry(C4th)
    for I = 1:3
        for J = 1:3
            for K = 1:3
                for L = 1:3
                    C4th[I, J, K, L] != C4th[K, L, I, J] && return false
                    C4th[I, J, K, L] != C4th[J, I, K, L] && return false
                    C4th[I, J, K, L] != C4th[I, J, L, K] && return false
                end
            end
        end
    end
    return true
end

"""
    tens4symmtto6x6t!(M::Matrix{T}, ST::Array{T, 4}) where {T}

Convert a symmetric 4th-order tensor to a 6 x 6 matrix.

!!! Note
The order corresponds to the arrangement of the components of stress (or
strain) tensor, symmetric, three-dimensional, into a 6-component 
vector.

# Example
```
J=tens4_ijkl(eye(3),eye(3))
produces the tracor:
T=rand(3); 
sum(diag(T))*eye(3)
t= tens4_dot_2(J,T)
M= tens4_symm_to_6(ST)
```
"""
function tens4symmtto6x6t!(M::Matrix{T}, ST::Array{T,4}) where {T}
    # This corresponds to the arrangement of the components of stress (or
    # strain) tensor, symmetric, three-dimensional, into a 6-component 
    # vector.
    ix = [1 1; 2 2; 3 3; 1 2; 1 3; 2 3]
    for j = 1:6
        for i = 1:6
            M[i, j] = ST[ix[i, 1], ix[i, 2], ix[j, 1], ix[j, 2]]
        end
    end
    return M
end

"""
    tens4symmt6x6tot!(ST::Array{T, 4}, M::Matrix{T}) where {T}

Convert a symmetric 6 x 6 matrix to a symmetric 4th-order tensor.

!!! Note
The order corresponds to the arrangement of the components of stress (or
strain) tensor, symmetric, three-dimensional, into a 6-component 
vector.
"""
function tens4symmt6x6tot!(ST::Array{T,4}, M::Matrix{T}) where {T}
    ix = [1 4 5; 4 2 6; 5 6 3]
    n = 3
    for i = 1:n
        for j = 1:n
            for k = 1:n
                for l = 1:n
                    ST[i, j, k, l] = M[ix[i, j], ix[k, l]]
                end
            end
        end
    end
    return ST
end

"""
    tens4dot2!(R::Array{T, 2}, F::Array{T, 4}, S::Array{T, 2}) where {T}

Compute the double contraction of a 4th-order and a 2nd-order tensors.

!!! note
The double contraction  of two second-order sensors is defined as 
`A:B = tr(A'*B) = A_ij B_ij`

The resulting second-order tensor is first zeroed out, and then the result is
accumulated.
"""
function tens4dot2!(R::Array{T,2}, F::Array{T,4}, S::Array{T,2}) where {T}
    R .= zero(T)
    for l = 1:3
        for k = 1:3
            for j = 1:3
                for i = 1:3
                    R[i, j] += F[i, j, k, l] * S[k, l]
                end
            end
        end
    end
    return R
end

"""
    tens4ijkl!(t::Array{T, 4}, A::FA, B::FB, op = :+) where {T, FA, FB}

Fill a 4th-order tensor as a dyadic product of two 2nd-order tensors.

The `i,j,k,l` component is given as `t[i,j,k,l]=A(i,j)*B(k,l)`.

!!! note
The tensor is accumulated to. It needs to be initialized to zero, if that is
desired as the initial state.

# Example
```
t = fill(0.0, 3, 3, 3, 3)
delta = (I, J) -> I == J ? 1.0 : 0.0
tens4ijkl!(t, delta, delta)
S = rand(3, 3)
@show tr(S) * I
tS = fill(0.0, 3, 3)
@show tens4dot2!(tS, t, S)
```
"""
function tens4ijkl!(t::Array{T,4}, A::FA, B::FB) where {T,FA,FB}
    for l = 1:3
        for k = 1:3
            for j = 1:3
                for i = 1:3
                    t[i, j, k, l] += A(i, j) * B(k, l)
                end
            end
        end
    end
    return t
end

"""
    tens4ikjl!(t::Array{T, 4}, A::FA, B::FB) where {T, FA, FB}

Fill a 4th-order tensor as a dyadic product of two 2nd-order tensors.

The `i,j,k,l` component is given as `t[i,j,k,l]=A(i,k)*B(j,l)`.

!!! note
The tensor is accumulated to. It needs to be initialized to zero, if that is
desired as the initial state.

# Example
```
t = fill(0.0, 3, 3, 3, 3)
delta = (I, J) -> I == J ? 1.0 : 0.0
tens4ikjl!(t, delta, delta)
S = rand(3, 3)
@show transpose(S) 
tS = fill(0.0, 3, 3)
@show transpose(S) - tens4dot2!(tS, t, S)
```
"""
function tens4ikjl!(t::Array{T,4}, A::FA, B::FB) where {T,FA,FB}
    for l = 1:3
        for k = 1:3
            for j = 1:3
                for i = 1:3
                    t[i, j, k, l] += A(i, k) * B(j, l)
                end
            end
        end
    end
    return t
end

"""
    tens4iljk!(t::Array{T, 4}, A::FA, B::FB) where {T, FA, FB}

Fill a 4th-order tensor as a dyadic product of two 2nd-order tensors.

The `i,j,k,l` component is given as `t[i,j,k,l]=A(i,l)*B(j,k)`.

!!! note
The tensor is accumulated to. It needs to be initialized to zero, if that is
desired as the initial state.

# Example
```
t = fill(0.0, 3, 3, 3, 3)
delta = (I, J) -> I == J ? 1.0 : 0.0
tens4iljk!(t, delta, delta)
S = rand(3, 3)
tS = fill(0.0, 3, 3)
@show S - tens4dot2!(tS, t, S)
```
"""
function tens4iljk!(t::Array{T,4}, A::FA, B::FB) where {T,FA,FB}
    for l = 1:3
        for k = 1:3
            for j = 1:3
                for i = 1:3
                    t[i, j, k, l] += A(i, l) * B(j, k)
                end
            end
        end
    end
    return t
end

"""
    tens4identity!(t::Array{T, 4}) where {T}

Compute 4th-order identity tensor.

# Example

The product of the identity tensor with the second-order tensor `S` is 
```
t = fill(0.0, 3, 3, 3, 3)
tens4identity!(t)
S = rand(3, 3)
tS = fill(0.0, 3, 3)
tens4dot2!(tS, t, S)
@show S - tS
```
"""
function tens4identity!(t::Array{T,4}) where {T}
    delta = (I, J) -> I == J ? 1.0 : 0.0
    t .= zero(T)
    return tens4ikjl!(t, delta, delta)
end

"""
    tens4transposor!(t::Array{T, 4}) where {T}

Compute 4th-order transposor tensor.

# Example

The product of the transposor tensor with the second-order tensor `S` is 
```
t = fill(0.0, 3, 3, 3, 3)
tens4transposor!(t)
S = rand(3, 3)
tS = fill(0.0, 3, 3)
tens4dot2!(tS, t, S)
@show S' - tS
```
"""
function tens4transposor!(t::Array{T,4}) where {T}
    delta = (I, J) -> I == J ? 1.0 : 0.0
    t .= zero(T)
    return tens4iljk!(t, delta, delta)
end

"""
    tens4tracor!(t::Array{T, 4}) where {T}

Compute 4th-order tracor tensor.

Double contraction of a second order tensor with this fourth-order tensor
produces the spherical part of the second order tensor.

# Example

The product of the tracor tensor with the second-order tensor `S` is 
```
t = fill(0.0, 3, 3, 3, 3)
tens4tracor!(t)
S = rand(3, 3)
tS = fill(0.0, 3, 3)
tens4dot2!(tS, t, S)
@show tr(S) * I - tS
```
"""
function tens4tracor!(t::Array{T,4}) where {T}
    delta = (I, J) -> I == J ? 1.0 : 0.0
    t .= zero(T)
    return tens4ijkl!(t, delta, delta)
end

"""
    tens4symmetrizor!(t::Array{T, 4}) where {T}

Compute 4th-order symmetrizor tensor.

Double contraction of a second order tensor with this fourth-order tensor
produces the symmetric part of the second order tensor.

# Example

The product of the symmetrizor tensor with the second-order tensor `S` is 
```
t = fill(0.0, 3, 3, 3, 3)
tens4symmetrizor!(t)
S = rand(3, 3)
tS = fill(0.0, 3, 3)
tens4dot2!(tS, t, S)
@show (S + S')/2 * I - tS
```
"""
function tens4symmetrizor!(t::Array{T,4}) where {T}
    delta = (I, J) -> I == J ? 1.0 : 0.0
    t .= zero(T)
    tens4ikjl!(t, delta, delta) # identity
    tens4iljk!(t, delta, delta) # transposor
    t .*= 0.5
    return t
end

"""
    tens4skewor!(t::Array{T, 4}) where {T}

Compute 4th-order skewor tensor.

Double contraction of a second order tensor with this fourth-order tensor
produces the skew part of the second order tensor.

# Example

The product of the skewor tensor with the second-order tensor `S` is 
```
t = fill(0.0, 3, 3, 3, 3)
tens4skewor!(t)
S = rand(3, 3)
tS = fill(0.0, 3, 3)
tens4dot2!(tS, t, S)
@show (S - S')/2 * I - tS
```
"""
function tens4skewor!(t::Array{T,4}) where {T}
    delta = (I, J) -> I == J ? 1.0 : 0.0
    t .= zero(T)
    tens4iljk!(t, delta, delta) # transposor
    t .= -t # subtract that part
    tens4ikjl!(t, delta, delta) # identity
    t .*= 0.5
    return t
end

"""
    tens4deviator!(t::Array{T, 4}) where {T}

Compute 4th-order deviator tensor.

Double contraction of a second order tensor with this fourth-order tensor
produces the deviator part of the second order tensor.

# Example

The product of the deviator tensor with the second-order tensor `S` is 
```
t = fill(0.0, 3, 3, 3, 3)
tens4deviator!(t)
S = rand(3, 3)
tS = fill(0.0, 3, 3)
tens4dot2!(tS, t, S)
@show tr((S - tr(S)/3*I) ), tr(tS)
```
"""
function tens4deviator!(t::Array{T,4}) where {T}
    delta = (I, J) -> I == J ? 1.0 : 0.0
    t .= zero(T)
    tens4tracor!(t) # tracor
    t .= -(1.0 / 3) .* t # subtract (1/3) that part
    tens4ikjl!(t, delta, delta) # identity
    t .*= 0.5
    return t
end
