################################################################################
# 3-D model

"""
    dett(::Type{DeforModelRed3D},  C::Matrix{T}) where {T}

Compute the determinant of a general square matrix.
"""
function dett(::Type{DeforModelRed3D},  C::Matrix{T}) where {T}
	return (C[1,1] * C[2,2] * C[3,3] + 
		C[1,2] * C[2,3] * C[3,1] + 
		C[1,3] * C[2,1] * C[3,2] - 
		C[1,3] * C[2,2] * C[3,1] - 
		C[1,2] * C[2,1] * C[3,3] - 
		C[1,1] * C[2,3] * C[3,2])
end

"""
    strainvdet(::Type{DeforModelRed3D},  Cv::Vector{T}) where {T}

Compute the determinant of a symmetric strain-like square matrix represented
as a vector. Remember that the shear strain components are twice the entries
of the matrix representation.
"""
function strainvdet(::Type{DeforModelRed3D},  Cv::Vector{T}) where {T}
	return (Cv[1] * Cv[2] * Cv[3] + 
		Cv[4]/2 * Cv[6]/2 * Cv[5]/2 + 
		Cv[5]/2 * Cv[4]/2 * Cv[6]/2 -  
		Cv[5]/2 * Cv[2] * Cv[5]/2 - 
		Cv[4]/2 * Cv[4]/2 * Cv[3] - 
		Cv[1] * Cv[6]/2 * Cv[6]/2)
end

"""
    strainvtr(::Type{DeforModelRed3D},  Cv::Vector{T}) where {T}

Compute the trace of a symmetric strain-like square matrix represented as a
vector.
"""
function strainvtr(::Type{DeforModelRed3D},  Cv::Vector{T}) where {T}
	return (Cv[1] + Cv[2] + Cv[3])
end

"""
    strainttov!(::Type{DeforModelRed3D}, v::Vector{T}, t::Matrix{T}) where {T}

Convert a symmetric matrix of 3x3 strain components  into a 6-component vector.
"""
function strainttov!(::Type{DeforModelRed3D}, v::Vector{T}, t::Matrix{T}) where {T}
	v[1] = t[1,1];
	v[2] = t[2,2];
	v[3] = t[3,3];
	v[4] = t[1,2] + t[2,1];
	v[5] = t[1,3] + t[3,1];
	v[6] = t[3,2] + t[2,3];
	return v
end

"""
    strainvtot!(::Type{DeforModelRed3D}, t::Matrix{T}, v::Vector{T}) where {T}

Convert a strain 3-vector to a  matrix of 2x2 strain components (symmetric tensor).
"""
function strainvtot!(::Type{DeforModelRed3D}, t::Matrix{T}, v::Vector{T}) where {T}
	t[1,1] = v[1];
	t[2,2] = v[2];
	t[3,3] = v[3];
	t[1,2] = v[4]/2.;
	t[2,1] = v[4]/2.;
	t[1,3] = v[5]/2.;
	t[3,1] = v[5]/2.;
	t[3,2] = v[6]/2.;
	t[2,3] = v[6]/2.;
    return t
end

"""
    stressttov!(::Type{DeforModelRed3D}, v::Vector{T}, t::Matrix{T}) where {T}

Convert a symmetric matrix of 3x3 stress components to a 6-component vector.
"""
function stressttov!(::Type{DeforModelRed3D}, v::Vector{T}, t::Matrix{T}) where {T}
	v[1] = t[1,1];
	v[2] = t[2,2];
	v[3] = t[3,3];
	v[4] = 1/2.0*(t[1,2] + t[2,1]);
	v[5] = 1/2.0*(t[1,3] + t[3,1]);
	v[6] = 1/2.0*(t[3,2] + t[2,3]);
	return v
end

"""
    stressvtot!(::Type{DeforModelRed3D}, t::Matrix{T}, v::Vector{T}) where {T}

Convert a 6-vector to a  matrix of 3x3 stress components (symmetric tensor).
"""
function stressvtot!(::Type{DeforModelRed3D}, t::Matrix{T}, v::Vector{T}) where {T}
	t[1,1] = v[1];
	t[2,2] = v[2];
	t[3,3] = v[3];
	t[1,2] = v[4];
	t[2,1] = v[4];
	t[1,3] = v[5];
	t[3,1] = v[5];
	t[3,2] = v[6];
	t[2,3] = v[6];
	return t
end

################################################################################
# 2-D plane strain model

"""
    dett(::Type{DeforModelRed2DStrain},  C::Matrix{T}) where {T}

Compute the determinant of a general square matrix.
"""
function dett(::Type{DeforModelRed2DStrain},  C::Matrix{T}) where {T}
	return (C[1,1] * C[2,2] - C[1,2] * C[2,1])
end

"""
    strainvdet(::Type{DeforModelRed2DStrain},  Cv::Vector{T}) where {T}

Compute the determinant of a symmetric strain-like square matrix represented
as a vector. Remember that the shear strain components are twice the entries
of the matrix representation.
"""
function strainvdet(::Type{DeforModelRed2DStrain},  Cv::Vector{T}) where {T}
	return (Cv[1] * Cv[2] - Cv[3]/2 * Cv[3]/2)
end

"""
    strainvtr(::Type{DeforModelRed2DStrain},  Cv::Vector{T}) where {T}

Compute the trace of a symmetric strain-like square matrix represented as a
vector.
"""
function strainvtr(::Type{DeforModelRed2DStrain},  Cv::Vector{T}) where {T}
	return (Cv[1] + Cv[2])
end

"""
    strainttov!(::Type{DeforModelRed2DStrain}, v::Vector{T}, t::Matrix{T}) where {T}

Convert a symmetric matrix of 2x2 strain components  into a 3-component vector.
"""
function strainttov!(::Type{DeforModelRed2DStrain}, v::Vector{T}, t::Matrix{T}) where {T}
    v[1] = t[1,1];
    v[2] = t[2,2];
    v[3] = t[1,2] + t[2,1];
    return v
end

"""
    strainvtot!(::Type{DeforModelRed2DStrain}, t::Matrix{T}, v::Vector{T}) where {T}

Convert a strain 3-vector to a  matrix of 2x2 strain components (symmetric tensor).
"""
function strainvtot!(::Type{DeforModelRed2DStrain}, t::Matrix{T}, v::Vector{T}) where {T}
    t[1,1] = v[1];
    t[2,2] = v[2];
    t[1,2] = v[3]/2.;
    t[2,1] = v[3]/2.;
    return t
end

"""
    stressttov!(::Type{DeforModelRed2DStrain}, v::Vector{T}, t::Matrix{T}) where {T}

Convert a symmetric matrix of 2x2 stress components to a 3-component vector.
"""
function stressttov!(::Type{DeforModelRed2DStrain}, v::Vector{T}, t::Matrix{T}) where {T}
	v[1] = t[1,1];
	v[2] = t[2,2];
	v[3] = 0.5*(t[1,2] + t[2,1]);
	return v
end

"""
    stressvtot!(::Type{DeforModelRed2DStrain}, t::Matrix{T}, v::Vector{T}) where {T}

Convert a vector to a  matrix of 2x2 stress components (symmetric tensor).

If `v` has 4 entries, also the `t[3,3]` matrix entry is set.

The stress vector components need to be ordered as:
    sigmax, sigmay, tauxy, sigmaz,
which is the ordering used for the plane-strain model reduction.
"""
function stressvtot!(::Type{DeforModelRed2DStrain}, t::Matrix{T}, v::Vector{T}) where {T}
    t[1,1] = v[1];
    t[2,2] = v[2];
    t[1,2] = v[3];
    t[2,1] = v[3];
    (length(v) == 4)  && (t[3,3] = v[4]);
    return t
end


################################################################################
# 2-D plane stress model

"""
    stressttov!(::Type{DeforModelRed2DStress}, v::Vector{T}, t::Matrix{T}) where {T}

Convert a symmetric matrix of 2x2 stress components to a 3-component vector.
"""
function stressttov!(::Type{DeforModelRed2DStress}, v::Vector{T}, t::Matrix{T}) where {T}
	v[1] = t[1,1];
	v[2] = t[2,2];
	v[3] = 0.5*(t[1,2] + t[2,1]);
	return v
end

"""
    stressvtot!(::Type{DeforModelRed2DStress}, t::Matrix{T}, v::Vector{T}) where {T}

Convert a 3-vector to a  matrix of 2x2 stress components (symmetric tensor).
"""
function stressvtot!(::Type{DeforModelRed2DStress}, t::Matrix{T}, v::Vector{T}) where {T}
    t[1,1] = v[1];
    t[2,2] = v[2];
    t[1,2] = v[3];
    t[2,1] = v[3];
    return t
end

################################################################################
# 2-D axially symmetric stress model

"""
    stressvtot!(::Type{DeforModelRed2DAxisymm}, t::Matrix{T}, v::Vector{T}) where {T}

Convert a 4-vector to a  matrix of 3x3 stress components (tensor).

Convert a 4-vector to a *symmetric* matrix of 3x3 stress components (tensor).

The stress vector components need to be ordered as:
    sigmax, sigmay, sigmaz, tauxy.
"""
function stressvtot!(::Type{DeforModelRed2DAxisymm}, t::Matrix{T}, v::Vector{T}) where {T}
    t[1,1] = v[1];
    t[2,2] = v[2];
    t[1,2] = v[4];
    t[2,1] = v[4];
    t[3,3] = v[3];
    return t
end
