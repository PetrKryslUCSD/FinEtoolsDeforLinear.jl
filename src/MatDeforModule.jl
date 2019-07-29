"""
    MatDeforModule

Module to support general operations for deformation material models.
"""
module MatDeforModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D, DeforModelRed2DStrain, DeforModelRed2DStress, DeforModelRed2DAxisymm, DeforModelRed1D
import FinEtools.MatModule: AbstractMat
using LinearAlgebra

_RotationMatrix = Union{Array{T, 2}, Transpose{T, Array{T, 2}}} where {T}

"""
    AbstractMatDefor

Abstract type that represents deformable materials.
"""
abstract type AbstractMatDefor <: AbstractMat end

"""
    strain2x2tto3v!(v::FVec{T}, t::FMat{T}) where {T}

Convert a matrix of 2x2 strain components  into a 3-component vector.
"""
function strain2x2tto3v!(v::FVec{T}, t::FMat{T}) where {T}
    v[1] = t[1,1];
    v[2] = t[2,2];
    v[3] = t[1,2] + t[2,1];
    return v
end

"""
    strain3vto2x2t!(t::FMat{T}, v::FVec{T}) where {T}

Convert a strain 3-vector to a  matrix of 2x2 strain components (symmetric tensor).
"""
function strain3vto2x2t!(t::FMat{T}, v::FVec{T}) where {T}
    t[1,1] = v[1];
    t[2,2] = v[2];
    t[1,2] = v[3]/2.;
    t[2,1] = v[3]/2.;
    return t
end

"""
    strain3x3tto6v!(v::FVec{T}, t::FMat{T}) where {T}

Convert a matrix of 3x3 strain components to a 6-component vector.
"""
function strain3x3tto6v!(v::FVec{T}, t::FMat{T}) where {T}
    v[1] = t[1,1];
    v[2] = t[2,2];
    v[3] = t[3,3];
    v[4] = t[1,2] + t[2,1];
    v[5] = t[1,3] + t[3,1];
    v[6] = t[3,2] + t[2,3];
    return v
end

"""
    strain6vto3x3t!(t::FMat{T}, v::FVec{T}) where {T}

Convert a strain 6-vector to a  matrix of 3x3 strain components (symmetric tensor)..
"""
function strain6vto3x3t!(t::FMat{T}, v::FVec{T}) where {T}
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
    strain3x3tto9v!(v::FVec{T}, t::FMat{T}) where {T}

Convert a matrix of 3x3 strain components to a 9-component vector.

The  strain components are in the order
     ex, ey, ez, gxy, gyx, gyz, gzy, gxz, gzx
"""
function strain3x3tto9v!(v::FVec{T}, t::FMat{T}) where {T}
    v[1] = t[1,1];
    v[2] = t[2,2];
    v[3] = t[3,3];
    v[4] = t[1,2];
    v[5] = t[2,1];
    v[6] = t[2,3];
    v[7] = t[3,2];
    v[8] = t[1,3];
    v[9] = t[3,1];
    return v
end

"""
    strain9vto3x3t!(v::FVec{T}, t::FMat{T}) where {T}

Convert a matrix of 3x3 strain components to a 9-component vector.

The  strain components are in the order
     ex, ey, ez, gxy, gyx, gyz, gzy, gxz, gzx
"""
function strain9vto3x3t!(t::FMat{T}, v::FVec{T}) where {T}
    t[1,1] = v[1];
    t[2,2] = v[2];
    t[3,3] = v[3];
    t[1,2] = v[4];
    t[2,1] = v[5];
    t[2,3] = v[6];
    t[3,2] = v[7];
    t[1,3] = v[8];
    t[3,1] = v[9];
    return t
end

"""
    stress2x2to3v!(v::FVec{T}, t::FMat{T})

Convert a symmetric matrix of 2x2 stress components to a 3-component vector.
"""
function stress2x2to3v!(v::FVec{T}, t::FMat{T}) where {T}
    v[1] = t[1,1];
    v[2] = t[2,2];
    v[3] = 1.0/2.0*(t[1,2] + t[2,1]);
    return v
end

"""
    stress3vto2x2t!(t::FMat{T}, v::FVec{T})

Convert a 3-vector to a  matrix of 2x2 stress components (symmetric tensor).
"""
function stress3vto2x2t!(t::FMat{T}, v::FVec{T}) where {T}
    t[1,1] = v[1];
    t[2,2] = v[2];
    t[1,2] = v[3];
    t[2,1] = v[3];
    return t
end

"""
    stress3vto3x3t!(t::FMat{T}, v::FVec{T}) where {T}

Convert a 3-vector to a matrix of 3x3 stress components (symmetric tensor).
"""
function stress3vto3x3t!(t::FMat{T}, v::FVec{T}) where {T}
    t[1,1] = v[1];
    t[2,2] = v[2];
    t[1,2] = v[3];
    t[2,1] = v[3];
    return t
end

"""
    stress4vto3x3t!(t::FMat{T}, v::FVec{T}) where {T}

Convert a 4-vector to a  matrix of 3x3 stress components (tensor).

Convert a 4-vector to a *symmetric*
matrix of 3x3 stress components (tensor).  This is
conversion routine that would be useful for plane strain or
axially symmetric conditions.
The stress vector components need to be ordered as:
    sigmax, sigmay, tauxy, sigmaz,
which is the ordering used for the plane-strain model reduction.
Therefore, for axially symmetric analysis the components need to be
reordered, as from the constitutive equation they come out
as sigmax, sigmay, sigmaz, tauxy.
"""
function stress4vto3x3t!(t::FMat{T}, v::FVec{T}) where {T}
    t[1,1] = v[1];
    t[2,2] = v[2];
    t[1,2] = v[3];
    t[2,1] = v[3];
    t[3,3] = v[4];
    return t
end

"""
    stress6vto3x3t!(t::FMat{T}, v::FVec{T}) where {T}

Convert a 6-vector to a  matrix of 3x3 stress components (symmetric tensor).
"""
function stress6vto3x3t!(t::FMat{T}, v::FVec{T}) where {T}
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

"""
    stress3x3tto6v!(v::FVec{T}, t::FMat{T}) where {T}

Convert a matrix of 3x3 stress components to a 6-component vector.
"""
function stress3x3tto6v!(v::FVec{T}, t::FMat{T}) where {T}
    v[1] = t[1,1];
    v[2] = t[2,2];
    v[3] = t[3,3];
    v[4] = 1/2.0*(t[1,2] + t[2,1]);
    v[5] = 1/2.0*(t[1,3] + t[3,1]);
    v[6] = 1/2.0*(t[3,2] + t[2,3]);
    return v
end

"""
    strain9vto6v!(t::FVec{T}, v::FVec{T}) where {T}

Convert a strain 9-vector to a  strain 6-vector components (tensor).
"""
function strain9vto6v!(t::FVec{T}, v::FVec{T}) where {T}
    t[1] = v[1];
    t[2] = v[2];
    t[3] = v[3];
    t[4] = v[4]+v[5];
    t[5] = v[8]+v[9];
    t[6] = v[6]+v[7];
    return t
end

"""
    strain6vto9v!(t::FVec{T}, v::FVec{T}) where {T}

Convert a strain 6-vector to a  strain 9-vector components (tensor).

The  strain components are in the order
     ex, ey, ez, gxy/2, gxy/2, gyz/2, gyz/2, gxz/2, gxz/2
"""
function strain6vto9v!(t::FVec{T}, v::FVec{T}) where {T}
    t[1] = v[1];
    t[2] = v[2];
    t[3] = v[3];
    t[4] = v[4]/2.;
    t[5] = v[4]/2.;
    t[6] = v[6]/2.;
    t[7] = v[6]/2.;
    t[8] = v[5]/2.;
    t[9] = v[5]/2.;
    return t
end

"""
    stress9vto6v!(t::FVec{T}, v::FVec{T}) where {T}

Convert a stress 9-vector (tensor) to a  stress 6-vector components.

The  stress components are in the order
     sigx, sigy, sigz, tauxy, tauxy, tauyz, tauyz, tauxz, tauxz
"""
function stress9vto6v!(t::FVec{T}, v::FVec{T}) where {T}
    t[1] = v[1];
    t[2] = v[2];
    t[3] = v[3];
    t[4] = v[4];
    t[5] = v[8];
    t[6] = v[6];
    return t
end

"""
    stress6vto9v!(t::FVec{T}, v::FVec{T}) where {T}

Convert a stress 6-vector to a  stress 9-vector components (tensor).

The  stress components are in the order
     sigx, sigy, sigz, tauxy, tauxy, tauyz, tauyz, tauxz, tauxz
"""
function stress6vto9v!(t::FVec{T}, v::FVec{T}) where {T}
    t[1] = v[1];
    t[2] = v[2];
    t[3] = v[3];
    t[4] = v[4];
    t[5] = v[4];
    t[6] = v[6];
    t[7] = v[6];
    t[8] = v[5];
    t[9] = v[5];
    return t
end


################################################################################
# 3-D model

"""
    rotstressvec!(::Type{DeforModelRed3D},  outstress::FVec{T},  instress::FVec{T}, Rm::_RotationMatrix) where {T}

Rotate the stress vector by the supplied rotation matrix.

Calculate the rotation of the stress vector to the 'bar' coordinate system given
by the columns of the rotation matrix `Rm`.

- `outstress` = output stress vector, overwritten inside
- `instress` = input stress vector
- `Rm` = columns are components of 'bar' basis vectors on the 'plain'
     basis vectors
"""
function  rotstressvec!(::Type{DeforModelRed3D},  outstress::FVec{T},  instress::FVec{T}, Rm::_RotationMatrix) where {T}
      # # Derivation of the transformation matrix [T]
    # #This is from Barbero''s  book Finite element analysis of composite
    # #materials  using Abaqus.  Note that his matrix "a"  is the transpose of
    # #the FinEALE matrix "Rm".
    # #  We also use the FinEALE numbering of the strains.
    #
    # syms T alpha R
    # syms a a11 a12 a13 a21 a22 a23 a31 a32 a33 real
    # a = [a11,a12,a13;
    #     a21,a22,a23;
    #     a31,a32,a33];
    # a = a';#his matrix "a"  is the transpose of the FinEALE matrix "Rm".
    # # # it can be done in terms of l,m,n's as well
    # # syms a l1 m1 n1 l2 m2 n2 l3 m3 n3
    # # a = [l1,m1,n1;l2,m2,n2;l3,m3,n3]
    # #  We also use the FinEALE numbering of the strains.
    # Numbering =[1,4,5;
    #     4,2,6;
    #     5,6,3];
    # T(1:6,1:6) = 0;
    # for i=1:1:3
    #     for j=1:1:3
    #         #if i==j; alpha = j; else alpha = 9-i-j; end
    #         alpha = Numbering(i,j);
    #         for p=1:1:3
    #             for q=1:1:3
    #                 #   if p==q beta = p; else beta = 9-p-q; end
    #                 beta = Numbering(p,q);
    #                 T(alpha,beta) = 0;
    #                 if alpha<=3 & beta<= 3; T(alpha,beta)=a(i,p)*a(i,p); end
    #                 if alpha> 3 & beta<= 3; T(alpha,beta)=a(i,p)*a(j,p); end
    #                 if alpha<=3 & beta>3; T(alpha,beta)=a(i,q)*a(i,p)+a(i,p)*a(i,q);end
    #                 if alpha>3 & beta>3; T(alpha,beta)=a(i,p)*a(j,q)+a(i,q)*a(j,p);end
    #             end
    #         end
    #     end
    # end
    # T
    # R = eye(6,6); R(4,4)=2; R(5,5)=2; R(6,6)=2; # Reuter matrix
    # Tbar = R*T*R^(-1)
    a11=Rm[1,1]; a12=Rm[1,2]; a13=Rm[1,3];
    a21=Rm[2,1]; a22=Rm[2,2]; a23=Rm[2,3];
    a31=Rm[3,1]; a32=Rm[3,2]; a33=Rm[3,3];

    outstress[1] =  (a11^2)*instress[1] + (a21^2)*instress[2] +
                    (a31^2)*instress[3] + (2*a11*a21)*instress[4] +
                    (2*a11*a31)*instress[5] + (2*a21*a31)*instress[6]
    outstress[2] =  (a12^2)*instress[1] +  (a22^2)*instress[2] +
                    (a32^2)*instress[3] + (2*a12*a22)*instress[4] +
                    (2*a12*a32)*instress[5] + (2*a22*a32)*instress[6]
    outstress[3] =  (a13^2)*instress[1] + (a23^2)*instress[2] +
                    (a33^2)*instress[3] + (2*a13*a23)*instress[4] +
                    (2*a13*a33)*instress[5] + (2*a23*a33)*instress[6]
    outstress[4] =  a11*a12*instress[1] +  a21*a22*instress[2] +
                    a31*a32*instress[3] + (a11*a22 + a12*a21)*instress[4] +
                    (a11*a32 + a12*a31)*instress[5] + (a21*a32 + a22*a31)*instress[6]
    outstress[5] =  (a11*a13)*instress[1] + (a21*a23)*instress[2] +
                    (a31*a33)*instress[3] + (a11*a23 + a13*a21)*instress[4] +
                    (a11*a33 + a13*a31)*instress[5] + (a21*a33 + a23*a31)*instress[6]
    outstress[6] =  (a12*a13)*instress[1] + (a22*a23)*instress[2] +
                    (a32*a33)*instress[3] + (a12*a23 + a13*a22)*instress[4] +
                    (a12*a33 + a13*a32)*instress[5] + (a22*a33 + a23*a32)*instress[6]
    return outstress
end



################################################################################
# 2-D plane strain model

"""
    rotstressvec!(::Type{DeforModelRed2DStrain},  outstress::FVec{T},  instress::FVec{T},  Rm::_RotationMatrix) where {T}

Rotate the stress vector by the supplied rotation matrix.

Calculate the rotation of the stress vector to the 'bar' coordinate system given
by the columns of the rotation matrix `Rm`.

- `outstress` = output stress vector, overwritten inside
- `instress` = input stress vector
- `Rm` = columns are components of 'bar' basis vectors on the 'plain'
     basis vectors
"""
function  rotstressvec!(::Type{DeforModelRed2DStrain},  outstress::FVec{T},  instress::FVec{T},  Rm::_RotationMatrix) where {T}
    a11=Rm[1,1]; a12=Rm[1,2]; a13=0.0;
    a21=Rm[2,1]; a22=Rm[2,2]; a23=0.0;
    a31=0.0; a32=0.0; a33=1.0;
    # Note the special arrangement  of the components for plane strain
    outstress[1] =  (a11^2)*instress[1] + (a21^2)*instress[2] + (a31^2)*instress[4] + (2*a11*a21)*instress[3]
    outstress[2] =  (a12^2)*instress[1] +  (a22^2)*instress[2] + (a32^2)*instress[4] + (2*a12*a22)*instress[3]
    outstress[4] =  (a13^2)*instress[1] + (a23^2)*instress[2] + (a33^2)*instress[4] + (2*a13*a23)*instress[3]
    outstress[3] =  a11*a12*instress[1] +  a21*a22*instress[2] + a31*a32*instress[4] + (a11*a22 + a12*a21)*instress[3]
    return outstress
end



################################################################################
# 2-D plane stress model

"""
    rotstressvec!(::Type{DeforModelRed2DStress},  outstress::FVec{T},  instress::FVec{T},  Rm::_RotationMatrix) where {T}

Rotate the stress vector by the supplied rotation matrix.

Calculate the rotation of the stress vector to the 'bar' coordinate system given
by the columns of the rotation matrix `Rm`.

- `outstress` = output stress vector, overwritten inside
- `instress` = input stress vector
- `Rm` = columns are components of 'bar' basis vectors on the 'plain'
     basis vectors
"""
function  rotstressvec!(::Type{DeforModelRed2DStress},  outstress::FVec{T},  instress::FVec{T},  Rm::_RotationMatrix) where {T}
    a11=Rm[1,1]; a12=Rm[1,2];
    a21=Rm[2,1]; a22=Rm[2,2];
    outstress[1] =  (a11^2)*instress[1] + (a21^2)*instress[2] + (2*a11*a21)*instress[3]
    outstress[2] =  (a12^2)*instress[1] +  (a22^2)*instress[2] + (2*a12*a22)*instress[3]
    outstress[3] =  (a11*a12)*instress[1] + (a21*a22)*instress[2] + (a11*a22+a12*a21)*instress[3]
    return outstress
end

################################################################################
# 2-D axially symmetric stress model

"""
    rotstressvec!(::Type{DeforModelRed2DAxisymm},  outstress::FVec{T},  instress::FVec{T},  Rm::_RotationMatrix) where {T}

Rotate the stress vector by the supplied rotation matrix.

Calculate the rotation of the stress vector to the 'bar' coordinate system given
by the columns of the rotation matrix `Rm`.

- `outstress` = output stress vector, overwritten inside
- `instress` = input stress vector
- `Rm` = columns are components of 'bar' basis vectors on the 'plain'
     basis vectors
"""
function  rotstressvec!(::Type{DeforModelRed2DAxisymm},  outstress::FVec{T},  instress::FVec{T},  Rm::_RotationMatrix) where {T}
    a11=Rm[1,1]; a12=Rm[1,2];
    a21=Rm[2,1]; a22=Rm[2,2];
    outstress[1] =  (a11^2)*instress[1] + (a21^2)*instress[2] + (0.0)*instress[3] + (2*a11*a21)*instress[4]
    outstress[2] =  (a12^2)*instress[1] + (a22^2)*instress[2] + (0.0)*instress[3] + (2*a12*a22)*instress[4]
    outstress[3] =  (0.0)*instress[1] + (0.0)*instress[2] + (1.0)*instress[3] + (0.0)*instress[4]
    outstress[4] =  (a11*a12)*instress[1] + (a21*a22)*instress[2] + (0.0)*instress[3] + (a11*a22 + a12*a21)*instress[4]
    return outstress
end


################################################################################
# 1-D stress model

"""
    rotstressvec!(::Type{DeforModelRed1D},  outstress::FVec{T},  instress::FVec{T},  Rm::_RotationMatrix) where {T}

Rotate the stress vector by the supplied rotation matrix.

Calculate the rotation of the stress vector to the 'bar' coordinate system given
by the columns of the rotation matrix `Rm`.

- `outstress` = output stress vector, overwritten inside
- `instress` = input stress vector
- `Rm` = columns are components of 'bar' basis vectors on the 'plain'
     basis vectors
"""
function  rotstressvec!(::Type{DeforModelRed1D},  outstress::FVec{T},  instress::FVec{T},  Rm::_RotationMatrix) where {T}
    copyto!(outstress, instress)
    return outstress
end

"""
    dett(::Type{DeforModelRed3D},  C::FMat{T}) where {T}

Compute the determinant of the general square matrix.
"""
function dett(::Type{DeforModelRed3D},  C::FMat{T}) where {T}
	return (C[1,1] * C[2,2] * C[3,3] + 
		C[1,2] * C[2,3] * C[3,1] + 
		C[1,3] * C[2,1] * C[3,2] - 
		C[1,3] * C[2,2] * C[3,1] - 
		C[1,2] * C[2,1] * C[3,3] - 
		C[1,1] * C[2,3] * C[3,2])
end

"""
    strain6vdet(::Type{DeforModelRed3D},  Cv::FVec{T}) where {T}

Compute the determinant of a symmetric strain-like square matrix represented
as a vector.
"""
function strain6vdet(::Type{DeforModelRed3D},  Cv::FVec{T}) where {T}
	return (Cv[1] * Cv[2] * Cv[3] + 
		Cv[4]/2 * Cv[6]/2 * Cv[5]/2 + 
		Cv[5]/2 * Cv[4]/2 * Cv[6]/2 -  
		Cv[5]/2 * Cv[2] * Cv[5]/2 - 
		Cv[4]/2 * Cv[4]/2 * Cv[3] - 
		Cv[1] * Cv[6]/2 * Cv[6]/2)
end

"""
    strain6vtr(::Type{DeforModelRed3D},  Cv::FVec{T}) where {T}

Compute the trace of a symmetric strain-like square matrix represented as a
vector.
"""
function strain6vtr(::Type{DeforModelRed3D},  Cv::FVec{T}) where {T}
	return (Cv[1] + Cv[2] + Cv[3])
end

"""
    tens4symmto6x6t!(M::FMat{T}, ST::FMat{T}) where {T}

Convert a symmetric 4th-order tensor to a 6 x 6 matrix.

# Example

J=tens4_ijkl(eye(3),eye(3))
produces the tracor:
T=rand(3); 
sum(diag(T))*eye(3)
t= tens4_dot_2(J,T)
M= tens4_symm_to_6(ST)
"""
function tens4symmto6x6t!(M::FMat{T}, ST::FMat{T}) where {T}
	# This corresponds to the arrangement of the components of stress (or
	# strain) tensor, symmetric, three-dimensional, into a 6-component 
	# vector.
	ix=[1 1; 2 2; 3 3; 1 2; 1 3; 2 3];
	for  j in 1:6
		for  i in 1:6
			M[i,j] = ST[ix[i,1],ix[i,2],ix[j,1],ix[j,2]];
		end
	end
	return M
end

end # module

# function stressvectorrotation{MR<:DeforModelRed2DStress}(::Type{MR},
#   Rm::FMat{T})
#
#     # Rm = columns are components of 'bar' basis vectors on the 'plain'
#     #      basis vectors
#     #
#     # Calculate the rotation of the 'plain' stress vector to the
#     # 'bar' coordinate system given by the columns of the rotation matrix Rm.
#     #
#     # Example:
#     # The stress vector "stress" is given in the material coordinate
#     # system defined by the orientation matrix Rm. The following two
#     # transformations are equivalent:
#     #
#     #         t = stress_6v_to_3x3t (mat,stress);
#     #         t = (Rm*t*Rm');# in global coordinate system
#     #         t = (outputRm'*t*outputRm);# in output coordinate system
#     #         stress =stress_3x3t_to_6v (mat,t);# in output coordinate system
#     #
#     #        stress =mat.stress_vector_rotation(outputRm)...
#     #                  *mat.stress_vector_rotation(Rm')...
#     #                      *stress;# in output coordinate system
#
#     a11=Rm[1,1]; a12=Rm[1,2];
#     a21=Rm[2,1]; a22=Rm[2,2];
#     # switch self.reduction
#     #     case {'axisymm','strain'}
#     #         T =[[   a11^2,   a21^2, 0,         2*a11*a21]
#     #             [   a12^2,   a22^2, 0,         2*a12*a22]
#     #             [       0,       0, 1,                 0]
#     #             [ a11*a12, a21*a22, 0, a11*a22 + a12*a21]];
#     #             case 'stress'
#     T =[[   a11^2   a21^2       2*a11*a21]
#         [   a12^2   a22^2       2*a12*a22]
#         [ a11*a12 a21*a22 a11*a22+a12*a21]];
#     return T
# end
#
# function strainvectorrotation{MR<:DeforModelRed2DStress}(::Type{MR},
#   Rm::FMat{T})
#     # Calculate the rotation matrix for a strain vector.
#     #
#     #   Tbar = strain_vector_rotation(self,Rm)
#     #
#     # Rm = columns are components of 'bar' basis vectors on the 'plain'
#     #      basis vectors
#     #
#     # Calculate the rotation of the 'plain' strain vector to the
#     # 'bar' coordinate system given by the columns of the rotation matrix Rm.
#     #
#
#     a11=Rm[1,1]; a12=Rm[1,2];
#     a21=Rm[2,1]; a22=Rm[2,2];
#     # switch self.reduction
#     #     case {'axisymm','strain'}
#     #         Tbar =[ [     a11^2,     a21^2, 0,           a11*a21]
#     #                 [     a12^2,     a22^2, 0,           a12*a22]
#     #                 [         0,         0, 1,                 0]
#     #                 [ 2*a11*a12, 2*a21*a22, 0, a11*a22 + a12*a21]];
#     #   case 'stress'
#     Tbar =[ [     a11^2,     a21^2,           a11*a21]
#            [     a12^2,     a22^2,           a12*a22]
#            [ 2*a11*a12, 2*a21*a22, a11*a22 + a12*a21]];
#     return Tbar
# end
#
#
# ################################################################################
# # 2-D plane axially symmetric model
#
# function stressvectorrotation(::Type{MR},
#   Rm::FMat{T}) where {MR<:DeforModelRed2DAxisymm}
#
#     # Rm = columns are components of 'bar' basis vectors on the 'plain'
#     #      basis vectors
#     #
#     # Calculate the rotation of the 'plain' stress vector to the
#     # 'bar' coordinate system given by the columns of the rotation matrix Rm.
#     #
#     # Example:
#     # The stress vector "stress" is given in the material coordinate
#     # system defined by the orientation matrix Rm. The following two
#     # transformations are equivalent:
#     #
#     #         t = stress_6v_to_3x3t (mat,stress);
#     #         t = (Rm*t*Rm');# in global coordinate system
#     #         t = (outputRm'*t*outputRm);# in output coordinate system
#     #         stress =stress_3x3t_to_6v (mat,t);# in output coordinate system
#     #
#     #        stress =mat.stress_vector_rotation(outputRm)...
#     #                  *mat.stress_vector_rotation(Rm')...
#     #                      *stress;# in output coordinate system
#
#     a11=Rm[1,1]; a12=Rm[1,2];
#     a21=Rm[2,1]; a22=Rm[2,2];
#     T =[[   a11^2    a21^2  0          2*a11*a21]
#         [   a12^2    a22^2  0          2*a12*a22]
#         [       0        0  1                  0]
#         [ a11*a12  a21*a22  0  a11*a22 + a12*a21]];
#     return T
# end
#
# function strainvectorrotation(::Type{MR},
#   Rm::FMat{T}) where {MR<:DeforModelRed2DAxisymm}
#     # Calculate the rotation matrix for a strain vector.
#     #
#     #   Tbar = strain_vector_rotation(self,Rm)
#     #
#     # Rm = columns are components of 'bar' basis vectors on the 'plain'
#     #      basis vectors
#     #
#     # Calculate the rotation of the 'plain' strain vector to the
#     # 'bar' coordinate system given by the columns of the rotation matrix Rm.
#     #
#
#     a11=Rm[1,1]; a12=Rm[1,2];
#     a21=Rm[2,1]; a22=Rm[2,2];
#     Tbar =[ [     a11^2      a21^2  0            a11*a21]
#            [     a12^2      a22^2  0            a12*a22]
#            [         0          0  1                  0]
#            [ 2*a11*a12  2*a21*a22  0  a11*a22 + a12*a21]];
#     return Tbar
# end
#
# ################################################################################
# # Generic versions of rotations of stiffness and compliance matrices
#
# function rotatestiffness!(::Type{MR}, D::FMat{T},
#                                                            Rm::FMat{T}) where {MR<:DeforModelRed}
#     # Rotate constitutive stiffness matrix of the material.
#     #
#     #         function D=transform_stiffness(self,D,Rm)
#     #
#     # Rotate constitutive stiffness matrix of the material to the
#     # coordinate system given by the columns of the rotation matrix Rm.
#     T =stressvectorrotation(MR,Rm);
#     D = T*D*T';
#     return D
# end

#
# function rotatecompliance!(::Type{MR}, C::FMat{T},
#   Rm::FMat{T}) where {MR<:DeforModelRed}
#     # Rotate constitutive compliance matrix of the material.
#     #
#     #   C = rotate_compliance(self,C,Rm)
#     #
#     # Rotate constitutive compliance matrix of the material to the
#     # coordinate system given by the columns of the rotation matrix Rm.
#     Tbar =strainvectorrotation(MR,Rm);
#     C = Tbar*C*Tbar';
#     return C
# end
# function stressvectorrotation(::Type{DeforModelRed2DStrain}, Rm::FMat{T})
#
#     # Rm = columns are components of 'bar' basis vectors on the 'plain'
#     #      basis vectors
#     #
#     # Calculate the rotation of the 'plain' stress vector to the
#     # 'bar' coordinate system given by the columns of the rotation matrix Rm.
#     #
#     # Example:
#     # The stress vector "stress" is given in the material coordinate
#     # system defined by the orientation matrix Rm. The following two
#     # transformations are equivalent:
#     #
#     #         t = stress_6v_to_3x3t (mat,stress);
#     #         t = (Rm*t*Rm');# in global coordinate system
#     #         t = (outputRm'*t*outputRm);# in output coordinate system
#     #         stress =stress_3x3t_to_6v (mat,t);# in output coordinate system
#     #
#     #        stress =mat.stress_vector_rotation(outputRm)...
#     #                  *mat.stress_vector_rotation(Rm')...
#     #                      *stress;# in output coordinate system
#
#     a11=Rm[1,1]; a12=Rm[1,2];
#     a21=Rm[2,1]; a22=Rm[2,2];
#     # sigmax, sigmay, sigmaz, tauxy
#     T =[[   a11^2    a21^2  0          2*a11*a21]
#         [   a12^2    a22^2  0          2*a12*a22]
#         [       0        0  1                  0]
#         [ a11*a12  a21*a22  0  a11*a22 + a12*a21]];
#     return T
# end
#
# function strainvectorrotation(::Type{DeforModelRed2DStrain},  Rm::FMat{T})
#     # Calculate the rotation matrix for a strain vector.
#     #
#     #   Tbar = strain_vector_rotation(self,Rm)
#     #
#     # Rm = columns are components of 'bar' basis vectors on the 'plain'
#     #      basis vectors
#     #
#     # Calculate the rotation of the 'plain' strain vector to the
#     # 'bar' coordinate system given by the columns of the rotation matrix Rm.
#     #
#
#     a11=Rm[1,1]; a12=Rm[1,2];
#     a21=Rm[2,1]; a22=Rm[2,2];
#     Tbar =[ [     a11^2      a21^2  0            a11*a21]
#            [     a12^2      a22^2  0            a12*a22]
#            [         0          0  1                  0]
#            [ 2*a11*a12  2*a21*a22  0  a11*a22 + a12*a21]];
#     return Tbar
# end
# function strainvectorrotation(::Type{DeforModelRed3D}, Rm::FMat{T})
#     # Calculate the rotation matrix for a strain vector.
#     #
#     #   Tbar = strain_vector_rotation(self,Rm)
#     #
#     # Rm = columns are components of 'bar' basis vectors on the 'plain'
#     #      basis vectors
#     #
#     # Calculate the rotation of the 'plain' strain vector to the
#     # 'bar' coordinate system given by the columns of the rotation matrix Rm.
#     #
#
#     a11=Rm[1,1]; a12=Rm[1,2]; a13=Rm[1,3];
#     a21=Rm[2,1]; a22=Rm[2,2]; a23=Rm[2,3];
#     a31=Rm[3,1]; a32=Rm[3,2]; a33=Rm[3,3];
#     Tbar =[
#            [     a11^2     a21^2     a31^2           a11*a21           a11*a31           a21*a31]
#            [     a12^2     a22^2     a32^2           a12*a22           a12*a32           a22*a32]
#            [     a13^2     a23^2     a33^2           a13*a23           a13*a33           a23*a33]
#            [ 2*a11*a12 2*a21*a22 2*a31*a32 a11*a22 + a12*a21 a11*a32 + a12*a31 a21*a32 + a22*a31]
#            [ 2*a11*a13 2*a21*a23 2*a31*a33 a11*a23 + a13*a21 a11*a33 + a13*a31 a21*a33 + a23*a31]
#            [ 2*a12*a13 2*a22*a23 2*a32*a33 a12*a23 + a13*a22 a12*a33 + a13*a32 a22*a33 + a23*a32]];
#     return Tbar
# end
