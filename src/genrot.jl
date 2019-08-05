
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
