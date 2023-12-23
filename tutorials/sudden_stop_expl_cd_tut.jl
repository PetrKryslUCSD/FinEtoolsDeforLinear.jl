# # Suddenly-stopped bar: Centered difference explicit

# Source code: [`sudden_stop_expl_cd_tut.jl`](sudden_stop_expl_cd_tut.jl)

# ## Description


# A bar is given an initial velocity and then at time 0.0 it is
# suddenly stopped by fixing one of its ends. This sends a wave down the bar.
# The output of the simulation is the velocity, which tends to reproduce 
# the rectangular pulses in which the velocity bounces back and forth.

# The beam is modeled as a solid. The classical centered difference 
# rule is used to integrate the
# equations of motion in time. No damping is present. 

# ## Goals

# - Show how to create the discrete model for explicit dynamics.
# - Demonstrate centered difference explicit time stepping.

##
# ## Definitions

# Basic imports.
using LinearAlgebra
using Arpack

# This is the finite element toolkit itself.
using FinEtools

# The linear stress analysis application is implemented in this package.
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule

# Input parameters
E = 205000*phun("MPa");# Young's modulus
nu = 0.3;# Poisson ratio
rho = 7850*phun("KG*M^-3");# mass density
L = 20*phun("mm");
W = 1*phun("mm");
H = 1*phun("mm");
tolerance = W/500;
vmag = 0.1*phun("m")/phun("SEC");
tend = 0.00005*phun("SEC");
    
##
# ## Create the discrete model

MR = DeforModelRed3D
fens,fes  = H8block(L,W,H, 80,4,4)

geom = NodalField(fens.xyz)
u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

nl = selectnode(fens, box=[0 0 -Inf Inf -Inf Inf], inflate=tolerance)
setebc!(u, nl, true, 1)
applyebc!(u)
numberdofs!(u)

corner = selectnode(fens, nearestto=[L 0 0])
cornerxdof = u.dofnums[corner[1], 1]

material = MatDeforElastIso(MR, rho, E, nu, 0.0)

femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,2)), material)
femm = associategeometry!(femm, geom)
K = stiffness(femm, geom, u)
# Assemble the mass matrix as diagonal. The HRZ lumping technique is 
# applied through the assembler of the sparse matrix.
hrzass = SysmatAssemblerSparseHRZLumpingSymm(0.0)
femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,3)), material)
M = mass(femm, hrzass, geom, u)

# Extract the free-free block of the matrices.
M_ff = matrix_blocked(M, nfreedofs(u))[:ff]
K_ff = matrix_blocked(K, nfreedofs(u))[:ff]

# Figure out the highest frequency in the model, and use a time step that is
# smaller than the period of the highest frequency.
evals, evecs = eigs(K_ff, M_ff; nev=1, which=:LM);
@show dt = 0.9 * 2/real(sqrt(evals[1]));

# The time stepping loop is protected by `let end` to avoid unpleasant surprises
# with variables getting clobbered by globals.
ts, cornervxs = let dt = dt
    # Initial displacement, velocity, and acceleration.
    U0 = gathersysvec(u)
    v = deepcopy(u)
    v.values[:, 1] .= -vmag
    V0 = gathersysvec(v)
    F1 = fill(0.0, length(V0))
    U1 = fill(0.0, length(V0))
    V1 = fill(0.0, length(V0))
    A0 = fill(0.0, length(V0))
    A1 = fill(0.0, length(V0))
    phi = 1.005
    # The times and displacements of the corner will be collected into two vectors
    ts = Float64[]
    cornervxs = Float64[]
    # Let us begin the time integration loop:
    t = 0.0; 
    step = 0;
    while t < tend
        push!(ts, t)
        push!(cornervxs, V0[cornerxdof])
        t = t+dt;
        step = step + 1;
        (mod(step,1000)==0) && println("Step $(step): $(t)")
        # Zero out the load
        fill!(F1, 0.0);
        # Initial acceleration
        if step == 1
            A0 = M_ff \ (F1)
        end
        # Update displacement.
        @. U1 = U0 + dt*V0 + (dt^2/2)*A0;
        # Compute updated acceleration.
        A1 .= M_ff \ (-K_ff*U1 + F1)
        # Update the velocities.
        @. V1 = V0 + (dt/2)*(A0 + A1)
        # Switch the temporary vectors for the next step.
        U0, U1 = U1, U0;
        V0, V1 = V1, V0;
        A0, A1 = A1, A0;
        if (t == tend) # Are we done yet?
            break;
        end
        if (t+dt > tend) # Adjust the last time step so that we exactly reach tend
            dt = tend-t;
        end
    end
    ts, cornervxs # return the collected results
end

##
# ## Plot the results

using Gnuplot

@gp "set terminal windows 7 " :-
@gp  :- ts cornervxs./phun("mm") "lw 2 lc rgb 'red' with lines title 'Displacement of the corner' " 
@gp  :- "set xlabel 'Time [s]'"
@gp  :- "set ylabel 'Velocity [mm/s]'"


# The end.
true
