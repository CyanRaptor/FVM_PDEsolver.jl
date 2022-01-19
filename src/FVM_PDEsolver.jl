module FVM_PDEsolver

using ForwardDiff
using DifferentialEquations
using LinearAlgebra

include("FVMPDE_grid.jl")
export FVMPDEGrid

include("FVMPDE_problem.jl")
export FVMPDEProblem

include("FVMPDE_methods.jl")
export uw1,uw2,uw3,uw4,cd,fr,kn,sb,mm,mu,ha,va1,va2,vl,op,hc,hq,cm,mc,sm,um,gg


include("FVMPDE_solver.jl")
export FVMPDESolve

end # module
