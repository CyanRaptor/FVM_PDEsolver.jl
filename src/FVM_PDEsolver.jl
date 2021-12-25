module FVM_PDEsolver

using ForwardDiff
using DifferentialEquations
using LinearAlgebra

include("FVMPDE_grid.jl")
export FVMPDEGrid

include("FVMPDE_problem.jl")
export FVMPDEProblem

include("FVMPDE_methods.jl")
export UpWind_F, UpWind_S, UpWind_Q, UpWind_C, SCD, Koren, SuperBee

include("FVMPDE_solver.jl")
export FVMPDESolve

end # module
