#include("HPDE_options.jl")
##########################################################
abstract type AbstractHPDEProblem end

abstract type Dimension end
abstract type _1D <: Dimension end
abstract type _2D <: Dimension end
abstract type _3D <: Dimension end

function checkDim(D::Type{<:Dimension})
    if D === _1D
        return 1
    elseif D === _2D
        return 2
    elseif D === _3D
        return 3
    end
    return 0
end

mutable struct FVMPDEProblem
    """
    Flux function F in x direction

    """
    F::Function

    """
    Flux function G in y direction
    !!! compat
        Will be added in future versions
    """
    G::Function

    """
    Flux function H in z direction
    !!! compat
        Will be added in future versions
    """
    H::Function

    """
    Source function
    !!! compat
        Will be added in future versions
    """
    S::Function

    nvars::Int8

    param::Any

    U0

    U_min
    U_max

    dimension::Type{<:Dimension}

    function FVMPDEProblem(func_F::Function, U0; kwargs...)

        zeroFunc(U) = zeros(size(U))
        F = func_F
        D = _1D
        (G,D) = haskey(kwargs,:func_G) ? (kwargs[:func_G],_2D) : (zeroFunc,D)
        (H,D) = haskey(kwargs,:func_G) && haskey(kwargs,:func_H) ? (kwargs[:func_H],_3D) : (zeroFunc,D)
        S = haskey(kwargs,:func_S) ? kwargs[:func_S] : zeroFunc

        @assert size(F(U0)) == size(U0)
        @assert size(G(U0)) == size(U0)
        @assert size(H(U0)) == size(U0)
        @assert size(S(U0)) == size(U0)

        θ = haskey(kwargs,:parameters) ? kwargs[:parameters] : nothing

        U_min = haskey(kwargs,:ul) ? kwargs[:ul] : -Inf * ones(size(U0))
        U_max = haskey(kwargs,:uh) ? kwargs[:uh] : Inf * ones(size(U0))

        @assert size(U_min) == size(U0)
        @assert size(U_max) == size(U0)

        nvars = size(U0,1)

        #@assert checkDim(D) == size(U0,)

        return new(F, G, H, S, nvars, θ, U0, U_min, U_max, D )
    end


end
