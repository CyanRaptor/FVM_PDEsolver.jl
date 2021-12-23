#include("HPDE_options.jl")
##########################################################

mutable struct FVMPDEProblem
    grid::FVMPDEGrid

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

    function FVMPDEProblem(grid::FVMPDEGrid, U0, func_F::Function; kwargs...)

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
        @assert D == grid.dimension
        @assert length(size(U0))-1 == checkDim(D)
        @assert size(U0,1) == length(grid.X)
        if D == _2D || D == _3D
            @assert size(U0,2) == length(grid.Y)
        end
        if D == _3D
            @assert size(U0,3) == length(grid.Z)
        end
        grid.U = copy(U0)
        θ = haskey(kwargs,:parameters) ? kwargs[:parameters] : nothing

        U_min = haskey(kwargs,:ul) ? kwargs[:ul] : -Inf * ones(size(U0))
        U_max = haskey(kwargs,:uh) ? kwargs[:uh] : Inf * ones(size(U0))

        @assert size(U_min) == size(U0)
        @assert size(U_max) == size(U0)

        nvars = size(U0,checkDim(D)+1)

        #@assert checkDim(D) == size(U0,)

        return new(grid,F, G, H, S, nvars, θ, U0, U_min, U_max, D )
    end


end

X = [0.0, 0.1, 0.3, 0.6, 1.0, 1.5]
grid = FVMPDEGrid(X,X)

myF(U) = zeros(size(U))
u0=zeros(6,6,1)

FVMPDEProblem(grid,u0,myF,func_G=myF)
