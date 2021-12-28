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

    jac

    function FVMPDEProblem(grid::FVMPDEGrid, U0, func_F::Function; kwargs...)
        jac = haskey(kwargs,:jacobian) ? kwargs[:jacobian] : false
        zeroFunc(U) = zeros(size(U))
        F = func_F
        @assert size(U0,1) == grid.nx
        if haskey(kwargs,:func_G)
            G = kwargs[:func_G]
            @assert size(U0,2) == grid.ny
            if haskey(kwargs,:func_H)
                H = kwargs[:func_H]
                nvars = size(U0,4)
                @assert grid.dimension == _3D
                @assert length(size(U0))-1 == 3
                @assert size(U0,3) == grid.nz
            else
                H = zeroFunc
                nvars = size(U0,3)
                @assert grid.dimension == _2D
                @assert length(size(U0))-1 == 2
            end
        else
            G = zeroFunc
            H = zeroFunc
            nvars = size(U0,2)
            @assert grid.dimension == _1D
            @assert length(size(U0))-1 == 1
        end

        S = haskey(kwargs,:func_S) ? kwargs[:func_S] : zeroFunc
        if !jac
            @assert size(F(U0[grid.indices[1],:])) == size(U0[grid.indices[1],:])
            @assert size(G(U0[grid.indices[1],:])) == size(U0[grid.indices[1],:])
            @assert size(H(U0[grid.indices[1],:])) == size(U0[grid.indices[1],:])
        else
            @assert size(F(U0[grid.indices[1],:]),1) == length(U0[grid.indices[1],:])
            @assert size(F(U0[grid.indices[1],:]),2) == length(U0[grid.indices[1],:])
            # TODO: add @assert for G and H
        end
        @assert size(S(U0)) == size(U0)

        #grid.U = reshape(U0,grid.nx,grid.ny,grid.nz,nvars)
        #grid.U = copy(U0)
        θ = haskey(kwargs,:parameters) ? kwargs[:parameters] : nothing

        U_min = haskey(kwargs,:ul) ? kwargs[:ul] : -Inf * ones(size(U0))
        U_max = haskey(kwargs,:uh) ? kwargs[:uh] : Inf * ones(size(U0))

        @assert size(U_min) == size(U0)
        @assert size(U_max) == size(U0)



        #@assert checkDim(D) == size(U0,)

        return new(grid,F, G, H, S, nvars, θ, U0, U_min, U_max,jac)
    end


end
