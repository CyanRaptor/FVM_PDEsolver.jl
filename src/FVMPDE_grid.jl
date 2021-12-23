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



mutable struct FVMPDEGrid
    X
    Y
    Z
    dX
    dY
    dZ
    U
    Uboundries
    dimension::Type{<:Dimension}
    function FVMPDEGrid(X,Y,Z,dX,dY,dZ,U,Uboundries,dimension)
        return new(X,Y,Z,dX,dY,dZ,U,Uboundries,dimension)
    end
end

function FVMPDEGrid(_X::AbstractVector, _Y::AbstractVector, _Z::AbstractVector)
    X = _X
    Y = _Y
    Z = _Z
    dX = _X[2:end] .- _X[1:end-1]
    dY = _Y[2:end] .- _Y[1:end-1]
    dZ = _Z[2:end] .- _Z[1:end-1]
    #Uboundries = zeros(length(_X),length(_Y))
    dimension = _2D

    return FVMPDEGrid(X,Y,Z,dX,dY,dZ,nothing,nothing,dimension)
end

function FVMPDEGrid(_X::AbstractVector, _Y::AbstractVector)
    X = _X
    Y = _Y
    dX = _X[2:end] .- _X[1:end-1]
    dY = _Y[2:end] .- _Y[1:end-1]
    #Uboundries = zeros(length(_X),length(_Y))
    dimension = _2D

    return FVMPDEGrid(X,Y,nothing,dX,dY,nothing,nothing,nothing,dimension)
end

function FVMPDEGrid(_X::AbstractVector)
    X = _X
    dX = _X[2:end] .- _X[1:end-1]
    #Uboundries = zeros(length(_X))
    dimension = _1D

    return FVMPDEGrid(X,nothing,nothing,dX,nothing,nothing,nothing,nothing,dimension)
end