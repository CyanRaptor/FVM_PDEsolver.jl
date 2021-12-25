
function FVMPDESolve(prob::FVMPDEProblem,tspan;kwargs...)
    scheme = haskey(kwargs,:scheme) ? kwargs[:scheme] : UpWind_F

    U0 = copy(prob.U0)
    ODE_Function(u,p,t) = FVMPDE_∂u∂t(prob,u,t,scheme)

    ODE_Problem = ODEProblem(ODE_Function,U0,tspan);
    sol = DifferentialEquations.solve(ODE_Problem)

    if checkDim(prob.grid.dimension) == 1
        for i = 1:length(sol.u)
            sol.u[i] = sol.u[i][:,1,1,:]
        end
    elseif checkDim(prob.grid.dimension) == 2
        for i = 1:length(sol.u)
            sol.u[i] = sol.u[i][:,:,1,:]
        end
    end

    return sol
end

function FVMPDE_∂u∂t(prob::FVMPDEProblem,u,t,scheme)
    θ = prob.param
    indices = prob.grid.indices
    Δx = prob.grid.dX[1]
    ∂F∂x = FVMPDE_∂F∂x(scheme,prob.F,u,Δx,indices,1,t,θ)

    if prob.grid.dimension == _2D
        Δy = prob.grid.dY[1]
        ∂G∂y = FVMPDE_∂F∂x(scheme,prob.G,u,Δy,indices,2,t,θ)
        ∂H∂z = 0
    elseif prob.grid.dimension == _3D
        Δy = prob.grid.dY[1]
        Δz = prob.grid.dZ[1]
        ∂G∂y = FVMPDE_∂F∂x(scheme,prob.G,u,Δy,indices,2,t,θ)
        ∂H∂z = FVMPDE_∂F∂x(scheme,prob.H,u,Δz,indices,3,t,θ)
    else
        ∂G∂y = 0
        ∂H∂z = 0
    end

    S = prob.S(u)

    ∂u∂t = S .- (∂F∂x .+ ∂G∂y .+ ∂H∂z)

    return ∂u∂t
end

function FVMPDE_∂F∂x(scheme,func_F,u,Δx,indices,_dim,t,p)
    Uₓᴸ, Uₓᴿ = scheme(u,Δx,indices,_dim)

    ∂f∂x = zeros(size(u))

    for i in indices
        A⁺ , A⁻ = LRJacobian(func_F,u[i,:],t,p)
        ∂f∂x[i,:] = (A⁻ * Uₓᴸ[i,:] .+ A⁺ * Uₓᴿ[i,:])
    end
    return ∂f∂x
end


function LRJacobian(func_F,u,t,p)
    A = ForwardDiff.jacobian(func_F, u)
    #A = (abs.(A) .> 1e-10) .* A
    #println(A)
    Λ = diagm(eigvals(A))
    Λ⁺ = 0.5 .* (Λ .+ abs.(Λ))
    Λ⁻ = 0.5 .* (Λ .- abs.(Λ))
    P = eigvecs(A)
    P⁻¹ = inv(P)
    A⁺ = P * Λ⁺ * P⁻¹
    A⁻ = P * Λ⁻ * P⁻¹
    return A⁺, A⁻
end
