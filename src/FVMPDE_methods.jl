function GeneralFluxLimiter(func_Ψ,u,Δx,indices,dim)
    ϵ = 1.0e-8

    func_r(Uₓ₁,Uₓ₂) = (Uₓ₁ .+ ϵ) ./ (Uₓ₂ .+ ϵ)
    func_U(Δxᵢ,r,uₓ) = (Δxᵢ / 2) .* func_Ψ(r) .* (uₓ .+ ϵ)

    d = length(size(u)) - 1
    #nx,ny,nz,nvars = size(u)

    minIndex = fill(Tuple(ones(Int8,d)),size(indices))
    maxIndex = fill(Tuple(size(u)[1:end-1]),size(indices))

    adderArray = zeros(Int8,d)
    adderArray[dim] = 1
    adder = CartesianIndex(Tuple(adderArray))

    ViewUᵢ₋₂ = CartesianIndex.(max.(Tuple.(indices .- adder .- adder),minIndex))
    ViewUᵢ₋₁ = CartesianIndex.(max.(Tuple.(indices .- adder         ),minIndex))
    ViewUᵢ₊₁ = CartesianIndex.(min.(Tuple.(indices .+ adder         ),maxIndex))
    ViewUᵢ₊₂ = CartesianIndex.(min.(Tuple.(indices .+ adder .+ adder),maxIndex))



    Uₓi₊¾ = (u[ViewUᵢ₊₂,:] .- u[ViewUᵢ₊₁,:]) ./ Δx    # ¾ : 3/2
    Uₓi₊½ = (u[ViewUᵢ₊₁,:] .- u) ./ Δx
    Uₓi₋½ = (u .- u[ViewUᵢ₋₁,:]) ./ Δx
    Uₓi₋¾ = (u[ViewUᵢ₋₁,:] .- u[ViewUᵢ₋₂,:]) ./ Δx

    Uᴸi₊½ = u[ViewUᵢ₊₁,:] .- func_U(Δx,func_r(Uₓi₊½,Uₓi₊¾),Uₓi₊¾)
    Uᴸi₋½ = u             .- func_U(Δx,func_r(Uₓi₋½,Uₓi₊½),Uₓi₊½)
    Uᴿi₊½ = u             .+ func_U(Δx,func_r(Uₓi₊½,Uₓi₋½),Uₓi₋½)
    Uᴿi₋½ = u[ViewUᵢ₋₁,:] .+ func_U(Δx,func_r(Uₓi₋½,Uₓi₋¾),Uₓi₋¾)

    Uₓᴸ = (Uᴸi₊½ .- Uᴸi₋½) ./ Δx
    Uₓᴿ = (Uᴿi₊½ .- Uᴿi₋½) ./ Δx
    #return Uᴿi₊½, Uᴿi₋½, Uᴸi₊½ , Uᴸi₋½
    #return Uₓi₊¾, Uₓi₊½, Uₓi₋½, Uₓi₋¾
    return Uₓᴸ, Uₓᴿ
end





################## LIMITERS ####################

function κ_Scheme(r,κ)
    α₁ = (1+κ)/2
    α₂ = (1-κ)/2
    return (α₁ .*  r .+ α₂)
end


function UpWind_F(u,Δx,indices,dim)
    f(r) = zeros(size(r))
    return GeneralFluxLimiter(f,u,Δx,indices,dim)
end

function UpWind_S(u,Δx,indices,dim)
    f(r) = κ_Scheme(r,-1)
    return GeneralFluxLimiter(f,u,Δx,indices,dim)
end

function SCD(u,Δx,indices,dim)
    f(r) = κ_Scheme(r,1)
    return GeneralFluxLimiter(f,u,Δx,indices,dim)
end

function UpWind_Q(u,Δx,indices,dim)
    f(r) = κ_Scheme(r,0.5)
    return GeneralFluxLimiter(f,u,Δx,indices,dim)
end

function UpWind_C(u,Δx,indices,dim)
    f(r) = κ_Scheme(r,1/3.0)
    return GeneralFluxLimiter(f,u,Δx,indices,dim)
end

function Koren(u,Δx,indices,dim)
    f(r) =  max.(0,min.(2 .* r,min.(1/3 .+ 2 /3 .* r,2)))
    return GeneralFluxLimiter(f,u,Δx,indices,dim)
end

function SuperBee(u,Δx,indices,dim)
    f(r) = max.(0,min.(r,1))
    return GeneralFluxLimiter(f,u,Δx,indices,dim)
end