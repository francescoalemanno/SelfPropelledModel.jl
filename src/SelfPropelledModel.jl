module SelfPropelledModel
using Random

mutable struct Particle{S}
    d::Int
    class::Int
    x::S
    v::S
    ox::S
    ov::S
    function Particle(x::A,v::A,ox::A,ov::A) where {A<:AbstractVector}
        d1=length(x)
        d2=length(v)
        d3=length(ox)
        d4=length(ov)
        (d1==d2 && d3==d4 && d2==d3) || error("all arguments must have the same length.")
        new{A}(d1,0,copy(x),copy(v),copy(ox),copy(ov))
    end
end

Particle(x,v)=Particle(x,v,x,v)

mutable struct Params{T}
    Ci::T
    Ca::T
    Cb::T
    w::T
    γ::T
    δt::T
end

function params(;Ci=0.01,Ca=0.01,Cb=0.2,w=-2.0,γ=0.5,δt=0.5)
    tpl=promote(Ci,Ca,Cb,w,γ,δt)
    T=eltype(tpl)
    Params{T}(tpl...)
end


using LinearAlgebra
import Base.rand
import Random.rand!
Base.rand(::Type{Particle},d::Integer)=Particle(rand(d)*1000,randn(d)./sqrt(d))
Base.rand(::Type{Particle})=rand(Particle,2)
function Random.rand!(P::Particle)
    randn!(P.v)
    rand!(P.x)
    P.x.*=1000
    P.v./=sqrt(P.d)
    P.ox.=P.x
    P.ov.=P.v
    P
end

function update_v!(params::Params,p::Particle,nc::AbstractVector{<:Particle})
    p.v .= 0
    ew=exp(-params.w)
    for q in eachindex(nc)
        p.x .= p.ox .- nc[q].ox
        w = 1 / (1 + ew*abs2(norm(p.x)))
        p.v .+= w .* ((nc[q].ov .- p.ov) .* params.Ci .+ (nc[q].ox .- p.ox) .* params.Ca)
    end
    sc=norm(p.v)
    on=norm(p.ov)
    randn!(p.x)
    p.v .= p.ov .+ (params.δt * params.γ/sc) .* p.v  .+ params.Cb .* p.x
    nn=norm(p.v)
    p.v.*=on/nn
end

function update_x!(params::Params,p::Particle)
    p.x .= p.ox .+ params.δt .* p.ov
    p.ov .= p.v
    p.ox .= p.x
end

function step!(params::Params,P::AbstractVector{<:Particle})
    for i in eachindex(P)
        update_v!(params,P[i],P)
    end
    for i in eachindex(P)
        update_x!(params,P[i])
    end
end

function getX(P::AbstractVector{<:Particle})
    collect(foldl(hcat,getfield(p,:x) for p in P)')
end

par=params()
par.Ca=0.05
par.Ci=-0.05
par.Cb=0.3

V=[rand(Particle) for i in 1:100]
step!(par,V)

function steptest(par,V)
    for i in 1:2000
        step!(par,V)
        f=getX(V)
        scatter(f[:,1],f[:,2],color="blue",s=2)
        sleep(0.1)
    end
end

using PyPlot
pygui(true)
@time steptest(par,V)


export sample_particles, spm_step!
end # module
