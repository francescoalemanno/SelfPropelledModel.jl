module SelfPropelledModel
using StructArrays
using Random

struct Particle{S}
    x::S
    v::S
    ox::S
    ov::S
end

function Particle(x::A,v::A,ox::A,ov::A) where {N,T,A<:NTuple{N,T}}
    Particle{A}(x,v,ox,ov)
end

function Particle(a_x,a_v)
    x=Tuple(a_x)
    v=Tuple(a_v)
    Particle(x,v,x,v)
end

function sqnorm(a)
    foldl((x,y)->x+y*y,a,init=zero(first(a)))
end

function θ(a,b)
    n1=a*a
    n2=b*b
    ns=sqrt(n1+n2)
    return a/ns,b/ns
end

function dist(v)
    N=length(v)
    d=zero(v.x[1][1])
    m=zero(v.x[1][1])
    id=0
    for i in 2:(N-1)
        nd=sqnorm(v[i].x.-v[i+1].x)+sqnorm(v[i].x.-v[i-1].x)
        if nd>=d
            id=i
            d=nd
        end
        m+=nd
    end
    d,id,m/N
end

function optimize!(v)
    for i in 1:(length(v)*length(v))
        D,j,mD=dist(v)
        r=rand(eachindex(v))
        v[j],v[r]=v[r],v[j]
        nD,nj,nmD=dist(v)
        if nmD>mD && nD>D
            v[j],v[r]=v[r],v[j]
        end
        r1=rand(eachindex(v))
        r2=rand(eachindex(v))
        v[j+1],v[r1]=v[r1],v[j+1]
        v[j-1],v[r2]=v[r2],v[j-1]
        nD,nj,nmD=dist(v)
        if nmD>mD && nD>D
            v[j-1],v[r2]=v[r2],v[j-1]
            v[j+1],v[r1]=v[r1],v[j+1]
        end
    end
end

function update_velocity(i::Int,P,α,γ,ρ)
    p=P[i]
    N=length(P)
    v0=sqrt(sqnorm(p.ov))
    p_n=P[ifelse(i<N,i+1,1)]
    p_b=P[ifelse(i>1,i-1,N)]
    η=randn(length(p.x))
    nvx=v0 .* θ(@.(α * p.ov + γ * (p_n.ov+p_b.ov)/2 + ρ*η)...)
    P[i]=Particle(p.x,nvx,p.ox,p.ov)
end

function update_position(i::Int,P,Δt)
    p=P[i]
    nx=p.x .+ p.v .* Δt
    P[i]=Particle(nx,p.v,p.x,p.v)
end

function spm_step!(P,α,γ,ρ,Δt)
    Threads.@threads for i in eachindex(P)
        update_velocity(i,P,α,γ,ρ)
    end
    Threads.@threads for i in eachindex(P)
        update_position(i,P,Δt)
    end
end

function sample_particles(N,x_sampler,v_sampler)
    sv=StructArray(Particle(x_sampler(),v_sampler()) for i in 1:N)
    optimize!(sv)
    optimize!(sv)
    sv
end

export sample_particles, spm_step!
end # module
