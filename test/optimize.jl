using SelfPropelledModel
using NaturalES
using LinearAlgebra

transpv(v)=collect(v')

function mysim(P)
    N=160
    sv=sample_particles(N,()->500.0.*(rand(2).-0.5),()->(rand(2).-0.5))
    frames=[]
    for i in 1:170
        X=foldl(hcat,[[p.x...] for p in sv])
        push!(frames,X)
        spm_step!(sv,P[1],P[2],P[3],P[4])
    end
    transpv.(frames)
end

function metrics(f,nf)
    Δ=nf.-f
    Np,Nc=size(Δ)
    m0=0.0
    m1=0.0
    m2=0.0
    m3=0.0
    m4=0.0
    for i in 3:Np-2, j in i-2:i+2
        s=normalize(Δ[i,:])⋅normalize(Δ[j,:])
        m0+=1
        m1+=s
        m2+=s*s
        m3+=s*s*s
        m4+=s*s*s*s
    end
    M1=m1/m0
    M2=m2/m0
    M3=m3/m0
    M4=m4/m0


    (M1,
    M2-M1*M1,
    M3-3M1*M2+2M1*M1*M1
    )
end

function metrics(frames)
    m=metrics(frames[1],frames[2])
    N=1
    for i in 2:length(frames)-1
        m=m.+metrics(frames[i],frames[i+1])
        N+=1
    end
    m./N
end
mysol=zeros(4)
Nmysol=0
function cost(x)
    global mysol.+=x
    global Nmysol+=1
    target=(0.20677718821137364, 0.5559917303075799, -0.14723289814743612)
    s=sqrt(sum(abs2,metrics(mysim(x)).-target))
    println(s," -> ",mysol./Nmysol)
    s
end

A=mysim([0.4,0.5,1.0,1.5])
metrics(A)

separable_nes(cost,[0.5,0.5,0.5,0.5],0.1)

x=mysol./mysol[3]
x[4]=1.5

frames=mysim(x)
using PyPlot
pts=PermutedDimsArray(cat(frames...,dims=3),(1,3,2))|>collect
pygui(true)
figure(figsize=(5,5))
for i in 1:160
    plot(pts[i,:,1],pts[i,:,2],color="black",lw=0.5)
end
xlim((-250,250))
ylim((-250,250))
tight_layout()
