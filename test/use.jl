using SelfPropelledModel
using Test
using PyPlot
include(dirname(@__FILE__)*"/mplstyle.jl")


N=160
sv=sample_particles(N,()->500.0.*(rand(2).-0.5),()->(rand(2).-0.5))

frames=[]
for i in 1:170
    X=foldl(hcat,[[x...] for x in sv.x])
    push!(frames,X)
    spm_step!(sv,0.4,0.2,0.2,1.5)
end

pts=PermutedDimsArray(cat(frames...,dims=3),(2,3,1))|>collect
pygui(false)
figure(figsize=(5,5))
for i in 1:N
    plot(pts[i,:,1],pts[i,:,2],color="black",lw=0.5)
end
xlim((-250,250))
ylim((-250,250))
tight_layout()
title("N=$(N) particles")
savefig("N$(N).pdf")
