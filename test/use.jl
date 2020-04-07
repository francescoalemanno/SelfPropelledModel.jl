using SelfPropelledModel
using Test
using PyPlot

using DelimitedFiles

pygui(true)
xlim((-250,250))
ylim((-250,250))
sv=sample_particles(100,()->500.0.*(rand(2).-0.5),()->(rand(2).-0.5))
scatter(getindex.(sv.x,1),getindex.(sv.x,2),color="black",s=0.1)

for i in 1:170
    scatter(getindex.(sv.x,1),getindex.(sv.x,2),color="black",s=0.1)
    X=collect(foldl(hcat,[[x...] for x in sv.x])')
    writedlm("results/$i.txt",X)
    spm_step!(sv,0.4,0.2,0.2,1.5)
end

using SPTrackingToolkit

frames=[readdlm("results/$i.txt")'|>collect for i in 1:170]

specs = SPT(maxtimegap=10.0,maxdist=10.0,dims=2,rest=0,verbose=true)

tr=track(specs,frames)
for i in 1:100
    plot(tr[i][2,:],tr[i][3,:])
end
