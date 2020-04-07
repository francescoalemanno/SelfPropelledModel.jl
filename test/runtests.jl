using SelfPropelledModel
using Test
using PyPlot

pygui(true)
xlim((-250,250))
ylim((-250,250))
sv=sample_particles(100,()->500.0.*(rand(2).-0.5),()->(rand(2).-0.5))
scatter(getindex.(sv.x,1),getindex.(sv.x,2),color="black",s=0.1)

for i in 1:170
    scatter(getindex.(sv.x,1),getindex.(sv.x,2),color="black",s=0.1)
    spm_step!(sv,0.4,0.2,0.2,1.5)
end
