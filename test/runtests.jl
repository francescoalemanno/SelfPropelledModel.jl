using SelfPropelledModel
using Test
using SPTrackingToolkit
using DelimitedFiles

mean(X)=foldl((x,y)->(x[1]+1,(x[2]*x[1]+y)/(x[1]+1)),X,init=(0,zero(first(X))/1))[2]

function simulation(Np,Nf)
    sv=sample_particles(Np,()->500.0.*(rand(2).-0.5),()->(rand(2).-0.5))
    frames=[]

    for i in 1:Nf
        X=foldl(hcat,[[x...] for x in sv.x])
        push!(frames,X)
        spm_step!(sv,0.4,0.2,0.2,1.5)
    end


    specs = SPT(maxtimegap=5.0,maxdist=5.0,dims=2,rest=0,verbose=true)

    tr=track(specs,frames)

    errors=zeros(length(tr))
    Threads.@threads for j in eachindex(tr)
        trace = tr[j]
        _ , Ntf=size(trace)
        op=nothing
        toterrors=-1
        for i in 1:Ntf
            t=trunc(Int,trace[1,i])
            X=trace[2:3,i]
            p=findfirst(isequal(X),eachcol(frames[t])|>collect)
            p==op || (toterrors+=1)
            op=p
        end
        errors[j]=toterrors/Ntf*100
    end

    (errors .> 0)|>sum,mean(errors)
end

function tohdsim(Np,Nf)
    R=[]
    for i in 1:300
        push!(R,simulation(Np,Nf))
    end
    AR=foldl(hcat,[x...] for x in R)'
    writedlm("results/P$(Np)_F$(Nf).dat",AR)
end
for p in [20,40,60,80,100,120,140,160]
    tohdsim(p,170)
end
