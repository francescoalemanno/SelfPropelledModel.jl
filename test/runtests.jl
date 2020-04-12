using SelfPropelledModel
using Test
using SPTrackingToolkit
using DelimitedFiles

mean(X)=foldl((x,y)->(x[1]+1,(x[2]*x[1]+y)/(x[1]+1)),X,init=(0,zero(first(X))/1))[2]

function simulation(Np,Nf)
    sv=sample_particles(Np,()->500.0.*(rand(2).-0.5),()->(rand(2).-0.5))
    frames=[]

    for i in 1:Nf
        X=foldl(hcat,[[p.x...] for p in sv])
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
Ns=[20,40,60,80,100,120,140,160]
for p in Ns
    tohdsim(p,170)
end

using DelimitedFiles
using PyPlot
using PyCall
using Statistics
cd(dirname(@__FILE__));
include("mplstyle.jl")
cd("../results")

pygui(false)
matplotlibstyle()
Ns=20:20:160
A=[(N,readdlm("P$(N)_F170.dat")) for N in Ns]

avgA=foldl(hcat,[N,Statistics.mean(x,dims=1)...,std(x,dims=1)...] for (N,x) in A)

py"setfonts(14)"
figure(figsize=py"mplfigsize(1)".*1.4,dpi=300)
errorbar(avgA[1,:],avgA[2,:],yerr=avgA[4,:]/10,fmt="o",lw=2,color="royalblue",
    markersize=7,markeredgewidth=2,marker="o",markerfacecolor="None",
    markeredgecolor="royalblue")
xticks(Ns)
xlabel(L"N^\circ \rm{particles}")
ylabel(L"N^\circ \rm{erroneous\,\,tracks}")
tight_layout()
savefig("n_errors_vs_n_particles.pdf")

figure(figsize=py"mplfigsize(1)".*1.4,dpi=300)
errorbar(avgA[1,:],avgA[3,:],yerr=avgA[5,:]/10,fmt="o",lw=2,color="royalblue",
    markersize=7,markeredgewidth=2,marker="o",markerfacecolor="None",
    markeredgecolor="royalblue")
xticks(Ns)
xlabel(L"N^\circ \rm{particles}")
ylabel(L"\rm{Percentual\,\,errors\,\,}\%")
tight_layout()
savefig("rel_errors_vs_n_particles.pdf")
