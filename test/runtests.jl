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

pygui(true)
matplotlibstyle()
Ns=20:20:160
A=[(N,readdlm("P$(N)_F170.dat")) for N in Ns]

avgA=foldl(hcat,[N,Statistics.mean(x,dims=1)...,std(x,dims=1)...] for (N,x) in A)

py"setfonts(12)"
figure(figsize=py"mplfigsize(1)".*1.8,dpi=300)
errorbar(avgA[1,:],avgA[2,:],yerr=avgA[4,:]/10,fmt="o",lw=2,color="royalblue",
    markersize=7,markeredgewidth=2,marker="o",markerfacecolor="None",
    markeredgecolor="royalblue")
xticks(Ns)
xlabel(L"N^\circ \rm{particles}")
ylabel(L"N^\circ \rm{erroneous\,\,tracks}")
tight_layout()
savefig("n_errors_vs_n_particles.pdf")

figure(figsize=py"mplfigsize(1)".*1.8,dpi=300)
errorbar(avgA[1,:],avgA[3,:],yerr=avgA[5,:]/10,fmt="o",lw=2,color="royalblue",
    markersize=7,markeredgewidth=2,marker="o",markerfacecolor="None",
    markeredgecolor="royalblue")
xticks(Ns)
xlabel(L"N^\circ \rm{particles}")
ylabel(L"\rm{Percentual\,\,errors\,\,}\%")
tight_layout()
savefig("rel_errors_vs_n_particles.pdf")



using LsqFit
using LinearAlgebra


wt = inv.(avgA[4,:]/10).^2
m(n, p) = @. p[1]*p[1]*n*n + p[2]*p[1]*n + p[3]
fit = curve_fit(m, avgA[1,:],avgA[2,:], wt, [0.01,0.01,0.01])
cov = estimate_covar(fit)
P=fit.param


diag(cov).^0.5

py"setfonts(12)"
figure(figsize=py"mplfigsize(1.4)".*1.8,dpi=300)
subplot(2,1,1)
errorbar(avgA[1,:],avgA[2,:],yerr=avgA[4,:]/10,fmt="o",lw=1,color="C0",
    markersize=6,markeredgewidth=1.333,marker="o",markerfacecolor="None",
    markeredgecolor="C0")
plot(0:180,m(0:180,P),lw=1.333)
xticks(Ns)
xlim(10,170)
xlabel(L"N^\circ \rm{particles}")
ylabel(L"N^\circ \rm{erroneous\,\,tracks}")


χ²=sum(abs2,(avgA[2,:]-m(avgA[1,:],P))./(avgA[4,:]/10))/(length(avgA[2,:]) -3 )

text(20, 9, "χ² = $(round(χ²*100)/100)", fontsize=12)

subplot(2,1,2)
scatter(avgA[1,:],(avgA[2,:]-m(avgA[1,:],P))./(avgA[4,:]/10))
xlabel(L"N^\circ \rm{particles}")
ylabel(L"\rm{Std. Residuals}")
xticks(Ns)
xlim(10,170)
ylim(-1,1)
tight_layout(rect=[0,0,0.95,0.95])
savefig("n_errors_vs_n_particles.pdf")






wt = inv.(avgA[5,:]/10).^2
m(n, p) = @. p[1]*p[1]*n*n + p[2]*p[1]*n + p[3]
fit = curve_fit(m, avgA[1,:],avgA[3,:], wt, [0.01,0.01,0.01])
cov = estimate_covar(fit)
P=fit.param

py"setfonts(12)"
figure(figsize=py"mplfigsize(1.4)".*1.8,dpi=300)
subplot(2,1,1)
errorbar(avgA[1,:],avgA[3,:],yerr=avgA[5,:]/10,fmt="o",lw=1,color="C0",
    markersize=6,markeredgewidth=1.333,marker="o",markerfacecolor="None",
    markeredgecolor="C0")
plot(0:180,m(0:180,P),lw=1.333)
xticks(Ns)
xlim(10,170)
xlabel(L"N^\circ \rm{particles}")
ylabel(L"\rm{Percentual\,\,errors\,\,}\%")

χ²=sum(abs2,(avgA[3,:]-m(avgA[1,:],P))./(avgA[5,:]/10))/(length(avgA[2,:]) -3 )

text(20, 0.04, "χ² = $(round(χ²*100)/100)", fontsize=12)

subplot(2,1,2)
scatter(avgA[1,:],(avgA[3,:]-m(avgA[1,:],P))./(avgA[5,:]/10))
xlabel(L"N^\circ \rm{particles}")
ylabel(L"\rm{Std. Residuals}")
xticks(Ns)
xlim(10,170)
ylim(-1,1)
tight_layout(rect=[0,0,0.95,0.95])

tight_layout()
savefig("rel_errors_vs_n_particles.pdf")




py"setfonts(12)"
figure(figsize=py"mplfigsize(1)".*1.4,dpi=300)



A=randn(10000).*2.2 .+ 16.3
B=randn(10000).*0.6 .+ 19.5
AB=[A;B]

samples=(2,2)
dhp=[abs(mean(rand(A,samples[1]))-mean(rand(B,samples[2]))) for i in 1:120000]
dnull=[abs(mean(rand(AB,samples[1]))-mean(rand(AB,samples[2]))) for i in 1:120000]
using Statistics
qnull=quantile(dnull,0.95)
V=[dhp dnull]


hist(V,60,density=true, histtype="step", stacked=false,label=["Truth", "Null"], fill=false,color=["C3","C0"])
vlines(qnull,0,0.29,linestyle="dashed")
text(qnull, 0.3, "95%", fontsize=10, horizontalalignment="center")
text(mean(B)-mean(A), 0.3, "obs", fontsize=10, horizontalalignment="center")
vlines(mean(B)-mean(A),0,0.29)
legend()
xlabel(L"ΔT = |Δτ_{mcf7} - Δτ_{mda}|")
ylabel(L"p.d.f.")
xlim(0,10)
savefig("timecomp.pdf")



sv=sample_particles(160,()->500.0.*(rand(2).-0.5),()->(rand(2).-0.5))
frames=[]

for i in 1:170
    X=foldl(hcat,[[p.x...] for p in sv])
    push!(frames,X)
    spm_step!(sv,0.4,0.2,0.2,1.5)
end


pygui(true)
figure(figsize=(6,6))
for f in frames
    scatter(f[1,:],f[2,:],color="black",s=0.125)
end
