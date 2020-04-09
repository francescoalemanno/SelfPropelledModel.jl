
D=zeros(0,4)
c=0
Nps=20:20:160
for i=Nps
    c=c+1
    M=dlmread(sprintf('P%d_F170.dat', i))
    avg1=mean(M(:,1))
    dev1=std(M(:,1))
    avg2=mean(M(:,2))
    dev2=std(M(:,2))
    D(c,1)=avg1
    D(c,2)=avg2
    D(c,3)=dev1
    D(c,4)=dev2
end

nerrors=D(:,1)
relerrors=D(:,2)
dev_nerrors=D(:,3)
dev_relerrors=D(:,4)

figure();
errorbar(Nps,nerrors,dev_nerrors/10,'Marker','o','LineStyle','none');
xlabel({'N° particles'});
ylabel({'N° errors'});
set(gcf,'color','w');
% Uncomment the following line to preserve the X-limits of the axes
xlim([16.5613490152539 164.357387784744]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim([-0.502182939714394 8.59480702684078]);
export_fig n_errors_vs_n_particles_matlab.pdf


figure();
errorbar(Nps,relerrors,dev_relerrors/10,'Marker','o','LineStyle','none');
xlabel({'N° particles'});
ylabel({'% errors'});
set(gcf,'color','w');
xlim([16.5613490152539 164.357387784744]);
export_fig rel_errors_vs_n_particles_matlab.pdf
