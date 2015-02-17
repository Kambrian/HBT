clear
% global halo_file history_file history_par_file scaleF_file 
global a Omega0 OmegaLambda
Omega0=0.3;OmegaLambda=0.7;
% pmass gmass G HUBBLE0 Omega0 OmegaLambda

runnum={'6702DM'};pmass=zeros(size(runnum));
% runnum=[6700,6701,6702,6500,6501,6506,6402,6404,6409,6600,6601,6602];
% pmass=[0.0130057,0.00834475,0.008848,7.44348e-05,5.8378e-5,6.45321e-05,...
%     5.45271e-6,5.29086e-06,5.23311e-06,0.000789798,0.000670876,0.000546019];


G=43007.1;HUBBLE0=0.1;Omega0=0.3;OmegaLambda=0.7;

ADMmax=[];AKtmax=[];
xbin=0:0.02:1.2;
kbin=logspace(-3,2.8,100);
AcountM=zeros(numel(xbin),6);
AcountKt=zeros(numel(kbin),1);
for i=1:numel(runnum)
scaleF_file=['/mnt/A4700/data/',runnum{i},'/subcat/Redshift.dat'];
tmp=load(scaleF_file);a=tmp(:,2);    
history_file=['/mnt/A4700/data/',runnum{i},'/subcat/anal/history_000_010.dat'];
history_par_file=['/mnt/A4700/data/',runnum{i},'/subcat/anal/historypar_000_010.dat'];

disp('loading files...')
[hist_bin,pmass(i)]=load_hist_bin(history_file);
[Snappar,Mratepar,Mhostpar,Rhostpar,Kpar,j2par,Kerrpar,j2errpar,Chostpar,Csatpar]=load_hist_par(history_par_file);
% clear Mhostpar Rhostpar Kpar j2par Kerrpar j2errpar Chostpar Csatpar

disp('mass statistics...')
[countM,countKt,Ktmax]=mass_start_stat(hist_bin,Snappar,pmass(i),xbin,kbin);
AcountM=AcountM+countM;
AcountKt=AcountKt+countKt;
AKtmax=[AKtmax;Ktmax];
disp('Distance statistics...')
DMmax=Dist_Mmax(hist_bin,Snappar,pmass(i),Mratepar);
ADMmax=[ADMmax;DMmax];
% clear hist_bin Snappar
clear countM countKt DMmax Ktmax
end
plot_DMmax(ADMmax,AKtmax);
plot_Mstart_stat(AcountM,xbin,AcountKt,kbin);


%% Discussion: effect of local accretion? ommitting this would yield
%% overestimation of D_{Mmax}
%% what about the position when subhalo first have SubRank>0? what's the
%% relation of this distance to R_{vir}?