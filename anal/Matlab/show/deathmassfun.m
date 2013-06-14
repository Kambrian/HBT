basedir='/mnt/A4700/data/6702/subcat/anal';
load([basedir,'/DeathMassfun_098']);
load([basedir,'/DeathMassfun_040']);

xmin=min(subinR);
xmax=max(subinR)*1.001;
x=logspace(log10(xmin),log10(xmax),nbin);
% subinR=sub_mass_099(find(Rsub>rmin*Rvir&Rsub<rmax*Rvir),1)/M0;
massfun2=histc(DeathMassfun_040(1:70,1)*partmass/M0,x);summass2=cumsum(massfun2(nbin:-1:1));

figure(1);%subplot(2,1,1);hold on;
loglog(x(2:nbin),summass2(nbin:-1:2),'-- oc');
