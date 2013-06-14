% datadir='/mnt/A4700/data/8213/subcat/anal/massfun/';
datadir='/mnt/A4700/data/6113/subcat/anal/massfun/';
fid=fopen([datadir,'grptab_080'],'r');
grptab=fread(fid,'int32');
fclose(fid);
fid=fopen([datadir,'sublen_080'],'r');
sublen=fread(fid,'int32');
fclose(fid);

Ngroups=numel(grptab)/3;
grptab=reshape(grptab,3,Ngroups)';
Nsubs=numel(sublen);

n=histc(grptab(:,1),logspace(1,log10(grptab(1,1)),50));
figure;loglog(logspace(1,log10(grptab(1,1)),50),n);xlabel('Halo Mass (particles)');ylabel('Counts');title('halo mass fun');

subfrac=zeros(1,Nsubs);
for i=1:Ngroups
    if grptab(i,3)>1
    subind=grptab(i,2)+(2:grptab(i,3));%do not include main-sub
    subfrac(subind)=sublen(subind)/grptab(i,1);
    end
end

nbin=50;
G=43007.1;
Hz=0.1;Omega=0.3;
% partmass=0.023380;
% partmass=0.0620373;
partmass=0.187042;
massbin=[1,5,10,50,100,400]*100/partmass;
masslabel={'(1~4)x10^{14}Msun/h','(5~10)x10^{13}Msun/h','(1~5)x10^{13}Msun/h','(5~10)x10^{12}Msun/h','(1~5)x10^{12}Msun/h','refline:slope=-0.9 (-1.9)'};%[1e12,5e12,1e13,5e13,1e14,4e14]
grpbin=histc(grptab(:,1),massbin);
grpbin=grpbin(end:-1:1);
grpind=cumsum(grpbin)+1;
subbin=[grptab(grpind(1:end-1),2)+1,grptab(grpind(2:end),2)];
stylecyc={'r- .','g: .','b: .','c-- .','k-- .'}
figure;
for i=1:5
    subind=subbin(i,1):subbin(i,2);
    subfrac_bin=subfrac(subind);
    subfrac_bin=subfrac_bin(logical(subfrac_bin));
    xmin=min(subfrac_bin);
%     xmax=0.01;
    xmax=max(subfrac_bin)*1.001;
    xbin=logspace(log10(xmin),log10(xmax),nbin);
    ybin=histc(subfrac(subind),xbin)/grpbin(i+1);
    ycumbin=cumsum(ybin(end:-1:1));
    subplot(2,1,1);
    loglog(xbin(2:end),ycumbin(end:-1:2),stylecyc{i});
    hold on;
    subplot(2,1,2);
    loglog(xbin(2:end),ybin(1:end-1)./diff(xbin),stylecyc{i});
    hold on;
end
subplot(2,1,1);plot([10^-3.5,10^-1.5],[10^2.4,10^0.6],'--');hold off;
legend(masslabel);ylabel('<N(>M_{sub}/M_{host})> per host');title('6102,z=0.71,ips=2908,(non-main)subhalo massfun of fof halos')
subplot(2,1,2);plot([10^-3,10^-1],[10^5.8,10^2],'--');hold off;
xlabel('M_{sub}/M_{host}');ylabel('<dN_{sub}/dM_{sub}> per host');

