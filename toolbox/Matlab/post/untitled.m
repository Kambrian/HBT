global G HUBBLE0 a Omega0 OmegaLambda
G=43007.1;HUBBLE0=0.1;Omega0=0.3;OmegaLambda=0.7;
runnum={'6702DM'};
scaleF_file=['/mnt/A4700/data/',runnum{1},'/subcat/Redshift.dat'];
tmp=load(scaleF_file);a=tmp(:,2);
history_file=['/mnt/A4700/data/',runnum{1},'/subcat/anal/history_000_010.dat'];
history_par_file=['/mnt/A4700/data/',runnum{1},'/subcat/anal/historypar_000_010.dat'];
disp('loading files...')
[hist_bin,pmass]=load_hist_bin(history_file);
[Snappar,Mratepar,Mhostpar,Rhostpar,Kpar,j2par,Kerrpar,j2errpar,Chostpar,Csatpar]=load_hist_par(history_par_file);
%%
nnode=zeros(size(hist_bin));
ninfall=nnode;
msat=nnode;
snapend=nnode;
mend=nnode;
virsub=zeros(size(hist_bin));
rmin=zeros(size(hist_bin))-1;
factor_reacc=nnode+1;
for i=1:numel(hist_bin)
[Ninfall,t,r,v,vt2,kt,mdm,mhost,chost,Hz,Rvir,virialF,virsubs]=get_DMstrp_history(hist_bin{i},pmass,Snappar(i,2),1,1,1);
dmdm=diff(mdm);
imin=0;imax=0;
if Ninfall>1
for j=1:Ninfall-1
    if dmdm(j)>0
        imin=j;
        break;
    end
end
if imin>0
for j=imin:Ninfall-1
    if dmdm(j)<0
        imax=j;
        break;
    end
end
end
end
if imax>0
   factor_reacc(i)=mdm(imax)/mdm(imin);
end  
rmin(i)=min(r./Rvir);
nnode(i)=hist_bin{i}.Nnode;
nodeinfall=Snappar(i,3)-Snappar(i,1);
if(nodeinfall>0)
virsub(i)=hist_bin{i}.node(nodeinfall).mass(4)/hist_bin{i}.node(nodeinfall).mass(1);
ninfall(i)=nnode(i)-nodeinfall;
msat(i)=hist_bin{i}.node(nodeinfall).mass(1);
end
snapend(i)=hist_bin{i}.node(end).Nsnap;
mend(i)=hist_bin{i}.node(end).mass(1);
end
%%
find(ninfall>20&virsub>0&virsub<1.3&msat>5000&rmin<0.5&factor_reacc>1.4&factor_reacc<3)
figure;x=1:0.1:5;y=histc(factor_reacc(logical(factor_reacc)),x);plot(x,y,'o-');set(gca,'yscale','log');