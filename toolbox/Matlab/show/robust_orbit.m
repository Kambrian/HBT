global G HUBBLE0 a Omega0 OmegaLambda
G=43007.1;HUBBLE0=0.1;Omega0=0.3;OmegaLambda=0.7;
runnum={'6702DM'};
scaleF_file=['/mnt/A4700/data/',runnum{1},'/subcat/Redshift.dat'];
tmp=load(scaleF_file);a=tmp(:,2);
history_file=['/mnt/A4700/data/',runnum{1},'/subcat/anal/history_000_050.dat'];
history_par_file=['/mnt/A4700/data/',runnum{1},'/subcat/anal/historypar_000_050.dat'];
disp('loading files...')
[hist_bin,pmass]=load_hist_bin(history_file);
[Snappar,Mratepar,Mhostpar,Rhostpar,Kpar,j2par,Kerrpar,j2errpar,Chostpar,Csatpar]=load_hist_par(history_par_file);
%%
% snapwant=43;subidwant=5562;
snapwant=51;subidwant=7990;
% snapwant=40;subidwant=53;
% snapwant=36;subidwant=5211;
Nhist=size(hist_bin,1);
for h=1:Nhist
    Nnode=hist_bin{h}.Nnode;
    nodeid=snapwant-hist_bin{h}.node(1).Nsnap+1;
    if nodeid>Nnode, continue, end
    if hist_bin{h}.node(nodeid).subid==subidwant
        h
        break
    end
end
%%

[Ninfall,t,r,v,vt2,kt,mdm,mhost,chost,Hz,Rvir,virialF,virsubs]=get_DMstrp_history(hist_bin{h},pmass,1,0,0,0);
figure;
x=hist_bin{h}.node(1).Nsnap+(0:Ninfall-1);
plot(x,mdm/max(mdm),'k-- o');
hold on;
plot(x,mhost/max(mhost),'r-. .');
plot([Snappar(h,2),Snappar(h,2)],[0.01,1]);
plot([Snappar(h,3),Snappar(h,3)],[0.01,1],'--');
plot(x,r./Rvir,'b- *');
set(gca,'yscale','log');
legend('DM','host','Dstrp','Merge','r');
