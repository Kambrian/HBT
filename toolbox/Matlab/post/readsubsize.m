clear;
runnum='8213';Nsnap=59;
datadir=['/mnt/A4700/data/',runnum,'/subcatS/profile/'];
basedir='logbin';
% basedir='';
sizefile=fullfile(datadir,basedir,['sat_size_',num2str(Nsnap,'%03d')]);
proffile=fullfile(datadir,basedir,['sat_prof_',num2str(Nsnap,'%03d')]);
fid=fopen(sizefile,'r');
% 	int nbin;
% 	int flag_badbins;//1 means badbin,-1 means bad host,0 normal
% 	float rmax;
% 	float rvir;
% 	float req_bk_1;
% 	float req_all_1;
% 	float req_bk_02;
% 	float req_all_02;
% 	float rtidal;
% 	float rcen;
nbin=fread(fid,'int32',9*4);
fseek(fid,1*4,'bof');
flag_badbins=fread(fid,'int32',9*4);
fseek(fid,2*4,'bof');
rmax=fread(fid,'float32',9*4);
fseek(fid,3*4,'bof');
rvir=fread(fid,'float32',9*4);
fseek(fid,4*4,'bof');
req_bk_1=fread(fid,'float32',9*4);
fseek(fid,5*4,'bof');
req_all_1=fread(fid,'float32',9*4);
fseek(fid,6*4,'bof');
req_bk_02=fread(fid,'float32',9*4);
fseek(fid,7*4,'bof');
req_all_02=fread(fid,'float32',9*4);
fseek(fid,8*4,'bof');
rtidal=fread(fid,'float32',9*4);
fseek(fid,9*4,'bof');
rcen=fread(fid,'float32',9*4);
fclose(fid);
nsubs=numel(nbin);

f=@(x) x.*(3*(log(1+x)-x./(1+x))-(x./(1+x)).^2).^(-1.0/3);%ritdal_sub/rvir_sub as a function of rcen_sub/rs_host
%%

figure;
set(gcf,'DefaultLineMarkerSize',3);
nplot=1:nsubs;
% nplot=2:632;%HBT,m>300
% nplot=2:572;%subfind,m>300
nplot=nplot(~flag_badbins(nplot)&logical(rcen(nplot)));%all the satellites
% nplot=nplot(~flag_badbins(nplot)&logical(rcen(nplot)<rvir(1)));
axes('position',[.1,.5,.4,.4]);loglog(rtidal(nplot),req_bk_1(nplot),'o');hold on;loglog([1,1e4],[1,1e4]);ylabel('req_{bk}-1');axis([1,1e4,1,1e4]);
axes('position',[.5,.5,.4,.4]);loglog(rtidal(nplot),req_all_1(nplot),'o');hold on;loglog([1,1e4],[1,1e4]);ylabel('req_{all}-1');set(gca,'yaxislocation','right');axis([1,1e4,1,1e4]);
axes('position',[.1,.1,.4,.4]);loglog(rtidal(nplot),req_bk_02(nplot),'o');hold on;loglog([1,1e4],[1,1e4]);ylabel('req_{bk}-02');axis([1,1e4,1,1e4]);xlabel('rtidal');
axes('position',[.5,.1,.4,.4]);loglog(rtidal(nplot),req_all_02(nplot),'o');hold on;loglog([1,1e4],[1,1e4]);ylabel('req_{all}-02');set(gca,'yaxislocation','right');axis([1,1e4,1,1e4]);
xlabel('rtidal');
% 
% figure;
% % nplot=1500:15333;
% % nplot=nplot(~flag_badbins(1:965)&logical(rcen(1:965)<rvir(1)));
% subplot(2,2,1);loglog(rmax(nplot),req_bk_1(nplot),'.');hold on;loglog([1,1e4],[1,1e4]);ylabel('req_{bk}-1');axis equal;
% subplot(2,2,2);loglog(rmax(nplot),req_all_1(nplot),'.');hold on;loglog([1,1e4],[1,1e4]);ylabel('req_{all}-1');axis equal;
% subplot(2,2,3);loglog(rmax(nplot),req_bk_02(nplot),'.');hold on;loglog([1,1e4],[1,1e4]);ylabel('req_{bk}-02');axis equal;
% subplot(2,2,4);loglog(rmax(nplot),req_all_02(nplot),'.');hold on;loglog([1,1e4],[1,1e4]);ylabel('req_{all}-02');axis equal;	
% xlabel('rmax');
% 
% figure;
% % nplot=1500:15333;
% % nplot=nplot(~flag_badbins(1:965)&logical(rcen(1:965)<rvir(1)));
% subplot(2,2,1);loglog(rvir(nplot),req_bk_1(nplot),'.');hold on;loglog([1,1e4],[1,1e4]);ylabel('req_{bk}-1');axis equal;
% subplot(2,2,2);loglog(rvir(nplot),req_all_1(nplot),'.');hold on;loglog([1,1e4],[1,1e4]);ylabel('req_{all}-1');axis equal;
% subplot(2,2,3);loglog(rvir(nplot),req_bk_02(nplot),'.');hold on;loglog([1,1e4],[1,1e4]);ylabel('req_{bk}-02');axis equal;
% subplot(2,2,4);loglog(rvir(nplot),req_all_02(nplot),'.');hold on;loglog([1,1e4],[1,1e4]);ylabel('req_{all}-02');axis equal;	
% xlabel('rvir');

%%
% for subind=1:15333;
subind=1;
%subid=subind-1;
fid2=fopen(proffile,'r');
nskip=sum(nbin(1:subind-1))*6*4;
fseek(fid2,nskip,'bof');
ns=fread(fid2,nbin(subind),'int32');
no=fread(fid2,nbin(subind),'int32');
nb=fread(fid2,nbin(subind),'int32');
rs=fread(fid2,nbin(subind),'float32');
ro=fread(fid2,nbin(subind),'float32');
rb=fread(fid2,nbin(subind),'float32');
fclose(fid2);
% if(sum(ns)<300&&nbin(subind)~=0)
%     break;
% end
% end
% if (nbin(subind)~=0)&&(sum(nb+no)==0)
%     [subind,req_all_1(subind),rvir(subind)]
% end
% end

if strcmp(basedir,'logbin')
vs=diff(logspace(log10(max(1.5,rmax(subind)*1e-2)),log10(rmax(subind)),nbin(subind)+1).^3)';
else
vs=diff((0:nbin(subind)).^3)';
end
svs=cumsum(vs);
na=ns+no+nb;
sna=cumsum(na);
sns=cumsum(ns);
figure; subplot(2,2,1);
loglog(rs/rtidal(subind),ns./vs,'r- .');
hold on;
loglog(ro/rtidal(subind),no./vs,'g: .');
loglog(rb/rtidal(subind),nb./vs,'b-- .');
% loglog(rs/rtidal(subind),(no+nb)./vs,'c-. .');
loglog(rs/rtidal(subind),na./vs,'k-. .');
hold off;
ylabel('density');
legend('sub','other','back','all');

subplot(2,2,2);
loglog(rs/rtidal(subind),ns./vs.*svs(end)/sum(nb+no),'r- .');grid on;ylabel('sub density/average background');
hold on;
loglog(rs/rtidal(subind),ns./vs.*svs(end)/sum(nb),'b- .');
legend('to all','to back');

subplot(2,2,3);
% figure(1);
semilogx(rs/rtidal(subind),ns./vs./(sns./svs),'r-.');hold on;semilogx(rs/rtidal(subind),na./vs./(sna./svs),'g:*');legend('sub','all');ylabel('\rho/<\rho>');hold off;
xlabel('r/rtidal');

% semilogx(rs/rtidal(subind),ns./vs,'r- .');