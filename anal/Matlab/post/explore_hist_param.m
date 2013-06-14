% clear;
% 
% halo=cell(100,1);
% for nsnap=0:99
% halo{nsnap+1}=readhalo_size('/mnt/A4700/data/6702/subcat/profile/logbin',nsnap,'halo');
% end
% save halo.mat halo
% load halo.mat
pmass=0.008848;gmass=0.00156133;
%%
hist_bin=load_hist_bin('/mnt/A4700/data/6702/subcat/anal/history_001_010.dat');
a=load('/home/kam/Documents/research/Galaxy/code/BoundTracing/v7.9/anal/Matlab/outputs_zoom.txt');
age=a;
for i=1:100
    age(i)=cosmo_time(0,a(i),0.3,0.7,0.1);
end
Nhist=size(hist_bin,1);
[Snappar,Mratepar,Mhostpar,Rhostpar,Kpar,j2par,Kerrpar,j2errpar]=load_hist_par('/mnt/A4700/data/6702/subcat/anal/historypar_000_010.dat');
%%
G=43007.1;HUBBLE0=0.1;Omega0=0.3;OmegaLambda=0.7;
 Hz=HUBBLE0 * sqrt(Omega0 /scaleF^3+ (1 -Omega0 -OmegaLambda) / scaleF^2 +OmegaLambda);
%%
orbpar=cell(100,1);
for nsnap=0:99
fid=fopen(['/mnt/A4700/data/6702/subcat/anal/orbpar/orbparam_',num2str(nsnap,'%03d')]);
orbpar{nsnap+1}=fread(fid,[5,inf],'float32');
fclose(fid);
end
%%
% j2=j.^2;
% x=logspace(-8,0,20);xm=(x(1:end-1)+x(2:end))/2;
% y=histc(j2(abs(j2)<1),-x(end:-1:1));y=y(1:end-1)'./diff(-x(end:-1:1));
% loglog(xm(end:-1:1),y,'- .');
filter=(1:Nhist)';
j2=j2par(:,2);
K=Kpar(:,2);
j=sqrt(j2);
% filter=logical(Nodeinfall>0&Ms>50);
figure;loglog(sqrt(1-j2),K,'.');xlabel('e');ylabel('\kappa=-2K_{vir}/U_{vir}');
figure;semilogx(sqrt(1-j2),(j2./(2-K)-1)./sqrt(1-j2),'.');xlabel('e');ylabel('cos(f)');
E=2-K;
x=logspace(-3,2.8,100);xm=(x(1:end-1)+x(2:end))/2;
y=histc(K(filter),x);y=y(1:end-1)'./diff(x);
% [ymax,imax]=max(y);xfit=xm(imax:end);yfit=y(imax:end);
xfit=xm(logical(xm>2));yfit=y(logical(xm>2));
xfit=xfit(find(yfit));yfit=yfit(find(yfit));
p=polyfit(log(xfit),log(yfit),1);
figure;axes('position',[.1,.7,.8,.26]);loglog(xm,y,'- .');hold on;
loglog(xm,exp(p(2))*xm.^p(1),'--');legend('data',['dN/d\kappa=',num2str(exp(p(2)),'%.1f'),'\kappa^{',num2str(p(1),'%.2f'),'}']);
xlabel('Energy Ratio \kappa=-2K_{vir}/U_{vir}=2-\epsilon=2-R_{vir}/a');ylabel('dN/d\kappa');
% e2=1-j.^2;
% x=logspace(-8,12,100);xm=(x(1:end-1)+x(2:end))/2;
% y=histc(e2(filter),x);y=y(1:end-1)'./diff(x);
% axes('position',[.1,.38,.8,.26]);loglog(xm,y,'- .');xlabel('Eccentricity^2 e^2=(c/a)^2=1-j^2');ylabel('dN/de^2');
e=sqrt(1-j2);
x=0:0.05:2;xm=(x(1:end-1)+x(2:end))/2;
y=histc(e(filter),x);y=y/(sum(y)+numel(find(e>2)));y=y(1:end-1)';y=y./diff(x);
axes('position',[.1,.38,.8,.26]);plot(xm,y,'- .');xlabel('Eccentricity e=c/a=sqrt(1-j^2)');ylabel('df/de');bar(xm,y);
x=0:0.05:1;xm=(x(1:end-1)+x(2:end))/2;
y=histc(j((j2>0)&filter),x);y=y(1:end-1)';y=y/sum(y);y=y./diff(x);
axes('position',[.1,.05,.8,.26]);plot(xm,y,'- .');xlabel('Circularity j=J/Jcir(E)=sqrt(1-e^2)');ylabel('df/dj');
hold on;%stairs(x,[y,0]);
loglog(xm,2.77*xm.^1.19.*(1.55-xm).^2.99,'--');
% j2=j2par(:,2)+j2errpar(:,2);j=sqrt(j2);y=histc(j((j2>0)&filter),x);y=y(1:end-1)';y=y/sum(y);y=y./diff(x);plot(xm,y,'.');
% j2=j2par(:,2)-j2errpar(:,2);j=sqrt(j2);y=histc(j((j2>0)&filter),x);y=y(1:end-1)';y=y/sum(y);y=y./diff(x);plot(xm,y,'.');
j2=j2errpar(:,2);j=sqrt(j2);y=histc(j((j2>0)&filter),x);y=y(1:end-1)';y=y/sum(y);y=y./diff(x);plot(xm,y,'.');
% hist(j(logical(j<1.1)&(~imag(j))&filter),0:0.05:1);xlabel('Circularity j=J/Jcir(E)=sqrt(1-e^2)');ylabel('N');
% figure;subplot(2,1,1);hist(log10(Mh));subplot(2,1,2);hist(log10(Mc));
% figure;loglog(Mi,K,'.');
% figure;subplot(2,1,1);loglog(Mi,K,'.','Markersize',3);xlabel('Msat');ylabel('\kappa');subplot(2,1,2);loglog(Mi./Mh,K,'.','markersize',3);xlabel('Msat/Mhost');ylabel('\kappa');
% figure;loglog(Mh,Mi./Mh,'.');
j2=j2par(:,2);
e2=1-j2;
x=logspace(-3,1,10);xm=(x(1:end-1)+x(2:end))/2;
y=histc(e2(filter),x);y=y(1:end-1)';%y=y/sum(y);y=y./diff(x);
figure;subplot(3,1,1);loglog(xm,y,'- .');xlabel('e2');ylabel('dN');
x=0:0.2:5;xm=(x(1:end-1)+x(2:end))/2;
y=histc(e2(filter),x);y=y(1:end-1)';%y=y/sum(y);y=y./diff(x);
subplot(3,1,2);plot(xm,y,'- .');xlabel('e2');ylabel('dN');
x=0:0.2:5;xm=(x(1:end-1)+x(2:end))/2;
y=histc(sqrt(e2(filter)),x);y=y(1:end-1)';%y=y/sum(y);y=y./diff(x);
subplot(3,1,3);plot(xm,y,'- .');xlabel('e');ylabel('dN');
%%
Dstrp=@(x) 3^(1/3)+x.^(1/3)/3+(2*x.^(2/3))/(9*3^(1/3))+x/(9*(3^(1/3))^2);
Dcomp=@(x) 3/2*3^(1/3)+x.^(1/3)/3+x.^(2/3)/3^(4/3)-7/9*x./3^(2/3);
Dmean=@(x) 5/4*3^(1/3)+x.^(1/3)/3+5/18*x.^(2/3)/3^(1/3)-x./3^(5/3);
   rgasmax=zeros(Nhist,1);
   rdmmax=zeros(Nhist,1);
   ktdm=zeros(Nhist,1);
   ktgas=zeros(Nhist,1);
   d_strp=zeros(Nhist,1);
   d_comp=zeros(Nhist,1);
   mrate=zeros(Nhist,1);
for h=1:Nhist
    Nnode=hist_bin{h}.Nnode;
    mdm=zeros(Nnode,1);
    mgas=zeros(Nnode,1);
    mhost=zeros(Nnode,1);
    Node_infall=-1;
    for s=1:Nnode
        mdm(s)=hist_bin{h}.node(s).mass(1);
        mgas(s)=hist_bin{h}.node(s).mass(2);
        mhost(s)=hist_bin{h}.node(s).mass(3);
        if hist_bin{h}.node(s).Nsnap==Snappar(h,3)
            if mdm(s)>500&&mdm(s)/mhost(s)<0.05
            Node_infall=s;
            end
        end
    end
    if Node_infall>0
%     nodestart=max(1,Nodebefore(h));
% nodestart=1;
% nodeend=Snappar(h,3)-Snappar(h,1)+1;
    [tmp,idm]=max(mdm(1:Node_infall));
    [tmp,igas]=max(mgas(1:Node_infall));
%     idm=idm+nodestart-1;igas=igas+nodestart-1;
    rgasmax(h)=norm(hist_bin{h}.node(igas).pos)/halo{hist_bin{h}.node(igas).Nsnap+1}.Rvir(hist_bin{h}.node(igas).HostID+1,1);
    rdmmax(h)=norm(hist_bin{h}.node(idm).pos)/halo{hist_bin{h}.node(idm).Nsnap+1}.Rvir(hist_bin{h}.node(idm).HostID+1,1);
    d_strp(h)=Dstrp(mdm(idm)/mhost(idm));
    d_comp(h)=Dcomp(mdm(idm)/mhost(idm));
    mrate(h)=mdm(idm)/mhost(idm);
    
    Vel=hist_bin{h}.node(idm).vel;
    Pos=hist_bin{h}.node(idm).pos*a(hist_bin{h}.node(idm).Nsnap+1);%physical
    r=sqrt(sum(Pos.^2,1));
    v=sqrt(sum(Vel.^2,1));
    mu=G*mhost(idm)*pmass;
    ktdm(h)=v^2*(1-(Pos'*Vel/r/v)^2)/(mu/r);
    Vel=hist_bin{h}.node(igas).vel;
    Pos=hist_bin{h}.node(igas).pos*a(hist_bin{h}.node(igas).Nsnap+1);%physical
    r=sqrt(sum(Pos.^2,1));
    v=sqrt(sum(Vel.^2,1));
    mu=G*mhost(igas)*pmass;
    ktgas(h)=v^2*(1-(Pos'*Vel/r/v)^2)/(mu/r);
    else
        rgasmax(h)=-1;
        rdmmax(h)=-1;
        ktdm(h)=-2;
        ktgas(h)=-2;
    end
end
rgasmax_theory=(2+ktgas).^(1/3)*1.25;
% rdmmax_theory=(2+ktdm).^(1/3);
% rdmmax_theory=d_strp;
% rgasmax_theory=d_comp;
rdmmax_theory=Diststrp;
% tidalrategas=rgasmax./rgasmax_theory;
tidalrategas=rgasmax./rdmmax_theory;
tidalratedm=rdmmax./rdmmax_theory;
x=0:0.2:10;
y1=histc(rdmmax(rdmmax>0),x);
y2=histc(rgasmax(rgasmax>0),x);
y11=histc(rdmmax_theory(rdmmax>0),x);
y22=histc(rgasmax_theory(rgasmax>0),x);
y3=histc(tidalratedm(rdmmax>0),x);
y4=histc(tidalrategas(rgasmax>0),x);
figure;subplot(2,1,1);stairs(x,y1./sum(y1),'r-');hold on;stairs(x,y2/sum(y2),'b-');plot([1.44,1.44],[0,0.1],'k-');
stairs(x,y11./sum(y11),'r:');stairs(x,y22/sum(y2),'b:');
xlabel('D_{Mmax}/R_{vir,host}');ylabel('Fraction');legend('DM sub','gas sub','D/R_{vir}=1.44','DM theo','gas theo');title('Distribution of Stripping Distance')
subplot(2,1,2);stairs(x,y3/sum(y3),'r-');hold on;stairs(x,y4/sum(y4),'b-.');
xlabel('D_{Mmax}/D_{Mmax,theory}');ylabel('Fraction');legend('dm','gas');
x=logspace(-5,0,10);[n,bin]=histc(mrate,x);
y=zeros(9,1);
for i=1:9
y(i)=median(rdmmax(bin==i&rdmmax>0));
end
figure;loglog(mrate(rdmmax>0),rdmmax(rdmmax>0),'.');
hold on;plot(mrate(rdmmax>0),d_strp(rdmmax>0),'r.');plot(mrate(rdmmax>0),d_comp(rdmmax>0),'g.');plot(mrate(rdmmax>0),rdmmax_theory(rdmmax>0),'ko');
plot(x(2:10),y,'-');plot(x,Dmean(x),'r-');

%% 
G=43007.1;
Mmax=zeros(Nhist,1);
Mmax0=zeros(Nhist,1);
Mmaxv=zeros(Nhist,1);
Mmax6=zeros(Nhist,1);
Mmax1=zeros(Nhist,1);
Mmaxc=zeros(Nhist,1);
Nodeinfall=zeros(Nhist,1)-1;
Diststrp=zeros(Nhist,1);
for h=1:Nhist
    Nnode=hist_bin{h}.Nnode;
    mdm=zeros(Nnode,1);
    mgas=zeros(Nnode,1);
    mhost=zeros(Nnode,1);
    kt=zeros(Nnode,1);
    Dist=zeros(Nnode,1);
    for s=1:Nnode
        node=hist_bin{h}.node(s);
        mdm(s)=node.mass(1);
        mgas(s)=node.mass(2);
        mhost(s)=node.mass(3);
        Vel=node.vel;       Pos=node.pos*a(node.Nsnap+1);%physical
        r=norm(Pos);        v=norm(Vel);
        kt(s)=v^2*(1-(Pos'*Vel/r/v)^2)/(G*(mhost(s)+mdm(s))*pmass/r);
        Dist(s)=r/halo{node.Nsnap+1}.Rvir(node.HostID+1,1)/a(node.Nsnap+1);
        if hist_bin{h}.node(s).Nsnap==Snappar(h,3)
            if mdm(s)>500&&mdm(s)/mhost(s)<0.05
            Nodeinfall(h)=s;
            end
        end
    end
    if Nodeinfall(h)>0
        dd=Dist./(2+kt).^(1/3)/1.25;
        dd0=Dist./2.^(1/3);
        dd1=Dist./3.^(1/3);
        dd6=Dist./1.8;
        ddv=Dist;
        ddc=Dist/1.5/3^(1/3);
        iDmax0=0;iDmax=0;iDmax6=0;iDmax1=0;iDmaxc=0;iDmaxv=Nodeinfall(h);
        flag=1;flag0=1;flag6=1;flagc=1;flag1=1;
        for s=1:Nodeinfall(h)
            if flag
            if dd(s)>1&&dd(s+1)<1
                if mdm(s)>mdm(s+1)
                iDmax=s;
                else
                iDmax=s+1;
                end
%                 flag=0;
                Diststrp(h)=Dist(iDmax);
            end
            end
            if flag0
            if dd0(s)>1&&dd0(s+1)<1
                if mdm(s)>mdm(s+1)
                iDmax0=s;
                else
                iDmax0=s+1;
                end
%                 flag0=0;
            end
            end
            if flag1
            if dd1(s)>1&&dd1(s+1)<1
                if mdm(s)>mdm(s+1)
                iDmax1=s;
                else
                iDmax1=s+1;
                end
%                 flag1=0;
            end
            end
            if flag6
            if dd6(s)>1&&dd6(s+1)<1
                if mdm(s)>mdm(s+1)
                iDmax6=s;
                else
                iDmax6=s+1;
                end
%                 flag6=0;
%                 Diststrp(h)=Dist(iDmax6);
            end
            end
            if flagc
            if ddc(s)>1&&ddc(s+1)<1
                if mdm(s)>mdm(s+1)
                iDmaxc=s;
                else
                iDmaxc=s+1;
                end
%                 flagc=0;
            end
            end
        end
        [tmp,iMmax]=max(mdm(1:Nodeinfall(h)));
        if iDmax0>0
            Mmax0(h)=mdm(iDmax0)/mdm(iMmax);
        else
            Mmax0(h)=-1;
        end
        if iDmax>0
            Mmax(h)=mdm(iDmax)/mdm(iMmax);
        else
            Mmax(h)=-1;
        end
        if iDmax1>0
            Mmax1(h)=mdm(iDmax1)/mdm(iMmax);
        else
            Mmax1(h)=-1;
        end
        if iDmax6>0
            Mmax6(h)=mdm(iDmax6)/mdm(iMmax);
        else
            Mmax6(h)=-1;
        end
        if iDmaxc>0
            Mmaxc(h)=mdm(iDmaxc)/mdm(iMmax);
        else
            Mmaxc(h)=-1;
        end
        Mmaxv(h)=mdm(iDmaxv)/mdm(iMmax);
    else
        Mmax0(h)=-1;Mmax(h)=-1;Mmax6(h)=-1;Mmax1(h)=-1;Mmaxc(h)=-1;Mmaxv(h)=-1;
    end
    
end

x=0:0.02:1.5;
y=histc(Mmax(Mmax>0),x);y=cumsum(y(end:-1:1));y=y(end:-1:1);
y0=histc(Mmax0(Mmax0>0),x);y0=cumsum(y0(end:-1:1));y0=y0(end:-1:1);
y6=histc(Mmax6(Mmax6>0),x);y6=cumsum(y6(end:-1:1));y6=y6(end:-1:1);
y1=histc(Mmax1(Mmax1>0),x);y1=cumsum(y1(end:-1:1));y1=y1(end:-1:1);
yc=histc(Mmaxc(Mmaxc>0),x);yc=cumsum(yc(end:-1:1));yc=yc(end:-1:1);
yv=histc(Mmaxv(Mmaxv>0),x);yv=cumsum(yv(end:-1:1));yv=yv(end:-1:1);
figure;stairs(x,y./y(1),'r-');
hold on;
stairs(x,y0./y0(1),'k-.');
stairs(x,y1./y1(1),'k:');
stairs(x,y6./y6(1),'b--');
stairs(x,yc./yc(1),'g');
stairs(x,yv./yv(1),'c');
legend('D=5/4*(2+\kappa_t)^{1/3}Rvir','D=2^{1/3}Rvir','D=3^{1/3}Rvir','D=1.8*Rvir','D=1.5*3^{1/3}Rvir','D=Rvir');
xlabel('Mtidal/Mmax');ylabel('Fraction(>Mtidal/Mmax)');title('m/M~[0.01,0.1]');
%%
for i=1:Nhist
if Nodeinfall(i)>0
dmi(i)=hist_bin{i}.node(Nodeinfall(i)).mass(1);
else
dmi(i)=0;
end
end
find(dmi>10000)
for i=1:Nhist
if hist_bin{i}.node(end).subid==34&&hist_bin{i}.node(end).Nsnap==99
i
end
end
%%
G=43007.1;pmass=0.008848;gmass=0.00156133;
% es=[];
% figure;
% for h=1:Nhist
h=36 %[26 88 176 244 371 396  627 724 767  940        1045        1751        1804        1981        2277        2660]
    Nnode=hist_bin{h}.Nnode;
    t=zeros(Nnode,1);
    rhist=zeros(Nnode,1);
    Rhist=zeros(Nnode,1);
    mdm=zeros(Nnode,1);
    mgas=zeros(Nnode,1);
    mhost=zeros(Nnode,1);
    mcen=zeros(Nnode,1);
    Ehist=zeros(Nnode,1);
    Ereal=zeros(Nnode,1);
    kin=zeros(Nnode,1);
    pot=zeros(Nnode,1);
    jhist=zeros(Nnode,1);
    Jhist=zeros(Nnode,1);
    Jreal=zeros(Nnode,1);
    skin=zeros(Nnode,1);
    spot=zeros(Nnode,1);
    sj=zeros(Nnode,1);
    muhist=zeros(Nnode,1);
    Thist=zeros(Nnode,1);
    kthist=zeros(Nnode,1);
    for s=1:Nnode
        snapnum=hist_bin{h}.node(s).Nsnap;
%         t(s)=age(snapnum+1);
        t(s)=log(a(snapnum+1));
        mdm(s)=hist_bin{h}.node(s).mass(1);
        mgas(s)=hist_bin{h}.node(s).mass(2);
        mhost(s)=hist_bin{h}.node(s).mass(3);
        mcen(s)=hist_bin{h}.node(s).mass(4);
        Vel=hist_bin{h}.node(s).vel;
        Pos=hist_bin{h}.node(s).pos*a(snapnum+1);%physical
        r=sqrt(sum(Pos.^2,1));
        v=sqrt(sum(Vel.^2,1));
        rhist(s)=r;
        Rhist(s)=r/halo{snapnum+1}.Rvir(hist_bin{h}.node(s).HostID+1,1)/a(snapnum+1);
        Mvir=halo{snapnum+1}.Mvir(hist_bin{h}.node(s).HostID+1,1);
        Rvir=halo{snapnum+1}.Rvir(hist_bin{h}.node(s).HostID+1,1)*a(snapnum+1);
        mu=G*(mdm(s)+mcen(s))*pmass;
        subid=hist_bin{h}.node(s).subid;
         kin(s)=orbpar{snapnum+1}(3,subid+1);
         pot(s)=orbpar{snapnum+1}(4,subid+1);
%          if kin(s)~=0
%          mu=-r*pot(s);
%          end
         muhist(s)=mu;
         Ereal(s)=kin(s)+pot(s);
%         Jreal(s)=sqrt(orbpar{snapnum+1}(5,subid+1));
%         Ehist(s)=0.5*v^2-mu/r;
Ehist(s)=0.5*v^2+pot(s);
        Ehost=-0.5*mu/r;%to be improved
        Ehist(s)=Ehist(s)/Ehost;
        Thist(s)=2*pi*sqrt((r/Ehist(s))^3/muhist(s));
        jhist(s)=r*v*sqrt(1-(Pos'*Vel/r/v)^2);
        Jhist(s)=jhist(s);
        jhost=sqrt(mu*r/Ehist(s));
        jhist(s)=jhist(s)/jhost;
        kthist(s)=v^2*(1-(Pos'*Vel/r/v)^2)/(mu/r);
        df_dt(s)=(12-Mratepar(h,2)^(-1./3)*v/r/sqrt(G*Mvir/2/Rvir^3));
    end
%     if Nnode>20&&isempty(find(Ehist<0, 1))&&mdm(1)>10000
%         h,Nnode
% s=Nodeinfall(h);
for s=1:Nnode
    if hist_bin{h}.node(s).Nsnap==Snappar(h,2)
       s=s-1;
       break;
    end
end   %pick strip as infall
figure('name','a');
% plot(t,Ehist./rhist.*(mdm+mcen)*pmass*G*(-0.5),'- .'); %set(gca,'xlim',[-0.4,0]);
% plot(t,Ehist./rhist.*pmass*G*(-0.5),'- .'); %set(gca,'xlim',[-0.4,0]);
% hold on; plot(t(s:end),Ehist(s:end)./rhist(s:end).*pmass*G*(-0.5),'ro'); %set(gca,'xlim',[-0.4,0]);
% plot(t,Ehist+skin+spot,'- .');
% hold on;plot(t,Ehist,'r-- .',t,skin+spot,'g-- .',t(s:end),Ehist(s:end),'ro');
% plot(t,Ehist./rhist,'r-- .',t(s:end),Ehist(s:end)./rhist(s:end),'ro');
plot(t,rhist./Ehist,'r-- .',t(s:end),rhist(s:end)./Ehist(s:end),'ro');
figure('name','stripping rate');
plot(df_dt,'o');
hold on;plot(diff(mgas/mgas(s))./diff(t),'.');
figure('name','cos(theta)');
cost=(Jhist.^2./muhist./rhist-1)./sqrt(1-jhist.^2);
plot(t,cost,'r-- .',t(s:end),cost(s:end),'ro');
figure('name','T');plot(t,Thist,'r-- .',t(s:end),Thist(s:end),'ro');
figure('name','J2');  
% plot(t,jhist.^2,'-- .');set(gca,'ylim',[-1,1]);
% hold on; plot(t(s:end),jhist(s:end).^2,'ro');
% plot(t,jhist+sj,'- .');
% hold on;plot(t,jhist,'r-- .',t,sj,'g-- .',t(s:end),jhist(s:end),'ro');
% hold on; plot(t,Jreal,'o');
plot(t,jhist.^2,'r-- .',t(s:end),jhist(s:end).^2,'ro');
figure('name','J');plot(t,Jhist,'r-- .',t(s:end),Jhist(s:end),'ro');
figure('name','M');
semilogy(t,mdm+mcen,'- .',t,mhost,'-- o',t,mdm,'-. *');%set(gca,'xlim',[-0.4,0]);
hold on; semilogy(t(s:end),mdm(s:end),'o');
figure('name','r');
plot(t,rhist,'- .');%set(gca,'xlim',[-0.4,0]);
hold on;plot(t(s:end),rhist(s:end),'ro');
figure('name','R/Rvir');
semilogy(t,Rhist,'- .');
hold on;semilogy(t(s:end),Rhist(s:end),'ro');
semilogy(t,(2+kthist).^(1/3),'--s');
grid on;
% figure('name','E');
% plot(t,Ereal,'- .');%set(gca,'xlim',[-0.4,0]);
% figure('name','a_real');
% plot(t,Ereal./(mdm+mcen),'- .');%set(gca,'xlim',[-0.4,0]);
% figure;
% plot(t,Ereal.*mdm,'- .');%set(gca,'xlim',[-0.4,0]);
figure;
 subplot(2,1,1); plot(t,mgas/mgas(s),'.');hold on;plot(t(s:end),mgas(s:end)/mgas(s),'ro');
 subplot(2,1,2); plot(t,mgas./mdm/(mgas(s)/mdm(s)),'.');hold on;plot(t(s:end),mgas(s:end)./mdm(s:end)/(mgas(s)/mdm(s)),'ro');
%     end
% end
% figure;
% s=Nodeinfall(h);
%  subplot(2,1,1); plot(cost,mgas/mgas(s),'- .');hold on;plot(cost(s:end),mgas(s:end)/mgas(s),'ro');
%  subplot(2,1,2); plot(cost,mgas./mdm/(mgas(s)/mdm(s)),'- .');hold on;plot(cost(s:end),mgas(s:end)./mdm(s:end)/(mgas(s)/mdm(s)),'ro');
%%
figure;
s=Nodeinfall(h);
 subplot(2,1,1); plot(theta,mgas(s:end)/mgas(s),'ro');
 subplot(2,1,2); plot(theta,mgas(s:end)./mdm(s:end)/(mgas(s)/mdm(s)),'ro');
%%
for h=1:Nhist
    Nnode=hist_bin{h}.Nnode;
    t=zeros(Nnode,1);
    rhist=zeros(Nnode,1);
    mdm=zeros(Nnode,1);
    mgas=zeros(Nnode,1);
    mhost=zeros(Nnode,1);
    mcen=zeros(Nnode,1);
    Ehist=zeros(Nnode,1);
    kin=zeros(Nnode,1);
    pot=zeros(Nnode,1);
    jhist=zeros(Nnode,1);
    Jhist=zeros(Nnode,1);
    muhist=zeros(Nnode,1);
    for s=1:Nnode
        t(s)=log(a(hist_bin{h}.node(s).Nsnap+1));
        mdm(s)=hist_bin{h}.node(s).mass(1);
        mgas(s)=hist_bin{h}.node(s).mass(2);
        mhost(s)=hist_bin{h}.node(s).mass(3);
        mcen(s)=hist_bin{h}.node(s).mass(4);
        snapnum=hist_bin{h}.node(s).Nsnap;
        Vel=hist_bin{h}.node(s).vel;
        Pos=hist_bin{h}.node(s).pos*a(snapnum+1);%physical
        r=sqrt(sum(Pos.^2,1));
        v=sqrt(sum(Vel.^2,1));
        rhist(s)=r;
        mu=G*(mdm(s)+mcen(s))*pmass;
        subid=hist_bin{h}.node(s).subid;
         kin(s)=orbpar{snapnum+1}(3,subid+1);
         pot(s)=orbpar{snapnum+1}(4,subid+1);
         if kin(s)~=0
         mu=-r*pot(s);
         end
         muhist(s)=mu;
        Ehist(s)=0.5*v^2-mu/r;
        Ehost=-0.5*mu/r;%to be improved
        Ehist(s)=Ehist(s)/Ehost;
        jhist(s)=r*v*sqrt(1-(Pos'*Vel/r/v)^2);
        Jhist(s)=jhist(s);
        jhost=sqrt(mu*r/Ehist(s));
        jhist(s)=jhist(s)/jhost;
    end

cost=(Jhist.^2./muhist./rhist-1)./sqrt(1-jhist.^2);
 s=Nodeinfall(h);
 smerge=Nnode+1;%host major merger happens at this snap
    s0=max(s,2);
    for sm=s0:Nnode
        if mhost(sm)>15*mhost(sm-1)
            smerge=sm;
            break;
        end
        if mgas(sm)<100
            smerge=sm;
            break;
        end
    end
%     smerge=Nnode+1;
    snapsmooth=s:smerge-1;
    if s>0&&~isempty(snapsmooth)
        if mhost(s)>100%&&isreal(j(h))&&j(h)>0.0&&j(h)<0.2
            subplot(2,1,1); plot(cost(snapsmooth),mgas(snapsmooth)/mgas(s),'.');hold on;
            subplot(2,1,2); plot(cost(snapsmooth),mgas(snapsmooth)./mdm(snapsmooth)/(mgas(s)/mdm(s)),'.');hold on;
        end
    end 
end 
hold off;

%%
    eJ=[];
    eE=[];
for h=1:Nhist
    Nnode=hist_bin{h}.Nnode;
    t=zeros(Nnode,1);
    rhist=zeros(Nnode,1);
    mdm=zeros(Nnode,1);
    mgas=zeros(Nnode,1);
    mhost=zeros(Nnode,1);
    mcen=zeros(Nnode,1);
    Ehist=zeros(Nnode,1);
    kin=zeros(Nnode,1);
    pot=zeros(Nnode,1);
    jhist=zeros(Nnode,1);
    Jhist=zeros(Nnode,1);
    muhist=zeros(Nnode,1);
    for s=1:Nnode
        t(s)=log(a(hist_bin{h}.node(s).Nsnap+1));
        mdm(s)=hist_bin{h}.node(s).mass(1);
        mgas(s)=hist_bin{h}.node(s).mass(2);
        mhost(s)=hist_bin{h}.node(s).mass(3);
        mcen(s)=hist_bin{h}.node(s).mass(4);
        snapnum=hist_bin{h}.node(s).Nsnap;
        Vel=hist_bin{h}.node(s).vel;
        Pos=hist_bin{h}.node(s).pos*a(snapnum+1);%physical
        r=sqrt(sum(Pos.^2,1));
        v=sqrt(sum(Vel.^2,1));
        rhist(s)=r;
        mu=G*(mdm(s)+mcen(s))*pmass;
        subid=hist_bin{h}.node(s).subid;
         kin(s)=orbpar{snapnum+1}(3,subid+1);
         pot(s)=orbpar{snapnum+1}(4,subid+1);
%          if kin(s)~=0
%          mu=-r*pot(s);
%          end
%          muhist(s)=mu;
        Ehist(s)=0.5*v^2+pot(s);
%         Ehost=-0.5*mu/r;%to be improved
%         Ehist(s)=Ehist(s)/Ehost;
        jhist(s)=r*v*sqrt(1-(Pos'*Vel/r/v)^2);
        Jhist(s)=jhist(s);
%         jhost=sqrt(mu*r/Ehist(s));
%         jhist(s)=jhist(s)/jhost;
    end
    if(Nodeinfall(h)>1&&Nodeinfall(h)<Nnode-1)
    eJ=[eJ,std(Jhist(Nodeinfall(h)-1:Nodeinfall(h)+2))/mean(Jhist(Nodeinfall(h)-1:Nodeinfall(h)+2))];
    eE=[eE,std(Ehist(Nodeinfall(h)-1:Nodeinfall(h)+2))/mean(Ehist(Nodeinfall(h)-1:Nodeinfall(h)+2))];
    end
end
figure;hist(eJ);
figure;hist(eE(find(eE>-0.2&eE<2)),-0.2:0.2:2);
%%
G=43007.1;pmass=0.008848;gmass=0.00156133;
% es=[];
% figure;
% for h=1:Nhist
h=27 %[26 88 176 244 371 396  627 724 767  940        1045        1751        1804        1981        2277        2660]
   
Nnode=hist_bin{h}.Nnode;
    for s=1:Nnode
        if hist_bin{h}.node(s).Nsnap==Snappar(h,2)
            Node_infall=s;
        end
    end
    
    t=zeros(Nnode,1);
    mdm=zeros(Nnode,1);
    mgas=zeros(Nnode,1);
    mhost=zeros(Nnode,1);
    mcen=zeros(Nnode,1);
    df_dt=zeros(Nnode,1);
    ft=ones(Nnode,1);
    v_r=zeros(Nnode,1); 
    for s=Node_infall:Nnode
        snapnum=hist_bin{h}.node(s).Nsnap;
%         t(s)=age(snapnum+1);
        t(s)=log(a(snapnum+1));
        mdm(s)=hist_bin{h}.node(s).mass(1);
        mgas(s)=hist_bin{h}.node(s).mass(2);
        mhost(s)=hist_bin{h}.node(s).mass(3);
        mcen(s)=hist_bin{h}.node(s).mass(4);
        Vel=hist_bin{h}.node(s).vel;
        Pos=hist_bin{h}.node(s).pos*a(snapnum+1);%physical
        r=sqrt(sum(Pos.^2,1));
        v=sqrt(sum(Vel.^2,1));
        scaleF=a(snapnum+1);
        Rvir=halo{snapnum+1}.Rvir(hist_bin{h}.node(s).HostID+1,1)*scaleF;
        Hz=HUBBLE0 * sqrt(Omega0 /scaleF^3+ (1 -Omega0 -OmegaLambda) / scaleF^2 +OmegaLambda);
        ft(s)=ft(s-1)+df_dt(s-1)*0.0116;
%         mrate=Mratepar(h,2);
        mrate=mdm(s)/mhost(s);
        df_dt(s)=(sqrt(178)-1.3*(ft(s))^0.5*mrate^(-1./3)*(v/r)/Hz)/2/pi/2.3;
%         df_dt(s)=-sqrt(178)+178/ft(s)*mrate^(1./3)*(r/v)*Hz;
%           df_dt(s)=sqrt(178)*(1-ft(s)*2*Rvir/r);  %tidal stripping
 %         v_r(s)=v/r;
        v_r(s)=mrate^(-1./3)*(v/r)/Hz;
%         if hist_bin{h}.node(s).Nsnap==Snappar(h,3)
%             Node_infall=s;
%         end
    end
% df_dt=12-v_r;
% df_dt=df_dt/6;
s=Node_infall;
t=t(s:end);
v_r=v_r(s:end);
mgas=mgas(s:end);
mdm=mdm(s:end);
df_dt=df_dt(s:end);
ft=ft(s:end);
figure;subplot(3,1,1);
% plot(v_r,'-');hold on;
plot(v_r,'o');
subplot(3,1,2);
% plot(df_dt,'-');hold on;
plot(diff(mgas/mgas(1))./diff(t),'.');hold on;
plot(df_dt,'o');
subplot(3,1,3);
% plot(cumsum(df_dt)*0.0116+1,'*');hold on;
plot(ft,'o');
hold on;plot(mgas/mgas(1),'- .');set(gca,'ylim',[0,1]);
plot(mgas./mdm/(mgas(1)/mdm(1)),'-*');
plot(mdm/mdm(1),'- p');
% ft=diff(mgas/mgas(s))./diff(t);
% filter=find(ft<-1);
% x=[x;v_r(filter)];y=[y;ft(filter)];
% cftool(v_r(filter),ft(filter))
%%
G=43007.1;pmass=0.008848;gmass=0.00156133;
% es=[];
% figure;
% for h=1:Nhist
h=27
   
Nnode=hist_bin{h}.Nnode;
    for s=1:Nnode
        if hist_bin{h}.node(s).Nsnap==Snappar(h,2)
            Node_infall=s;
        end
    end
    
    t=zeros(Nnode,1);
    mdm=zeros(Nnode,1);
    mgas=zeros(Nnode,1);
    mhost=zeros(Nnode,1);
    mcen=zeros(Nnode,1);
    df_dt=zeros(Nnode,1);
    ft=ones(Nnode,1);
    v_r=zeros(Nnode,1); 
    for s=Node_infall:Nnode
        snapnum=hist_bin{h}.node(s).Nsnap;
%         t(s)=age(snapnum+1);
        t(s)=log(a(snapnum+1));
        mdm(s)=hist_bin{h}.node(s).mass(1);
        mgas(s)=hist_bin{h}.node(s).mass(2);
        mhost(s)=hist_bin{h}.node(s).mass(3);
        mcen(s)=hist_bin{h}.node(s).mass(4);
        Vel=hist_bin{h}.node(s).vel;
        Pos=hist_bin{h}.node(s).pos*a(snapnum+1);%physical
        r=sqrt(sum(Pos.^2,1));
        v=sqrt(sum(Vel.^2,1));
        scaleF=a(snapnum+1);
        Rvir=halo{snapnum+1}.Rvir(hist_bin{h}.node(s).HostID+1,1)*scaleF;
        Hz=HUBBLE0 * sqrt(Omega0 /scaleF^3+ (1 -Omega0 -OmegaLambda) / scaleF^2 +OmegaLambda);
        ft(s)=ft(s-1)+df_dt(s-1)*0.0116;
%         mrate=Mratepar(h,2);
        mrate=mdm(s)/mhost(s);
%         df_dt(s)=(sqrt(178)-ft(s)*mrate^(-1./3)*(v/r)/Hz)/2/pi;
%         df_dt(s)=-sqrt(178)+178/ft(s)*mrate^(1./3)*(r/v)*Hz;
% c=5;
% rc=r/Rvir*c;
% gamma=(3-rc^2/(1+rc)^2/(log(1+rc)-rc/(1+rc)))^(1/3)*2/1.2;
gamma=1.3;
           df_dt(s)=sqrt(178)*(1-gamma*(ft(s)*Rvir/r)^(2/3))/2/pi/0.7;  %tidal stripping
 %         v_r(s)=v/r;
        v_r(s)=mrate^(-1./3)*(v/r)/Hz;
%         if hist_bin{h}.node(s).Nsnap==Snappar(h,3)
%             Node_infall=s;
%         end
    end
% df_dt=12-v_r;
% df_dt=df_dt/6;
s=Node_infall;
t=t(s:end);
v_r=v_r(s:end);
mdm=mdm(s:end);
df_dt=df_dt(s:end);
ft=ft(s:end);
figure;subplot(3,1,1);
% plot(v_r,'-');hold on;
plot(v_r,'o');
subplot(3,1,2);
% plot(df_dt,'-');hold on;
plot(diff(mdm/mdm(1))./diff(t),'.');hold on;
plot(df_dt,'o');
subplot(3,1,3);
% plot(cumsum(df_dt)*0.0116+1,'*');hold on;
plot(ft,'o');
hold on;plot(mdm/mdm(1),'- .');set(gca,'ylim',[0,1]);
% ft=diff(mgas/mgas(s))./diff(t);
% filter=find(ft<-1);
% x=[x;v_r(filter)];y=[y;ft(filter)];
% cftool(v_r(filter),ft(filter))
%%
figure('name','M');
semilogy(t,mdm+mcen,'- .',t,mhost,'-- o',t,mdm,'-. *');%set(gca,'xlim',[-0.4,0]);
hold on; semilogy(t(s:end),mdm(s:end),'o');
figure('name','r');
semilogy(t,rhist,'- .');%set(gca,'xlim',[-0.4,0]);
hold on;semilogy(t(s:end),rhist(s:end),'ro');grid on;
figure;
 subplot(2,1,1); plot(t,mgas/mgas(s),'.');hold on;plot(t(s:end),mgas(s:end)/mgas(s),'ro');
 subplot(2,1,2); plot(t,mgas./mdm/(mgas(s)/mdm(s)),'.');hold on;plot(t(s:end),mgas(s:end)./mdm(s:end)/(mgas(s)/mdm(s)),'ro');
