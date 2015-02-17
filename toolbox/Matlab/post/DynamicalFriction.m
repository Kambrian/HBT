close all;
datadir='/mnt/A4700/data/6702/subcat/profile/logbin';
G=43007.1;pmass=0.008848;gmass=0.00156133;
HUBBLE0=0.1;Omega0=0.3;OmegaLambda=0.7;
h=23 %[26 88 176 244 371 396  627 724 767  940        1045        1751        1804        1981        2277        2660]
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
    muhist=zeros(Nnode,1);
    DFJ=zeros(Nnode,1);
    for s=Nodeinfall(h):Nnode
        snapnum=hist_bin{h}.node(s).Nsnap;
        scaleF=a(snapnum+1);
        t(s)=log(scaleF);
        mdm(s)=hist_bin{h}.node(s).mass(1);
        mgas(s)=hist_bin{h}.node(s).mass(2);
        mhost(s)=hist_bin{h}.node(s).mass(3);
        mcen(s)=hist_bin{h}.node(s).mass(4);
        Vel=hist_bin{h}.node(s).vel;
        Pos=hist_bin{h}.node(s).pos*scaleF;%physical
        r=sqrt(sum(Pos.^2,1));
        v=sqrt(sum(Vel.^2,1));
        Rhist(s)=r/halo{snapnum+1}.Rvir(hist_bin{h}.node(s).HostID+1,1)/a(snapnum+1);
        grpind=hist_bin{h}.node(s).HostID+1;
        haloprof=readhalo_prof_single(datadir,snapnum,'halo',halo{snapnum+1},grpind);
        dynamprof_sph=readdynam_prof_sph_single(datadir,snapnum,halo{snapnum+1},grpind);
        vbin=diff(logspace(log10(max(1.5,halo{snapnum+1}.rmax(grpind)*1e-2)),log10(halo{snapnum+1}.rmax(grpind)),halo{snapnum+1}.nbin(grpind)+1).^3)';
        vbin=vbin*(4*pi/3)*scaleF^3;
        [tmp,ibin]=min(abs(haloprof.rs*scaleF-r));
        row=(haloprof.ns+haloprof.no+haloprof.nb)./vbin*pmass;
        rowbin=row(ibin);
        [tmp,ibin2]=min(abs(log(row/rowbin/10)));
        [tmp,ibin3]=min(abs(log(row*10/rowbin)));
        bmax=sqrt((r-haloprof.rs(ibin2)*scaleF)*(haloprof.rs(ibin3)*scaleF-r));
        bmax=max(abs(r-haloprof.rs(ibin2)*scaleF),abs(haloprof.rs(ibin3)*scaleF-r));
        if imag(bmax)~=0
            bmax=r;
        end
        sigma2=mean(dynamprof_sph.vd_rtf(:,ibin))*scaleF;
        Hz=HUBBLE0 * sqrt(Omega0 /scaleF^3+ (1 -Omega0 -OmegaLambda) / scaleF^2 +OmegaLambda);
        Hratio=Hz/HUBBLE0;
        OmegaZ=Omega0/(scaleF*scaleF*scaleF)/Hratio/Hratio;
        x=1-OmegaZ;
        virialF=18.0*3.1416*3.1416+82.0*x-39.0*x*x;%<Rho_vir>/Rho_cri,tophat
        Rsat=(2*G*mdm(s)*pmass/virialF/Hz^2)^(1./3);
        Rhost=halo{snapnum+1}.Rvir(grpind,1)*scaleF;
        Lambda=r/1.4/Rsat
        Lambda=mhost(s)/mdm(s)*r/Rhost
        Lambda=mhost(s)/mdm(s)*bmax/Rhost
        x=v/sqrt(2*sigma2)
        DFa=-4*pi*log(Lambda)*G^2*rowbin*(mdm(s)*pmass+mgas(s)*gmass)/v^3*(erf(x)-2*x/sqrt(pi)*exp(-x^2)); %dv/dt=DFa*v
        if isempty(DFa)
            DFJ(s)=NaN;
        else
        DFJ(s)=DFa/Hz;
        end
        rhist(s)=r;
        mu=G*(mdm(s)+mcen(s))*pmass;
        subid=hist_bin{h}.node(s).subid;
         kin(s)=orbpar{snapnum+1}(3,subid+1);
         pot(s)=orbpar{snapnum+1}(4,subid+1);
         if kin(s)~=0
         mu=-r*pot(s);
         end
         muhist(s)=mu;
         Ereal(s)=kin(s)+pot(s);
%         Jreal(s)=sqrt(orbpar{snapnum+1}(5,subid+1));
%         Ehist(s)=0.5*v^2-mu/r;
        Ehist(s)=0.5*v^2+pot(s);
        Ehost=-0.5*mu/r;%to be improved
        Ehist(s)=Ehist(s)/Ehost;
        jhist(s)=r*v*sqrt(1-(Pos'*Vel/r/v)^2);
        Jhist(s)=jhist(s);
        jhost=sqrt(mu*r/Ehist(s));
        jhist(s)=jhist(s)/jhost;
    end

s=Nodeinfall(h);
figure('name','J2');  
plot(t,jhist.^2,'r-- .',t(s:end),jhist(s:end).^2,'ro');
figure('name','M');
semilogy(t,mdm+mcen,'- .',t,mhost,'-- o',t,mdm,'-. *');%set(gca,'xlim',[-0.4,0]);
hold on; semilogy(t(s:end),mdm(s:end),'o');
% figure;
% s=Nodeinfall(h);
%  subplot(2,1,1); plot(t,mgas/mgas(s),'.');hold on;plot(t(s:end),mgas(s:end)/mgas(s),'ro');
%  subplot(2,1,2); plot(t,mgas./mdm/(mgas(s)/mdm(s)),'.');hold on;plot(t(s:end),mgas(s:end)./mdm(s:end)/(mgas(s)/mdm(s)),'ro');

figure('name','r');
plot(t,rhist,'- .');%set(gca,'xlim',[-0.4,0]);
hold on;plot(t(s:end),rhist(s:end),'ro');
figure('name','R/Rvir');
plot(t,Rhist,'- .');
hold on;plot(t(s:end),Rhist(s:end),'ro');
figure('name','J');plot(t,Jhist,'r-- .',t(s:end),Jhist(s:end),'ro');
figure;
subplot(2,1,1);plot(s:Nnode,Jhist(s:end),'- o');
subplot(2,1,2);plot(s+1:Nnode,diff(log(Jhist(s:Nnode)))./diff(t(s:Nnode)),'- .',s:Nnode,DFJ(s:end),'-- o');

