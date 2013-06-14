function [countM,countkt,kt_max]=mass_start_stat(hist_bin,Snappar,pmass,xbin,kbin)
global a
G=43007.1;

Nhist=size(hist_bin,1);
Mmax=zeros(Nhist,1);
Mmax0=zeros(Nhist,1);
Mmaxv=zeros(Nhist,1);
Mmax6=zeros(Nhist,1);
Mmax1=zeros(Nhist,1);
Mmaxc=zeros(Nhist,1);
kt_enter=zeros(Nhist,1);
kt_max=zeros(Nhist,1);
Nodeinfall=zeros(Nhist,1)-1;
% Diststrp=zeros(Nhist,1);
for h=1:Nhist
    Nnode=hist_bin{h}.Nnode;
    mdm=zeros(Nnode,1);
%     mgas=zeros(Nnode,1);
    mhost=zeros(Nnode,1);
    kt=zeros(Nnode,1);
    Dist=zeros(Nnode,1);
    for s=1:Nnode
        node=hist_bin{h}.node(s);
        mdm(s)=node.mass(1);
%         mgas(s)=node.mass(5);
        mhost(s)=node.mass(2);
        Vel=node.vel;       Pos=node.pos*a(node.Nsnap+1);%physical
        r=norm(Pos);        v=norm(Vel);
        kt(s)=v^2*(1-(Pos'*Vel/r/v)^2)/(G*(mhost(s)+mdm(s))*pmass/r);
        Dist(s)=r/comoving_virial_radius(mhost(s)*pmass,a(node.Nsnap+1))/a(node.Nsnap+1);
        if hist_bin{h}.node(s).Nsnap==Snappar(h,3)
            if mdm(s)>300
            Nodeinfall(h)=s;
            end
        end
    end
    if Nodeinfall(h)>0
        dd=Dist./(2+kt).^(1/3)/1.25;
%         dd=Dist./(2+kt).^(1/3);
        dd0=Dist./2.^(1/3);
        dd1=Dist./3.^(1/3);
        dd6=Dist./1.8;
        ddv=Dist;
        ddc=Dist/1.5/3^(1/3);
        iDmax0=0;iDmax=0;iDmax6=0;iDmax1=0;iDmaxc=0;iDmaxv=Nodeinfall(h);
        flag=1;flag0=1;flag6=1;flagc=1;flag1=1;
        for s=1:Nodeinfall(h)-1
            if flag
            if dd(s)>1&&dd(s+1)<1
                if mdm(s)>mdm(s+1)
                iDmax=s;
                else
                iDmax=s+1;
                end
                kt_enter(h)=kt(iDmax);
%                 flag=0;
%                 Diststrp(h)=Dist(iDmax);
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
        kt_max(h)=kt(iMmax);
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

x=xbin;
% x=0:0.02:1.2;
y=histc(Mmax(Mmax>0),x);
y0=histc(Mmax0(Mmax0>0),x);
y6=histc(Mmax6(Mmax6>0),x);
y1=histc(Mmax1(Mmax1>0),x);
yc=histc(Mmaxc(Mmaxc>0),x);
yv=histc(Mmaxv(Mmaxv>0),x);
ykt=histc(kt_enter(Mmax>0),kbin);
countM=[y,y0,y6,y1,yc,yv];
countkt=ykt;
y=cumsum(y(end:-1:1));y=y(end:-1:1);
y0=cumsum(y0(end:-1:1));y0=y0(end:-1:1);
y6=cumsum(y6(end:-1:1));y6=y6(end:-1:1);
y1=cumsum(y1(end:-1:1));y1=y1(end:-1:1);
yc=cumsum(yc(end:-1:1));yc=yc(end:-1:1);
yv=cumsum(yv(end:-1:1));yv=yv(end:-1:1);
figure;stairs(x,y./y(1),'r-');
hold on;
stairs(x,y0./y0(1),'k-.');
stairs(x,y1./y1(1),'k:');
stairs(x,y6./y6(1),'b--');
stairs(x,yc./yc(1),'g');
stairs(x,yv./yv(1),'c');
legend('D=5/4*(2+\kappa_t)^{1/3}Rvir','D=2^{1/3}Rvir','D=3^{1/3}Rvir','D=1.8*Rvir','D=1.5*3^{1/3}Rvir','D=Rvir');
% legend('D=(2+\kappa_t)^{1/3}Rvir','D=2^{1/3}Rvir','D=3^{1/3}Rvir','D=1.8*Rvir','D=1.5*3^{1/3}Rvir','D=Rvir');
xlabel('Mtidal/Mmax');ylabel('Fraction(>Mtidal/Mmax)');title('m/M~[0.01,0.1]');
ykt=ykt(1:end-1)'./diff(kbin);
figure;loglog(kbin(1:end-1),ykt,'. -');