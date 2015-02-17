function xi2=xi2_GASstrp_history(pmass,csat,r,v,t,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,par_gas,flag_plot,flag_redistr)

Ninfall=numel(r);
df_dt_gas=zeros(Ninfall,1);
ft_gas=ones(Ninfall,1);
dr_dt=ones(Ninfall,1);rt=ones(Ninfall,1);
mass=@(x) log(1+x)-x./(1+x);
% vc2=@(x) (log(1+x)-x./(1+x))./x;
    function y=vc2(x)
        if x<=0
            y=0;
        else
            y=(log(1+x)-x./(1+x))./x;
        end
    end
           
rhox=@(x) x./(1+x).^2;
rhoc=@(c) c./mass(c);

nfw_mass=@(x) log(1+x)-x./(1+x);
nfw_rho=@(x) x.*(1+x).^2;
% vs2=@(x) nfw_mass(x).^1.3./nfw_rho(x).^(2/3);
alp=1.9;
vs2=@(x) x.^(2/3*(alp-1)).*(1+x).^(-4/3);
G=43007.1;
rv=comoving_virial_radius(mdm(1)*pmass,exp(t(1)))*exp(t(1));
% coef=G*mdm(1)*pmass/rv/vs2(csat);
coef=5/9*Umsat(1)/(quadl(vs2,0,csat)/csat);
vs_inv=@(x) 1./sqrt(coef*vs2(x));

for i=1:Ninfall
    if i>1
        ft_gas(i)=ft_gas(i-1)+df_dt_gas(i-1)*(t(i)-t(i-1));
    end
    switch flag_redistr
        case -3
             if ft_gas(i)<0
%                 for ii=i:Ninfall
%                     ft_gas(ii)=0;
%                 end
%                 break;
                ft_gas(i)=0;
            end
            if i>1
%                 rt(i)=rt(i-1)+dr_dt(i-1)*(t(i)-t(i-1));
            rt(i)=radius_fraction(ft_gas(i),csat);
%             rt0=ft_gas(i)*vc2(csat)/vc2(ft_gas(i)*csat);
%             rt(i)=fzero(@(x) vc2(x*csat)/vc2(csat)*x-ft_gas(i),rt0);
%             if isnan(rt(i))
%                 rt(i)=rt0;
%             end
            end
            grate=mgas(i)/mdm(i)/fhost(i);
%             grate=1;
            rb0=par_gas(2)*sqrt(grate)*(mdm(i)/mhost(i))^(1/3)*r(i)/Rvir(i)*sqrt(G*mhost(i)*pmass/Rvir(i))/v(i);
            rb=gas_binding_radius(rb0,csat,rhox(r(i)/Rvir(i)*chost(i))*rhoc(chost(i)));
%             rb0=rb0*sqrt(vc2(rb0*csat)/vc2(csat))*sqrt(rhox(rb0*csat)*rhoc(csat)/rhox(r(i)/Rvir(i)*chost(i))/rhoc(chost(i)));
%             rb=fzero(@(x) (mdm(1)/mhost(i))^(1/3)*r(i)/Rvir(i)*sqrt(G*mhost(i)*pmass/Rvir(i))/v(i)*sqrt(vc2(x*csat)/vc2(csat))*sqrt(rhox(x*csat)*rhoc(csat)/rhox(r(i)/Rvir(i)*chost(i))/rhoc(chost(i)))-x,rb0);
%             if isnan(rb)
%                 rb=rb0;
%             end
%             timescale=2*pi*rb/sqrt(virialF(i)/2); % isothermal azimuthal sound crossing time
              timescale=rt(i)*rv*vs_inv(rt(i)*csat)*Hz(i)*par_gas(1);
%               timescale=2*pi*rb*rv*vs_inv(rb*csat)*Hz(i);
            df_dt_gas(i)=(vc2(rb*csat)/vc2(csat)*rb-vc2(rt(i)*csat)/vc2(csat)*rt(i))/timescale;
         case -2 %par_gas=[1.5,1]
             if ft_gas(i)<0
%                 for ii=i:Ninfall
%                     ft_gas(ii)=0;
%                 end
%                 break;
                ft_gas(i)=0;
            end
            if i>1
%                 rt(i)=rt(i-1)+dr_dt(i-1)*(t(i)-t(i-1));
            rt(i)=radius_fraction(ft_gas(i),csat);
%             rt0=ft_gas(i)*vc2(csat)/vc2(ft_gas(i)*csat);
%             rt(i)=fzero(@(x) vc2(x*csat)/vc2(csat)*x-ft_gas(i),rt0);
%             if isnan(rt(i))
%                 rt(i)=rt0;
%             end
            end
            grate=mgas(i)/mdm(i)/fhost(i);
%             grate=1;
            rb0=par_gas(2)*sqrt(grate)*(mdm(i)/mhost(i))^(1/3)*r(i)/Rvir(i)*sqrt(G*mhost(i)*pmass/Rvir(i))/v(i);
            rb=gas_binding_radius(rb0,csat,rhox(r(i)/Rvir(i)*chost(i))*rhoc(chost(i)));
%             rb0=rb0*sqrt(vc2(rb0*csat)/vc2(csat))*sqrt(rhox(rb0*csat)*rhoc(csat)/rhox(r(i)/Rvir(i)*chost(i))/rhoc(chost(i)));
%             rb=fzero(@(x) (mdm(1)/mhost(i))^(1/3)*r(i)/Rvir(i)*sqrt(G*mhost(i)*pmass/Rvir(i))/v(i)*sqrt(vc2(x*csat)/vc2(csat))*sqrt(rhox(x*csat)*rhoc(csat)/rhox(r(i)/Rvir(i)*chost(i))/rhoc(chost(i)))-x,rb0);
%             if isnan(rb)
%                 rb=rb0;
%             end
            timescale=quad(vs_inv,0.001,rt(i)*csat)*rv/csat*Hz(i)*par_gas(1);             %integrated sound crossing time
            df_dt_gas(i)=(vc2(rb*csat)/vc2(csat)*rb-vc2(rt(i)*csat)/vc2(csat)*rt(i))/timescale;
        case -1 %par=[2,1]
            if ft_gas(i)<0
%                 for ii=i:Ninfall
%                     ft_gas(ii)=0;
%                 end
%                 break;
                ft_gas(i)=0;
            end
            if i>1
%                 rt(i)=rt(i-1)+dr_dt(i-1)*(t(i)-t(i-1));
            rt(i)=radius_fraction(ft_gas(i),csat);
%             rt0=ft_gas(i)*vc2(csat)/vc2(ft_gas(i)*csat);
%             rt(i)=fzero(@(x) vc2(x*csat)/vc2(csat)*x-ft_gas(i),rt0);
%             if isnan(rt(i))
%                 rt(i)=rt0;
%             end
            end
            grate=mgas(i)/mdm(i)/fhost(i);
%             grate=1;
            rb0=par_gas(2)*sqrt(grate)*(mdm(i)/mhost(i))^(1/3)*r(i)/Rvir(i)*sqrt(G*mhost(i)*pmass/Rvir(i))/v(i);
            rb=gas_binding_radius(rb0,csat,rhox(r(i)/Rvir(i)*chost(i))*rhoc(chost(i)));
%             rb0=rb0*sqrt(vc2(rb0*csat)/vc2(csat))*sqrt(rhox(rb0*csat)*rhoc(csat)/rhox(r(i)/Rvir(i)*chost(i))/rhoc(chost(i)));
%             rb=fzero(@(x) (mdm(1)/mhost(i))^(1/3)*r(i)/Rvir(i)*sqrt(G*mhost(i)*pmass/Rvir(i))/v(i)*sqrt(vc2(x*csat)/vc2(csat))*sqrt(rhox(x*csat)*rhoc(csat)/rhox(r(i)/Rvir(i)*chost(i))/rhoc(chost(i)))-x,rb0);
%             if isnan(rb)
%                 rb=rb0;
%             end
            timescale=2*pi*rb/sqrt(virialF(i)/2)/sqrt(vc2(rb*csat)/vc2(csat))*par_gas(1); % dynamical time
            df_dt_gas(i)=(vc2(rb*csat)/vc2(csat)*rb-vc2(rt(i)*csat)/vc2(csat)*rt(i))/timescale;
        case 0
            % Both gas and DM sat do not redistribute
            mrate=mdm(1)/mhost(i);
            grate=mgas(1)/mdm(1)/fhost(i);
            df_dt_gas(i)=(sqrt(virialF(i))*par_gas(1)-par_gas(2)*(ft_gas(i))/sqrt(grate)*mrate^(-1./3)*(v(i)/r(i)/Hz(i)))/2/pi;
        case 1
            % gas redistribute, but not DM
            mrate=mdm(1)/mhost(i);
            grate=mgas(i)/mdm(1)/fhost(i);
            df_dt_gas(i)=(sqrt(virialF(i))*par_gas(1)-par_gas(2)*(ft_gas(i))/sqrt(grate)*mrate^(-1./3)*(v(i)/r(i)/Hz(i)))/2/pi;
        case 2
            % Both redistribute
            mrate=mdm(i)/mhost(i);
            grate=mgas(i)/mdm(i)/fhost(i);
%             grate=1;
%             grate=mgas(1)/mdm(1)*ft_gas(i)/fhost(i);
            df_dt_gas(i)=(sqrt(virialF(i))*par_gas(1)-par_gas(2)*ft_gas(i)/sqrt(grate)*mrate^(-1./3)*(v(i)/r(i)/Hz(i)))/2/pi;
%             df_dt_gas(i)=(sqrt(virialF(i))*par_gas(1)-par_gas(2)*4*(ft_gas(i))/sqrt(grate)*mhost(i)^(0.2)*mdm(i)^(-0.3)*(v(i)/r(i)/Hz(i)))/2/pi;
%             df_dt_gas(i)=(sqrt(virialF(i))*par_gas(1)-par_gas(2)*5*(ft_gas(i))/sqrt(grate)*mhost(i)^(0.2)*mdm(i)^(-1./3)*(v(i)/r(i)/Hz(i)))/2/pi;
%             df_dt_gas(i)=(sqrt(virialF(i))*2.4*mdm(i)^(-0.1)*par_gas(1)-par_gas(2)*5*(ft_gas(i))/sqrt(grate)*mhost(i)^(0.2)*mdm(i)^(-1./3)*(v(i)/r(i)/Hz(i)))/2/pi;
%             df_dt_gas(i)=(sqrt(virialF(i))*par_gas(1)-par_gas(2)*3*(ft_gas(i))/sqrt(grate)*mhost(1)^(-0.1)*mrate^(-1./3)*(v(i)/r(i)/Hz(i)))/2/pi;
%                 df_dt_gas(i)=(sqrt(virialF(i))*par_gas(1)-par_gas(2)*(ft_gas(i)^par_gas(3))/sqrt(grate)*mrate^(-1./3)*(v(i)/r(i)/Hz(i)))/2/pi;
%                 df_dt_gas(i)=(sqrt(virialF(i))*par_gas(1)-par_gas(2)*(ft_gas(i))/(grate)^par_gas(3)*mrate^(-1./3)*(v(i)/r(i)/Hz(i)))/2/pi;
            % Complete redistribution
        case 3
            mrate=mdm(i)/mhost(i);
            grate=mgas(i)/mdm(i)/fhost(i);
            df_dt_gas(i)=(sqrt(virialF(i))*par_gas(1)-par_gas(2)/sqrt(grate)*mrate^(-1./3)*(v(i)/r(i)/Hz(i)))*ft_gas(i)/2/pi;
        case 4
            frate=mgas(i)/mdm(i)/(mgas(1)/mdm(1));
%             df_dt_gas(i)=sqrt(virialF(i))*(1*par_gas(1)-par_gas(2)*frate^(1/3)*(ft_gas(i)*Rvir(i)/r(i))^(2/3))/2/pi;  %tidal stripping
%             df_dt_gas(i)=sqrt(virialF(i))*(1*par_gas(1)-par_gas(2)*(ft_gas(i)*Rvir(i)/r(i))^(2/3))/2/pi;  %tidal stripping
%             df_dt_gas(i)=sqrt(virialF(i))*ft_gas(i)*(par_gas(1)-par_gas(2)*(Rvir(i)/r(i))^(2/3))/2/pi;  %tidal stripping, with redistribution
            df_dt_gas(i)=sqrt(virialF(i))*(1*par_gas(1)-par_gas(2)*ft_gas(i)*Rvir(i)/r(i))/2/pi;  %tidal stripping
%             df_dt_gas(i)=sqrt(virialF(i))*ft_gas(i)*(1*par_gas(1)-par_gas(2)*Rvir(i)/r(i))/2/pi;  %tidal stripping
%             df_dt_gas(i)=sqrt(virialF(i))*ft_gas(i)*(1-2*Rvir(i)/r(i))/2/pi/par_gas(2);  %tidal stripping
        case 5
            mrate=mdm(i)/mhost(i);
            grate=mgas(i)/mdm(i)/fhost(i);
            df_dt_gas(i)=(sqrt(virialF(i))*par_gas(1)-...
                sqrt(par_gas(2)/grate+3*virialF(i)*mrate^(2/3)*(Rvir(i)/r(i))^2*(r(i)*Hz(i)/v(i))^2)...
                *(ft_gas(i))*mrate^(-1./3)*(v(i)/r(i)/Hz(i)))/2/pi;
%             df_dt_gas(i)=(sqrt(virialF(i))*par_gas(1)-...
%                 sqrt(par_gas(2)/grate+3*virialF(i)*mrate^(2/3)*(Rvir(i)/r(i))^2*(r(i)*Hz(i)/v(i))^2)...
%                 *mrate^(-1./3)*(v(i)/r(i)/Hz(i)))*(ft_gas(i))/2/pi;       %redistribution
        otherwise
            error('wrong flag_redistr');
    end
end
frac_gas=mgas(1:Ninfall)/mgas(1);
% weight=r(1:Ninfall)./Rvir(1:Ninfall);
weight=ones(Ninfall,1);
% weight=frac_gas;
% weight=1./mgas;
filter=logical(abs(ft_gas)<inf&(~isnan(ft_gas))&(imag(ft_gas)==0));mscale_gas=sum(frac_gas(filter).*ft_gas(filter).*weight(filter))/sum(frac_gas(filter).^2.*weight(filter));
if mscale_gas>1.01
    mscale_gas=1.01;
else if mscale_gas<0.99
        mscale_gas=0.99;
    end
end
% mscale_gas=1;
frac_gas(1)=1/mscale_gas;
xi2=sum((frac_gas*mscale_gas-ft_gas).^2.*weight)/sum(weight);

da=t-t(1);
switch flag_plot
    case 1
    figure;
    plot(da,frac_gas*mscale_gas,'*');hold on;
    plot(da,ft_gas,'--');
    xlabel('\Delta ln(a)');ylabel('remained fraction');
    set(gca,'ylim',[0,1]);
    text(.05,.1,num2str([r(1:3)./Rvir(1:3);mscale_gas]));
    legend('gas fraction',['ram pressure model',num2str(sqrt(xi2)),sprintf('\n'),num2str(par_gas)],'location','best');
    hold off;
    case 2
    plot(frac_gas,ft_gas,'o','markersize',2,'markerfacecolor','k');hold on;
    case 3
    subplot(2,1,1);plot(da,frac_gas*mscale_gas,'o');hold on;
    subplot(2,1,2);plot(da,ft_gas,'.');hold on;
    otherwise
end
end