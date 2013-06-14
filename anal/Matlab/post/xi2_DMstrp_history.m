function xi2=xi2_DMstrp_history(Hz,csat,r,vt2,kt,t,mdm,chost,Rvir,virialF,par_dm,flag_plot,flag_redistr)

Ninfall=numel(r);
df_dt_dm=zeros(Ninfall,1);
ft_dm=ones(Ninfall,1);
dr_dt=ones(Ninfall,1);rt=ones(Ninfall,1);
mass=@(x) log(1+x)-x./(1+x);
p=@(x,k) (2+k).*mass(x)-(x./(1+x)).^2;
% vc2=@(x) (log(1+x)-x./(1+x))./x;
function y=vc2(x)
        if x<=0
            y=0;
        else
            y=(log(1+x)-x./(1+x))./x;
        end
    end
mgrad=@(x) x.^2./(1+x).^2./mass(x); %dlnM/dlnx, x=r/rs

for i=1:Ninfall
    if i>1
        ft_dm(i)=ft_dm(i-1)+df_dt_dm(i-1)*(t(i)-t(i-1));
    end
    switch flag_redistr
            case -6  %accurate sat+host density prof, orbital timescale, mass rather than radius stripping 
                % par=[0.6,0.4]
            if ft_dm(i)<0
                ft_dm(i)=0;
            end    
            if i>1
%                 rt(i)=rt(i-1)+dr_dt(i-1)*(t(i)-t(i-1));
            rt(i)=radius_fraction(ft_dm(i),csat);
            end
%             rb=(r(i)/Rvir(i))^(3/2)*sqrt(mass(chost(i)))/sqrt(p(r(i)/Rvir(i)*chost(i),kt(i)));
%             rb=rb*vc2(rb*csat)/vc2(csat)/par_dm(2);
            rb0=par_dm(2)*r(i)/Rvir(i)/sqrt(vc2(r(i)/Rvir(i)*chost(i))/vc2(chost(i)))/sqrt(3-mgrad(r(i)/Rvir(i)*chost(i)));
            rb=dm_binding_radius(rb0,csat);
            timescale=par_dm(1)*(2*pi*r(i)*Hz(i)/sqrt(vt2(i)));
            df_dt_dm(i)=(vc2(rb*csat)/vc2(csat)*rb-vc2(rt(i)*csat)/vc2(csat)*rt(i))/timescale;
        case -5  %accurate sat+host density prof, self timescale, mass rather than radius stripping 
            % par=[1,0.5]
            if ft_dm(i)<0
                ft_dm(i)=0;
            end
            if i>1
%                 rt(i)=rt(i-1)+dr_dt(i-1)*(t(i)-t(i-1));
            rt(i)=radius_fraction(ft_dm(i),csat);
            end
%             rb=(r(i)/Rvir(i))^(3/2)*sqrt(mass(chost(i)))/sqrt(p(r(i)/Rvir(i)*chost(i),kt(i)));
%             rb=rb*vc2(rb*csat)/vc2(csat)/par_dm(2);
            rb0=par_dm(2)*r(i)/Rvir(i)/sqrt(vc2(r(i)/Rvir(i)*chost(i))/vc2(chost(i)))/sqrt(3-mgrad(r(i)/Rvir(i)*chost(i)));
            rb=dm_binding_radius(rb0,csat);
            timescale=par_dm(1)*2*pi*rb/(sqrt(vc2(rb*csat)/vc2(csat))*sqrt(virialF(i)/2));
            df_dt_dm(i)=(vc2(rb*csat)/vc2(csat)*rb-vc2(rt(i)*csat)/vc2(csat)*rt(i))/timescale;
        case -4  %accurate sat+host density prof, orbital timescale
            if i>1
                rt(i)=rt(i-1)+dr_dt(i-1)*(t(i)-t(i-1));
            end
%             rb=(r(i)/Rvir(i))^(3/2)*sqrt(mass(chost(i)))/sqrt(p(r(i)/Rvir(i)*chost(i),kt(i)));
%             rb=rb*vc2(rb*csat)/vc2(csat)/par_dm(2);
            rb0=par_dm(2)*r(i)/Rvir(i)/sqrt(vc2(r(i)/Rvir(i)*chost(i))/vc2(chost(i)))/sqrt(3-mgrad(r(i)/Rvir(i)*chost(i)));
            rb=dm_binding_radius(rb0,csat);
%             rb=0;
            timescale=par_dm(1)*(2*pi*r(i)*Hz(i)/sqrt(vt2(i)));
            dr_dt(i)=(rb-rt(i))/timescale;
        case -3  %accurate sat+host density prof
            if i>1
                rt(i)=rt(i-1)+dr_dt(i-1)*(t(i)-t(i-1));
            end
%             rb=(r(i)/Rvir(i))^(3/2)*sqrt(mass(chost(i)))/sqrt(p(r(i)/Rvir(i)*chost(i),kt(i)));
% %             rb=(r(i)/Rvir(i))/1.4/sqrt(2);
%             rb=rb*vc2(rb*csat)/vc2(csat)/par_dm(2);
            rb0=par_dm(2)*r(i)/Rvir(i)/sqrt(vc2(r(i)/Rvir(i)*chost(i))/vc2(chost(i)))/sqrt(3-mgrad(r(i)/Rvir(i)*chost(i)));
            rb=dm_binding_radius(rb0,csat);
%             rd=min(rb,rt(i));
            timescale=par_dm(1)*2*pi*rb/(sqrt(vc2(rb*csat)/vc2(csat))*sqrt(virialF(i)/2));
            dr_dt(i)=(rb-rt(i))/timescale;
        case -2  %accurate host density prof
            df_dt_dm(i)=sqrt(virialF(i))*(1-par_dm(2)*ft_dm(i)*(Rvir(i)/r(i))^(3/2)*sqrt(p(r(i)/Rvir(i)*chost(i),kt(i)))/sqrt(mass(chost(i))))/2/pi/sqrt(2)/par_dm(1);  %tidal stripping
        case -1  %kt
            df_dt_dm(i)=sqrt(virialF(i))*(1-par_dm(2)*1.4*sqrt(1+kt(i)).*(ft_dm(i)*Rvir(i)/r(i)))/2/pi/sqrt(2)/par_dm(1);  %tidal stripping
        case 0
%             par_dm=[0.36,1.5];
%                 df_dt_dm(i)=sqrt(virialF(i))*(1*par_dm(1)-par_dm(2)*(ft_dm(i)*Rvir(i)/r(i))^(2/3))/2/pi;  %tidal stripping
%             df_dt_dm(i)=sqrt(virialF(i))*(0-par_dm(2)*(ft_dm(i)*Rvir(i)/r(i)))/2/pi;%/par_dm(1);  %tidal stripping
%            par_dm(1)=1;
            df_dt_dm(i)=sqrt(virialF(i))*(par_dm(1)-par_dm(2)*(ft_dm(i)*Rvir(i)/r(i)))/2/pi;  %tidal stripping
        case 1
%                 df_dt_dm(i)=sqrt(virialF(i))*ft_dm(i)*(par_dm(1)-par_dm(2)*(Rvir(i)/r(i))^(2/3))/2/pi;  %tidal stripping, with redistribution
% par_dm(1)=1;
            df_dt_dm(i)=sqrt(virialF(i))*ft_dm(i)*(1-par_dm(2)*Rvir(i)/r(i))/2/pi/par_dm(1);  %tidal stripping, with redistribution
            %       df_dt_dm(i)=sqrt(virialF(i))*ft_dm(i)*(1-2*Rvir(i)/r(i))/2/pi/par_dm(2);  %tidal stripping
        case 2    
%             df_dt_dm(i)=sqrt(virialF(i))*(sqrt(ft_dm(i))-ft_dm(i)*par_dm(2)*Rvir(i)/r(i))/2/pi/par_dm(1);  %tidal stripping, with redistribution
            df_dt_dm(i)=sqrt(virialF(i))*((ft_dm(i)*Rvir(i)/r(i))^2-ft_dm(i)*par_dm(2)*Rvir(i)/r(i))/2/pi/par_dm(1);  %tidal stripping, with redistribution
        case 3
%             par_dm=[1,2];
            df_dt_dm(i)=sqrt(virialF(i))*(r(i)/Rvir(i)/ft_dm(i)/par_dm(2)-1)/2/pi/par_dm(1);  %tidal stripping, 
        case 4
            df_dt_dm(i)=sqrt(virialF(i))*(r(i)/Rvir(i)/ft_dm(i)/par_dm(2)-ft_dm(i)*par_dm(2)*Rvir(i)/r(i))/2/pi/par_dm(1)/2;  %tidal stripping, 
        case 5
             df_dt_dm(i)=sqrt(virialF(i))*(par_dm(1)-(par_dm(2)-1/par_dm(2))*(ft_dm(i)*Rvir(i)/r(i))-(ft_dm(i)*Rvir(i)/r(i))^2)/2/pi;
        otherwise
            error('wrong flag_redistr');
    end
end
if flag_redistr==-3||flag_redistr==-4
    ft_dm=vc2(rt*csat)/vc2(csat).*rt;
end
frac_dm=mdm(1:Ninfall)/mdm(1);
weight=ones(Ninfall,1);
% weight=frac_dm;
% weight=1./mdm;
filter=logical(abs(ft_dm)<inf&(0==imag(ft_dm))&(~isnan(ft_dm)));mscale_dm=sum(frac_dm(filter).*ft_dm(filter).*weight(filter))/sum(frac_dm(filter).^2.*weight(filter));
if mscale_dm>1.01
    mscale_dm=1.01;
else if mscale_dm<0.99
        mscale_dm=0.99;
    end
end
% mscale_dm=1;
frac_dm(1)=1/mscale_dm;
xi2=sum((frac_dm*mscale_dm-ft_dm).^2.*weight)/sum(weight);

da=t-t(1);
switch flag_plot
    case 1
    figure;  
    plot(da,frac_dm*mscale_dm,'o');hold on;
    plot(da,ft_dm,'-');
    xlabel('\Delta ln(a)');ylabel('remained fraction');
    set(gca,'ylim',[0,1]);
    text(.05,.1,num2str([r(1:3)./Rvir(1:3);mscale_dm]));hold off;
    legend('DM fraction',['tidal strip model',num2str(sqrt(xi2)),sprintf('\n'),num2str(par_dm)],'location','best');
    case 2
    plot(frac_dm,ft_dm,'o','markersize',2,'markerfacecolor','k');hold on;
    case 3
    subplot(2,1,1);plot(da,frac_dm*mscale_dm,'o');hold on;
    subplot(2,1,2);plot(da,ft_dm,'.');hold on;
    otherwise
end
end