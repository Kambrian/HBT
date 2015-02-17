function par_gas=gas_stripping_rate_stat(hist_bin,pmass_bin,Snappar,Chostpar,Csatpar)

Ninfall_Min=10;Minfall_Min=1000;

Nhist=size(hist_bin,1);
par_gas=zeros(Nhist,3);%alpha,beta,xi
flag_plot=0;
flag_redistr=2;

for h=1:Nhist
    [Ninfall,t,r,v,vt2,kt,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,virsubs]=get_strp_history(hist_bin{h},pmass_bin(h),Snappar(h,2),1,1,2);
    if Ninfall>Ninfall_Min&&mdm(1)>Minfall_Min&&all(chost>0)&&Csatpar(h,2)>0&&virsubs(1)<1.2
        par=strip_coeff(Chostpar(h,2),Csatpar(h,2));
        [x,xi2]=fmincon(@(params) xi2_GASstrp_history(pmass_bin(h),Csatpar(h,2),r,v,t,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,[par(3)/params(1),par(4)/params(1)/params(2)],flag_plot,flag_redistr),[2,1],[],[],[],[],[1,0.1],[3,2]);
%        [x,xi2]=fmincon(@(params) xi2_GASstrp_history(pmass_bin(h),Csatpar(h,2),r,v,t,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,params,0,-1),[2,1],[],[],[],[],[1,0.1],[3,2]);
%         [x,xi2]=fmincon(@(params) xi2_GASstrp_history(pmass_bin(h),Csatpar(h,2),r,v,t,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,params,0,-2),[1.5,1],[],[],[],[],[0.5,0.1],[3,2]);
        if isreal(xi2)&&(~isnan(xi2))
            par_gas(h,:)=[x,sqrt(xi2)];
        end
    end
end

% alpha=0;
% % alpha=logspace(-6,1,200);
% betag=0:0.01:5;
% 
% xig=zeros(numel(alpha),numel(betag));
% for h=1:Nhist
%     [Ninfall,t,r,v,vt2,kt,mdm,mgas,mhost,fhost,chost,Hz,Rvir,virialF,virsubs]=get_strp_history(hist_bin{h},pmass_bin(h),Snappar(h,2),1,1,2);
%     if Ninfall>Ninfall_Min&&mgas(1)>Minfall_Min
% %         mrate=mdm./mhost;
% %         grate=mgas./mdm./fhost;
% %         f=mgas/mgas(1);y=diff(f)./diff(t)*2*pi;x=f./sqrt(grate).*mrate.^(-1/3).*v./r./Hz;x=(x(1:end-1)+x(2:end))/2;  
% %         c=sqrt(virialF);c=(c(1:end-1)+c(2:end))/2;y=y./c;x=x./c;
% %         [p,S]=polyfit(x,y,1);
% %         par_gas(h,:)=[p(2),-p(1),S.normr/sqrt(S.df)];
%         for i=1:numel(alpha)
%             for j=1:numel(betag)
%                     xig(i,j)=xi2_GASstrp_history(r,v,t,mdm,mgas,mhost,fhost,Hz,Rvir,virialF,[alpha(i),betag(j)],flag_plot,flag_redistr);
%             end
%         end
%         [ximin,indmin]=min(abs(xig(:)));
%         if isreal(xig(indmin))&&(~isnan(xig(indmin)))
%             [ig,jg]=ind2sub(size(xig),indmin);
%             par_gas(h,:)=[alpha(ig),betag(jg),sqrt(ximin)];
%         end
%     end
% end

% par=par_gas(:,3);
% figure;subplot(2,1,1);hist((par(logical(par))),20);xlabel('<\Delta f_{gas}>');ylabel('Counts');
% subplot(2,1,2);hist((par_dm(logical(par),3)),20);xlabel('<\Delta f_{dm}>');ylabel('Counts');
% x=linspace(0,0.4,40);
% figure;h1=subplot(2,1,1);h2=subplot(2,1,2);
% filter=logical(par>0&par_gas(:,4)<0.7);
%     y=histc(par_gas(filter,3),x);
%     subplot(h1);plot(x,y,'k- .');hold on;
%     subplot(h2);plot(x,cumsum(y),'k- .');hold on;
