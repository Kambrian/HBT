function par_dm=DM_stripping_rate_stat(hist_bin,pmass_bin,Snappar,Chostpar,Csatpar)

Ninfall_Min=10;Minfall_Min=1000;

Nhist=size(hist_bin,1);

par_dm=zeros(Nhist,3);%alpha,beta,xi2,Kself

flag_plot=0;
flag_redistr=0;

for h=1:Nhist
    [Ninfall,t,r,v,vt2,kt,mdm,mhost,chost,Hz,Rvir,virialF,virsubs]=get_DMstrp_history(hist_bin{h},pmass_bin(h),Snappar(h,2),1,1,1);
    if Ninfall>Ninfall_Min&&mdm(1)>Minfall_Min&&all(chost>0)&&Csatpar(h,2)>0&&virsubs(1)<1.2
        kt=ones(size(t));
        par=strip_coeff(Chostpar(h,2),Csatpar(h,2));
[x,xi2]=fmincon(@(params) xi2_DMstrp_history(Hz,Csatpar(h,2),r,vt2,kt,t,mdm,chost,Rvir,virialF,[par(1)/params(1),par(2)/params(1)/params(2)],0,flag_redistr),[1,0.5],[],[],[],[],[0.01,0.01],[10,1.5]);
% [x,xi2]=fmincon(@(params) xi2_DMstrp_history(Hz,Csatpar(h,2),r,vt2,kt,t,mdm,chost,Rvir,virialF,params,0,-5),[1,0.4],[],[],[],[],[0.01,0.01],[10,1.5]);
% [x,xi2]=fmincon(@(params) xi2_DMstrp_history(Hz,Csatpar(h,2),r,vt2,kt,t,mdm,chost,Rvir,virialF,params,0,-6),[0.6,0.4],[],[],[],[],[0.01,0.01],[10,1.5]);
        if isreal(xi2)&&(~isnan(xi2))
            par_dm(h,:)=[x,sqrt(xi2)];
        end
    end
end

% 
% alpha=0;
% betad=0:0.01:5;
% 
% xid=zeros(numel(alpha),numel(betad));
% for h=1:Nhist
%     [Ninfall,t,r,v,kt,mdm,mgas,mhost,fhost,chost,Hz,Rvir,virialF,virsubs]=get_strp_history(hist_bin{h},pmass_bin(h),Snappar(h,2),1,1,1);
%     if Ninfall>Ninfall_Min&&mdm(1)>Minfall_Min&&all(chost>0)
% %         f=mdm/mdm(1);y=diff(f)./diff(t)*2*pi./sqrt(virialF(1:end-1));x=f.*Rvir./r;x=(x(1:end-1)+x(2:end))/2;  
% %         [p,S]=polyfit(x,y,1);
% %         par_dm(h,:)=[p(2),-p(1),S.normr/sqrt(S.df)];
%         for i=1:numel(alpha)
%             for j=1:numel(betad)
%                 xid(i,j)=xi2_DMstrp_history(r,kt,t,mdm,chost,Rvir,virialF,[alpha(i),betad(j)],flag_plot,flag_redistr);
%             end
%         end
%         [ximin,indmin]=min(abs(xid(:)));
%         if isreal(xid(indmin))&&(~isnan(xid(indmin)))
%             [id,jd]=ind2sub(size(xid),indmin);
%             par_dm(h,:)=[alpha(id),betad(jd),sqrt(ximin)];
%         end
%     end
% end

% pard=par_dm(:,3);
% figure;subplot(2,1,1);hist((par(logical(par))),20);xlabel('<\Delta f_{gas}>');ylabel('Counts');
% subplot(2,1,2);hist((par_dm(logical(par),3)),20);xlabel('<\Delta f_{dm}>');ylabel('Counts');
% x=linspace(0,0.4,40);
% figure;h1=subplot(2,1,1);h2=subplot(2,1,2);
% filter=logical(pard>0&pard<0.05&par_dm(:,4)<0.7);
% y=histc(par_dm(filter,3),x);
% subplot(h1);plot(x,y,'k- .');hold on;
% subplot(h2);plot(x,cumsum(y),'k- .');
% figure;subplot(2,1,1);hist(par_dm(filter,1),numel(alpha));subplot(2,1,2);hist(par_dm(filter,2),numel(betad));