function [Msubs,Mhost,Rvir,Rsub,IDsub]=Msublist_in_radii(datadir,virtype,Nsnap,grpid,rmin,rmax,datatype,plottype)
% list subhalo masses inside radius [rmin*Rvir,rmax*Rvir]
% return Msubs and host virial mass in units 10^10Msun/h
% also Rsub

switch datatype
    case 'textdata'
%         datadir=['/mnt/A4700/data/',RunName,'/subcat/anal/massfun/'];
        massdata=load([datadir,'sub_mass_',num2str(Nsnap,'%03d')]);
        groupdata=load([datadir,'group_offset_',num2str(Nsnap,'%03d')]);
        halo=readhalo_size(['/mnt/A4700/data/',num2str(RunNum),'/subcat/profile/logbin'],Nsnap,'halo');

        partmass=groupdata(grpid+1,3)/halo.mass(grpid+1);
        Mhost=halo.Mvir(grpid+1,virtype+1)*partmass;
        Rvir=halo.Rvir(grpid+1,virtype+1);           %almost all subs within Rvir are within the
        %fof halo rather than from other halos
        center_id=groupdata(grpid+1,1)+1;
        Rsub=sqrt(sum((massdata(:,2:4)-repmat(massdata(center_id,2:4),size(massdata,1),1)).^2,2));
        massdata(center_id,:)=[];Rsub(center_id,:)=[];%exclude central sub
        Msubs=massdata(Rsub>rmin*Rvir&Rsub<rmax*Rvir,1);
        Rsub=Rsub(Rsub>rmin*Rvir&Rsub<rmax*Rvir);
    case 'bindata'
%         datadir=['/mnt/A4700/data/',RunName,'/subcat/anal/massfun/'];
          submass=read_submass([datadir,'submass_',num2str(Nsnap,'%03d')]);
%         tmp=importdata('/mnt/charon/HBT/data/AqA2/subcatmore/anal/steller/SnapInfall_162',',',1);%replace with infall mass
%         submass=tmp.data(:,1);%./tmp.data(:,1);
%         clear tmp
        subcom=read_subcom([datadir,'subcom_',num2str(Nsnap,'%03d')]);
        massdata=[submass,subcom];
        clear submass subcom
        cenid=read_cid([datadir,'cid_',num2str(Nsnap,'%03d')]);
        switch virtype
            case 0
                grpsize=read_grpsize([datadir,'grpsizeVIR_',num2str(Nsnap,'%03d')]);
            case 1
                grpsize=read_grpsize([datadir,'grpsizeC200_',num2str(Nsnap,'%03d')]);
            case 2
                grpsize=read_grpsize([datadir,'grpsizeB200_',num2str(Nsnap,'%03d')]);
        end
        Mhost=grpsize(grpid+1,1);
        Rvir=grpsize(grpid+1,2);
        clear grpsize
        center_id=cenid(grpid+1)+1;
        Rsub=sqrt(sum((massdata(:,2:4)-repmat(massdata(center_id,2:4),size(massdata,1),1)).^2,2));
        massdata(center_id,:)=[];Rsub(center_id,:)=[];%exclude central sub
        Msubs=massdata(Rsub>rmin*Rvir&Rsub<rmax*Rvir,1);
        if nargout>4
            IDsub=0:size(massdata,1);
            IDsub(center_id)=[];
            IDsub=IDsub(Rsub>rmin*Rvir&Rsub<rmax*Rvir);
        end
        Rsub=Rsub(Rsub>rmin*Rvir&Rsub<rmax*Rvir);
%         Msubs(Msubs<5*min(Msubs))=[];%exclude subs with less than 100 particles
    otherwise
        error('wrong datafile type');
end
switch plottype
    case 'norm'
%uncomment for normalized mass function:
         Msubs=Msubs/Mhost;
    otherwise
end