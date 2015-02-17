function halo=readhalo_size(basedir,Nsnap,halotype)

sizefile=fullfile(basedir,[halotype,'_size_',num2str(Nsnap,'%03d')])

fid=fopen(sizefile,'r');
if isempty(fread(fid)) 
    halo.ngrps=0;
    fclose(fid);
    return;
else
    frewind(fid);
end
switch halotype
    case 'halo'
        %   int nbin;
        % 	float rmax;
        % 	float rvir;
        % 	float req_bk_1;//set to 0 when no bin,to rmax when no back
        % 	float req_all_1;//set to 0 when no bin,to rmax when no back+other
        % 	float req_bk_02;
        % 	float req_all_02;
        % 	float CoM[3];//halo CoM,for comparison
        % 	float Cen[3];//Center of mainsub,the profiles are w.r.t this center
        % 	int mass;//fof mass
        % 	int Mvir[3];
        % 	float Rvir[3];//[tophat,c200,b200]
        % 	int flag_badvir[3];
        %   int flag_fakehalo;
        skip=(24-1)*4;
        halo.nbin=fread(fid,'int32',skip);
        fseek(fid,1*4,'bof');
        halo.rmax=fread(fid,'float32',skip);
        fseek(fid,2*4,'bof');
        halo.rvir=fread(fid,'float32',skip);
        fseek(fid,3*4,'bof');
        halo.req_bk_1=fread(fid,'float32',skip);
        fseek(fid,4*4,'bof');
        halo.req_all_1=fread(fid,'float32',skip);
        fseek(fid,5*4,'bof');
        halo.req_bk_02=fread(fid,'float32',skip);
        fseek(fid,6*4,'bof');
        halo.req_all_02=fread(fid,'float32',skip);
        fseek(fid,7*4,'bof');
        halo.CoM(:,1)=fread(fid,'float32',skip);
        fseek(fid,8*4,'bof');
        halo.CoM(:,2)=fread(fid,'float32',skip);
        fseek(fid,9*4,'bof');
        halo.CoM(:,3)=fread(fid,'float32',skip);
        fseek(fid,10*4,'bof');
        halo.Cen(:,1)=fread(fid,'float32',skip);
        fseek(fid,11*4,'bof');
        halo.Cen(:,2)=fread(fid,'float32',skip);
        fseek(fid,12*4,'bof');
        halo.Cen(:,3)=fread(fid,'float32',skip);
        fseek(fid,13*4,'bof');
        halo.mass=fread(fid,'int32',skip);
        fseek(fid,14*4,'bof');
        halo.Mvir(:,1)=fread(fid,'int32',skip);
        fseek(fid,15*4,'bof');
        halo.Mvir(:,2)=fread(fid,'int32',skip);
        fseek(fid,16*4,'bof');
        halo.Mvir(:,3)=fread(fid,'int32',skip);
        fseek(fid,17*4,'bof');
        halo.Rvir(:,1)=fread(fid,'float32',skip);
        fseek(fid,18*4,'bof');
        halo.Rvir(:,2)=fread(fid,'float32',skip);
        fseek(fid,19*4,'bof');
        halo.Rvir(:,3)=fread(fid,'float32',skip);
        fseek(fid,20*4,'bof');
        halo.flag_badvir(:,1)=fread(fid,'int32',skip);
        fseek(fid,21*4,'bof');
        halo.flag_badvir(:,2)=fread(fid,'int32',skip);
        fseek(fid,22*4,'bof');
        halo.flag_badvir(:,3)=fread(fid,'int32',skip);
        fseek(fid,23*4,'bof');
        halo.flag_fakehalo=fread(fid,'int32',skip);
    case 'ghalo'
        %   int nbin;
        % 	float rmax;
        % 	float rvir;//estimated vir for DMhalo
        % 	float req_bk_1;//set to 0 when no bin,to rmax when no back
        % 	float req_all_1;//set to 0 when no bin,to rmax when no back+other
        % 	float req_bk_02;
        % 	float req_all_02;
        % 	float CoM[3];
        % 	int Mvir_this[3];//halo gas mass inside virial radii
        % 	int Mvir_all[3];//all gas mass inside virial radii
        % 	int mass;
        skip=(17-1)*4;
        halo.nbin=fread(fid,'int32',skip);
        fseek(fid,1*4,'bof');
        halo.rmax=fread(fid,'float32',skip);
        fseek(fid,2*4,'bof');
        halo.rvir=fread(fid,'float32',skip);
        fseek(fid,3*4,'bof');
        halo.req_bk_1=fread(fid,'float32',skip);
        fseek(fid,4*4,'bof');
        halo.req_all_1=fread(fid,'float32',skip);
        fseek(fid,5*4,'bof');
        halo.req_bk_02=fread(fid,'float32',skip);
        fseek(fid,6*4,'bof');
        halo.req_all_02=fread(fid,'float32',skip);
        fseek(fid,7*4,'bof');
        halo.CoM(:,1)=fread(fid,'float32',skip);
        fseek(fid,8*4,'bof');
        halo.CoM(:,2)=fread(fid,'float32',skip);
        fseek(fid,9*4,'bof');
        halo.CoM(:,3)=fread(fid,'float32',skip);
        fseek(fid,10*4,'bof');
        halo.Mvir_this(:,1)=fread(fid,'int32',skip);
        fseek(fid,11*4,'bof');
        halo.Mvir_this(:,2)=fread(fid,'int32',skip);
        fseek(fid,12*4,'bof');
        halo.Mvir_this(:,3)=fread(fid,'int32',skip);
        fseek(fid,13*4,'bof');
        halo.Mvir_all(:,1)=fread(fid,'int32',skip);
        fseek(fid,14*4,'bof');
        halo.Mvir_all(:,2)=fread(fid,'int32',skip);
        fseek(fid,15*4,'bof');
        halo.Mvir_all(:,3)=fread(fid,'int32',skip);
        fseek(fid,16*4,'bof');
        halo.mass=fread(fid,'int32',skip);
    otherwise
        error('wrong halotype, must be "halo" or "ghalo"');
end
fclose(fid);
halo.ngrps=numel(halo.nbin);
