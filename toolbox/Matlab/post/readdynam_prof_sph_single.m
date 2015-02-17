function dynamprof=readdynam_prof_sph_single(basedir,Nsnap,halo,grpind)

proffile=fullfile(basedir,['dynam_prof_',num2str(Nsnap,'%03d'),'.sph'])

nbin=halo.nbin;
fid2=fopen(proffile,'r');
fseek(fid2,sum(nbin(1:grpind-1))*25*4,'bof');
dynamprof.vm_xyz=fread(fid2,[3,nbin(grpind)],'float32');
dynamprof.vd_xyz=fread(fid2,[3,nbin(grpind)],'float32');
dynamprof.vm_rtf=fread(fid2,[3,nbin(grpind)],'float32');
dynamprof.vd_rtf=fread(fid2,[3,nbin(grpind)],'float32');
dynamprof.Um=fread(fid2,nbin(grpind),'float32');
dynamprof.gvm_xyz=fread(fid2,[3,nbin(grpind)],'float32');
dynamprof.gvd_xyz=fread(fid2,[3,nbin(grpind)],'float32');
dynamprof.gvm_rtf=fread(fid2,[3,nbin(grpind)],'float32');
dynamprof.gvd_rtf=fread(fid2,[3,nbin(grpind)],'float32');
fclose(fid2);
