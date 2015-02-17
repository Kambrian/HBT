function dynamprof=readdynam_prof_sph(basedir,Nsnap,halo)

proffile=fullfile(basedir,['dynam_prof_',num2str(Nsnap,'%03d'),'.sph'])

nbin=halo.nbin;
fid2=fopen(proffile,'r');
for grpind=1:halo.ngrps;
dynamprof.vm_xyz{grpind}=fread(fid2,[3,nbin(grpind)],'float32');
dynamprof.vd_xyz{grpind}=fread(fid2,[3,nbin(grpind)],'float32');
dynamprof.vm_rtf{grpind}=fread(fid2,[3,nbin(grpind)],'float32');
dynamprof.vd_rtf{grpind}=fread(fid2,[3,nbin(grpind)],'float32');
dynamprof.Um{grpind}=fread(fid2,nbin(grpind),'float32');
dynamprof.gvm_xyz{grpind}=fread(fid2,[3,nbin(grpind)],'float32');
dynamprof.gvd_xyz{grpind}=fread(fid2,[3,nbin(grpind)],'float32');
dynamprof.gvm_rtf{grpind}=fread(fid2,[3,nbin(grpind)],'float32');
dynamprof.gvd_rtf{grpind}=fread(fid2,[3,nbin(grpind)],'float32');
end
fclose(fid2);
