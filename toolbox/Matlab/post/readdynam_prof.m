function dynamprof=readdynam_prof(basedir,Nsnap,halo)

proffile=fullfile(basedir,['dynam_prof_',num2str(Nsnap,'%03d')])

nbin=halo.nbin;
fid2=fopen(proffile,'r');
for grpind=1:halo.ngrps;
dynamprof.vm2s{grpind}=fread(fid2,nbin(grpind),'float32');
dynamprof.vm2o{grpind}=fread(fid2,nbin(grpind),'float32');
dynamprof.vm2b{grpind}=fread(fid2,nbin(grpind),'float32');
dynamprof.vsms{grpind}=fread(fid2,nbin(grpind),'float32');
dynamprof.vsmo{grpind}=fread(fid2,nbin(grpind),'float32');
dynamprof.vsmb{grpind}=fread(fid2,nbin(grpind),'float32');
dynamprof.Ums{grpind}=fread(fid2,nbin(grpind),'float32');
dynamprof.Umo{grpind}=fread(fid2,nbin(grpind),'float32');
dynamprof.Umb{grpind}=fread(fid2,nbin(grpind),'float32');
dynamprof.gvm2s{grpind}=fread(fid2,nbin(grpind),'float32');
dynamprof.gvm2o{grpind}=fread(fid2,nbin(grpind),'float32');
dynamprof.gvm2b{grpind}=fread(fid2,nbin(grpind),'float32');
dynamprof.gvsms{grpind}=fread(fid2,nbin(grpind),'float32');
dynamprof.gvsmo{grpind}=fread(fid2,nbin(grpind),'float32');
dynamprof.gvsmb{grpind}=fread(fid2,nbin(grpind),'float32');
end
fclose(fid2);
