function haloprof=readhalo_prof(basedir,Nsnap,halotype,halo)

proffile=fullfile(basedir,[halotype,'_prof_',num2str(Nsnap,'%03d')])

nbin=halo.nbin;
fid2=fopen(proffile,'r');
for grpind=1:halo.ngrps;
haloprof.ns{grpind}=fread(fid2,nbin(grpind),'int32');
haloprof.no{grpind}=fread(fid2,nbin(grpind),'int32');
haloprof.nb{grpind}=fread(fid2,nbin(grpind),'int32');
haloprof.rs{grpind}=fread(fid2,nbin(grpind),'float32');
haloprof.ro{grpind}=fread(fid2,nbin(grpind),'float32');
haloprof.rb{grpind}=fread(fid2,nbin(grpind),'float32');
end
fclose(fid2);
