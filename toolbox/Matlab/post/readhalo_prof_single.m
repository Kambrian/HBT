function haloprof=readhalo_prof_single(basedir,Nsnap,halotype,halo,grpind)

proffile=fullfile(basedir,[halotype,'_prof_',num2str(Nsnap,'%03d')])

nbin=halo.nbin;
fid2=fopen(proffile,'r');
fseek(fid2,sum(nbin(1:grpind-1))*6*4,'bof');
haloprof.ns=fread(fid2,nbin(grpind),'int32');
haloprof.no=fread(fid2,nbin(grpind),'int32');
haloprof.nb=fread(fid2,nbin(grpind),'int32');
haloprof.rs=fread(fid2,nbin(grpind),'float32');
haloprof.ro=fread(fid2,nbin(grpind),'float32');
haloprof.rb=fread(fid2,nbin(grpind),'float32');
fclose(fid2);
