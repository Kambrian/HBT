global G HUBBLE0 snapdir fofdir subcatdir outputdir
HUBBLE0=0.1;
G=43007.1;
snapdir='/home/kambrain/data/6702/simu'
fofdir='/home/kambrain/data/6702/fof'
subcatdir='/home/kambrain/data/6702/subcat'
outputdir=[subcatdir,'/anal']

global Pdat header;
Pdat=struct('PID',[],'Pos',[],'Vel',[],'snap',-1);
header=struct('npart',zeros(6,1),'mass',zeros(6,1),'time',0,'z',0,'flag_sfr','0','flag_feedback',0,...
    'nall',zeros(6,1),'flag_cooling',0,'num_files',0,'boxsize',0,'Omega0',0,'OmegaL',0,'h',0,'Hz',0,'snap',-1);

global rtidal;
rtidal=[];