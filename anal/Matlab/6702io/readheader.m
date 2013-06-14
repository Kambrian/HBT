function readheader(snapnum)
global HUBBLE0 snapdir
global  header

snapfile=fullfile(snapdir,['snapshot_',num2str(snapnum,'%03d')]);
%
% typedef	struct 
% {
%   int      npart[6];
%   double   mass[6];
%   double   time;
%   double   redshift;
%   int      flag_sfr;
%   int      flag_feedback;
%   int      npartTotal[6];
%   int      flag_cooling;
%   int      num_files;
%   double   BoxSize;
%   double   Omega0;
%   double   OmegaLambda;
%   double   HubbleParam; 
%   double Hz;//current Hubble param in internal units, BT extension
%   char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8 -8];  /* fills to 256 Bytes */
% } IO_HEADER;
%
fid=fopen(snapfile,'r');
dummy=fread(fid,1,'int32');
header.npart=fread(fid,6,'int32');
header.mass=fread(fid,6,'double');
header.time=fread(fid,1,'double');
header.z=fread(fid,1,'double');
header.flag_sfr=fread(fid,1,'int32');
header.flag_feedback=fread(fid,1,'int32');
header.nall=fread(fid,6,'int32');
header.flag_cooling=fread(fid,1,'int32');
header.num_files=fread(fid,1,'int32');
header.boxsize=fread(fid,1,'double');
header.Omega0=fread(fid,1,'double');
header.OmegaL=fread(fid,1,'double');
header.h=fread(fid,1,'double');
fseek(fid,256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8,'cof');
dummy2=fread(fid,1,'int32');
if dummy~=dummy2
    error(['error reading header file:',snapfile,';\nbrackets mismatch:',num2str(dummy),',',num2str(dummy2)]);
end
fclose(fid);
header.Hz=HUBBLE0*sqrt(header.Omega0/(header.time^3)+(1-header.Omega0-header.OmegaL)/(header.time ^2)+header.OmegaL);%Hubble param for the current catalogue;
header.snap=snapnum;