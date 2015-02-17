function readheader(snapnum)
global HUBBLE0 G snapdir RUN_NUM 
global  header snaplist

snapfile=fullfile(snapdir,['pos',num2str(RUN_NUM),'.',num2str(snaplist(snapnum+1),'%04d')]);
%
% typedef struct
% {
%    int Np;   //number of particles in the simulation
%    int ips;  //time step number, or snapshot num
%    float ztp; // current redshift
%    float Omegat;  // current omega_m, mass density
%    float Lambdat; // current Omega_lambda, dark energy density
%    float rLbox;  // boxsize in unit of Mpc/h
%    float xscale;
%    float vscale;
%    float Hz; //Hubble Param at ztp
%    float vunit; //velocity unit to get physical peculiar velocity in km/s
%    /*==extension needed for calculating binding energy:==*/
%    float time;//current reduced scale factor 
%    float mass[2];//Gas (mass[0]) and DM (mass[1]) particle masses, in units of 10^10Msun/h
% }IO_HEADER; //header of jing's data structure


fid=fopen(snapfile,'r','ieee-be');
dummy=fread(fid,1,'int32');
header.Np=fread(fid,1,'int32');
header.ips=fread(fid,1,'int32');
header.z=fread(fid,1,'float32');
header.Omegat=fread(fid,1,'float32');
header.Lambdat=fread(fid,1,'float32');
header.rLbox=fread(fid,1,'float32');
header.xscale=fread(fid,1,'float32');
header.vscale=fread(fid,1,'float32');

dummy2=fread(fid,1,'int32');
if dummy~=dummy2
    error(['error reading header file:',snapfile,';\nbrackets mismatch:',num2str(dummy),',',num2str(dummy2)]);
end
fclose(fid);

OMEGAL0=0.732;OMEGA0=0.268;BOXSIZE=150000.0;Redshift_INI=144.0;Hratio=sqrt(OMEGAL0/header.Lambdat);scale_reduced=1./(1.+header.z);scale0=1+Redshift_INI;
header.Hz=HUBBLE0*Hratio;
header.time=scale_reduced;
header.mass(1)=0;
header.mass(2)=OMEGA0*3*HUBBLE0*HUBBLE0/8/pi/G*BOXSIZE*BOXSIZE*BOXSIZE/header.Np;
header.vunit=100.*header.rLbox*Hratio*scale_reduced*scale_reduced*scale0;
    
header.snap=snapnum;