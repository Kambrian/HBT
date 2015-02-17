function load_particle_data(snapnum,blocks)
% function load_particle_data(snapnum,blocks)
%        blocks='pvi','vi','pi','i','v'... 
%                (the first char of each block to be read) 
%                indicating whether each block should be read or not 
%               if omitted, set to 'pvi' by default
if nargin<2, blocks='ipv';end
blocks=lower(blocks);
skip_pos=isempty(strfind(blocks,'p'));
skip_vel=isempty(strfind(blocks,'v'));
skip_id=isempty(strfind(blocks,'i'));

global snapdir
global Pdat header


if header.snap~=snapnum %read header if haven't
    readheader(snapnum);
end

snapfile=fullfile(snapdir,['snapshot_',num2str(snapnum,'%03d')]);
fid=fopen(snapfile,'r');
fseek(fid,256+8,'bof');%skip header

pre_len=header.npart(1)*4*3;
tail_len=sum(header.npart(3:6),1)*4*3;


dummy=fread(fid,1,'int32');
fseek(fid,pre_len,'cof');
if skip_pos 
    fseek(fid,header.npart(2)*3*4,'cof');
else
    disp('reading particle positions...')
    Pdat.Pos=fread(fid,[3,header.npart(2)],'float32');
    Pdat.Pos=Pdat.Pos';
end
fseek(fid,tail_len,'cof');
dummy2=fread(fid,1,'int32');
if dummy~=dummy2
    error(['error reading particle pos:',snapfile,';\nbrackets mismatch:',num2str(dummy),',',num2str(dummy2)]);
end


dummy=fread(fid,1,'int32');
fseek(fid,pre_len,'cof');
if skip_vel
    fseek(fid,header.npart(2)*3*4,'cof');
else
    disp('reading particle velocities...')
    Pdat.Vel=fread(fid,[3,header.npart(2)],'float32');
    Pdat.Vel=Pdat.Vel';
end
fseek(fid,tail_len,'cof');
dummy2=fread(fid,1,'int32');
if dummy~=dummy2
    error(['error reading particle vel:',snapfile,';\nbrackets mismatch:',num2str(dummy),',',num2str(dummy2)]);
end

pre_len=header.npart(1)*4;
tail_len=sum(header.npart(3:6),1)*4;


dummy=fread(fid,1,'int32');
fseek(fid,pre_len,'cof');
if skip_id
    fseek(fid,header.npart(2)*4,'cof');
else
    disp('reading particle ids...')
    Pdat.PID=fread(fid,header.npart(2),'int32');
end
fseek(fid,tail_len,'cof');
dummy2=fread(fid,1,'int32');
if dummy~=dummy2
    error(['error reading particle ids:',snapfile,';\nbrackets mismatch:',num2str(dummy),',',num2str(dummy2)]);
end

fclose(fid);

Pdat.snap=snapnum;
