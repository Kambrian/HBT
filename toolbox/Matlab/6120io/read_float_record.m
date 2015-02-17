function arr=read_float_record(nel,fid)

flag_subrec=0;
s=zeros(4);
arr=zeros(1:nel);
dummy=fread(fid,1,'int32');
if dummy<0
    dummy=-1*dummy;
    flag_subrec=1;
end
nread=floor(dummy/4);
nrem=mod(dummy,4);
arr(1:nread)=fread(fid,nread,'float32');
offset=nread;
if nrem~=0
    s(1:nrem)=fread(fid,nrem,'uint8');
end
dummy2=fread(fid,1,'int32');
if dummy2<0
    error('subrecord reading error: wrong first trailing rec');
end
check_record_len(dummy,dummy2);
[filename, permission, machineformat] = fopen(fid);%get information about the file being read
while flag_subrec
    flag_subrec=0;
    dummy=fread(fid,1,'int32');
    if dummy<0
        dummy=-1*dummy;
        flag_subrec=1;
    end
    if nrem~=0
        s(nrem+1:4)=fread(fid,4-nrem,'unit8');
        tmpfile=fopen('_s.tmp','w');
        fwrite(tmpfile,s);
        fclose(tmpfile);
        tmpfile=fopen('_s.tmp','r',machineformat);
        arr(offset+1)=fread(tmpfile,1,'float32');
        fclose(tmpfile);
        offset=offset+1;
        nread=floor((dummy-(4-nrem))/4);
        nrem=mod(dummy-(4-nrem),4);
    else
        nread=floor(dummy/4);
        nrem=mod(dummy,4);
    end
    arr(offset+(1:nread))=fread(fid,nread,'float32');
    offset=offset+nread;
    if nrem~=0
        s(1:nrem)=fread(fid,nrem,'uint8');
    end
    dummy2=fread(fid,1,'int32');
    if dummy2>0
        error('subrecord reading error: wrong trailing rec');
    end
    dummy2=-1*dummy2;
    check_record_len(dummy,dummy2);
end
if offset~=nel
    error(['error reading record:wanted=',num2str(nel),'read=',num2str(offset)]);
end
if nrem~=0
    error(['error reading record:record does not contain integer number of elements',num2str(nrem)]);
end
   



function check_record_len(len1,len2)
if len1~=len2
    error(['error,record brackets mismatch:',num2str(len1),',',num2str(len2)]);
end