function A=load_abundance(file)
fid=fopen(file);
A.grpids=fread(fid,inf,'int32',8*4);
fseek(fid,1*4,'bof');
A.Ns=fread(fid,inf,'int32',8*4);
fseek(fid,2*4,'bof');
A.status=fread(fid,inf,'int32',8*4);
fseek(fid,3*4,'bof');
A.rhos=fread(fid,inf,'float32',8*4);
fseek(fid,4*4,'bof');
A.rs=fread(fid,inf,'float32',8*4);
fseek(fid,5*4,'bof');
A.c=fread(fid,inf,'float32',8*4);
fseek(fid,6*4,'bof');
A.dM=fread(fid,inf,'float32',8*4);
fseek(fid,7*4,'bof');
A.Mvir=fread(fid,inf,'float32',8*4);
fseek(fid,8*4,'bof');
A.Rvir=fread(fid,inf,'float32',8*4);
fclose(fid);

