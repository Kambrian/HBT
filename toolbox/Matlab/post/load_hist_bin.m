function [hist_bin,pmass]=load_hist_bin(filename)

fid=fopen(filename,'r');
pmass=fread(fid,1,'float32');% dark matter particle mass
gas_flag=fread(fid,1,'int32');
Nhist=fread(fid,1,'int32');
hist_bin=cell(Nhist,1);
for i=1:Nhist
%     hist_bin{i}.SID=fread(fid,1,'int32');%infall subid
    hist_bin{i}.Nnode=fread(fid,1,'int32');
    for j=1:hist_bin{i}.Nnode
        hist_bin{i}.node(j).Nsnap=fread(fid,1,'int32');
        hist_bin{i}.node(j).HostID=fread(fid,1,'int32');
        hist_bin{i}.node(j).subid=fread(fid,1,'int32');
        if gas_flag==1
        hist_bin{i}.node(j).mass=fread(fid,6,'int32');
        hist_bin{i}.node(j).Umsat=fread(fid,1,'float32');
        else
        hist_bin{i}.node(j).mass=fread(fid,4,'int32');    
        end
        hist_bin{i}.node(j).pos=fread(fid,3,'float32');
        hist_bin{i}.node(j).vel=fread(fid,3,'float32');
        hist_bin{i}.node(j).chost=fread(fid,1,'float32');
%         hist_bin{i}.node(j).pot=fread(fid,1,'float32');
%         hist_bin{i}.node(j).kin=fread(fid,1,'float32');
%         hist_bin{i}.node(j).AM=fread(fid,3,'float32');
%         hist_bin{i}.node(j).Hpot=fread(fid,1,'float32');
%         hist_bin{i}.node(j).Hkin=fread(fid,1,'float32');
%         hist_bin{i}.node(j).HAM=fread(fid,3,'float32');
    end
end
Nhist2=fread(fid,1,'int32');
if(Nhist~=Nhist2)
    error(['reading error: Nhist=',Nhist,',or ',Nhist2]);
end

        
      