function a=load_scaleF(scaleF_file)
% scaleF_file=['/mnt/A4700/data/',runnum,'/subcat/Redshift.dat'];
% scaleF_file=['/mnt/charon/HBT/data/',RunName,'/subcat/Redshift.dat'];
tmp=load(scaleF_file);a=tmp(:,2);    
