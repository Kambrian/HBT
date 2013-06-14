function halo=get_halo_mat(RunNum)
halo=cell(100,1);
for nsnap=0:99
halo{nsnap+1}=readhalo_size(['/mnt/A4700/data/',num2str(RunNum),'/subcat/profile/logbin'],nsnap,'halo');
end
save(['/home/kam/Documents/research/Galaxy/code/BoundTracing/data/halo',num2str(RunNum),'.mat'],'halo');
