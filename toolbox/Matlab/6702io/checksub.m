BTinit
snapnum=44;
load_particle_data(snapnum,'ip');
subcat=load_sub_catalogue(snapnum);
srccat=load_src_catalogue(snapnum);
[cat,subcat,srccat]=fresh_id2index(cat,subcat,srccat);

%subid=4430;
subid=find(subcat.HaloChains.ProSubID==4430)
grpid=subcat.HaloChains.HostID(subid+1)
dispsub(subcat,subid)
subarr=subcat.PSubArr{subid+1};
grparr=cat.PIDorIndex(cat.Offset(grpid+1)+(1:cat.Len(grpid+1)));
figure;plot3(Pdat.Pos(subarr,1),Pdat.Pos(subarr,2),Pdat.Pos(subarr,3),'o');
figure;plot3(Pdat.Pos(grparr,1),Pdat.Pos(grparr,2),Pdat.Pos(grparr,3),'o');