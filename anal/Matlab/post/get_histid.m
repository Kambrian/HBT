function h=get_histid(snap,subid,hist_bin)
% h=get_histid(snap,subid,hist_bin)
% return the hist_id for the subhalo whose history have subid at snap

% subid=2291;
% snap=70;
Nhist=numel(hist_bin);
for h=1:Nhist
    Nnode=hist_bin{h}.Nnode;
    for i=1:Nnode
        if snap==hist_bin{h}.node(i).Nsnap
            if subid==hist_bin{h}.node(i).subid
                return;
            end
        end
    end
end
h=-1; %not found

