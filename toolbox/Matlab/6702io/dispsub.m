function dispsub(subcat,subid)
%function dispsub(subcat,subid)

subind=subid+1;

switch subcat.property.type
    case 'subcat'
    fprintf('\tsubid:\t%d\tSubLen:\t\t%d\n', subid,subcat.SubLen(subind));
    hostid=subcat.HaloChains.HostID(subind);
    if hostid<0 %quasi halo
        hostlen=subcat.NQuasi;
    else
        hostlen=subcat.GrpLen_Sub(hostid+1);
    end
    fprintf('\tRank:\t%d\tGrpLen_Sub:\t%d\tHostID:\t%d\tProSubID:\t%d\n',subcat.SubRank(subind),hostlen,...
        hostid,subcat.HaloChains.ProSubID(subind));
    fprintf('\tnibs:\t%d\tpre:\t\t%d\tnext:\t%d\tson:\t\t%d\n',subcat.sub_hierarchy.nibs(subind),subcat.sub_hierarchy.pre(subind),...
        subcat.sub_hierarchy.next(subind),subcat.sub_hierarchy.sub(subind));
    fprintf('\tCoM:\t%f,\t%f,\t%f\n',subcat.SubProp.CoM(subind,:));
    fprintf('\tVCoM:\t%f,\t%f,\t%f\n',subcat.SubProp.VCoM(subind,:));
    fprintf('\tKin:\t%f,\tPot:%f,\tVirRatio:%f\n',subcat.SubProp.Kin(subind),subcat.SubProp.Pot(subind),...
        -2*subcat.SubProp.Kin(subind)/subcat.SubProp.Pot(subind));
    fprintf('\tAM:\t%f,\t%f,\t%f\n',subcat.SubProp.AM(subind,:));
    
    case 'srccat'
        fprintf('\tsubid:\t%d\tSrcLen:\t%d\tSrcLen2:\t%d\tCoreFrac:\t%f\n', subid,subcat.SubLen(subind),subcat.SubLen2(subind),subcat.CoreFrac(subind));
    otherwise
        error('input SUBCAT or SRCCAT to display!')
end

