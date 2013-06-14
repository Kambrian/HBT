function [cat1_new,cat2_new,cat3_new]=fresh_id2index(cat1,cat2,cat3)
%[cat1_new,cat2_new,cat3_new]=fresh_id2index(cat1,cat2,..)

global header Pdat
error(nargchk(1, 3, nargin));
if nargout~=nargin, error('input-output mismatch');end
if isempty(Pdat.PID), error('load particle ids first!'); end

global PInd  %seems octave does not support nested function's sharing of var...
              %use global as a work-around
NP=sum(header.npart);
PInd=zeros(NP,1);
input('making indexes for query:?');
PInd(Pdat.PID)=(1:header.npart(2))';

input('freshing indexes:?');
cat1_new=id2ind(cat1);
if nargin>1, cat2_new=id2ind(cat2); end
if nargin>2, cat3_new=id2ind(cat3); end

    function cat_new=id2ind(cat)
        global PInd
        if isstruct(cat)
            if strcmp(cat.property.particles,'id')
                switch cat.property.type
                    case 'cat'
                        cat.PIDorIndex=PInd(cat.PIDorIndex);
                    case 'subcat'
                        cat.PSubArr=mat2cell(PInd(cell2mat(cat.PSubArr)),cat.SubLen,1);
                    case 'srccat'
                        cat.PSubArr=mat2cell(PInd(cell2mat(cat.PSubArr)),cat.SubLen,1);
                        sublen2=cat.SubLen2;
                        sublen2(sublen2<0)=0;
                        if sum(sublen2)>0
                        cat.PSubArr2=mat2cell(PInd(cell2mat(cat.PSubArr2)),sublen2,1);
                        end
                    otherwise
                        error('unknown catalogue type');
                end
                cat.property.particles='matlab-index';
            else
                warning('cat already being index, no need to refresh') %#ok<WNTAG>
            end
        else
            cat(:)=PInd(cat(:));
        end
        cat_new=cat;
    end
clear PInd
end