function J=bin_image(I, binsz)
% to bin the image I with size binsz

if numel(binsz)==1
    binsz=repmat(binsz,1,ndims(I));
end

nx=size(I,1);
ny=size(I,2);

mx=nx/binsz(1);
my=ny/binsz(2);

J=zeros(mx,my);

for i=1:mx
    for j=1:my
        c=I((i-1)*binsz(1)+1:i*binsz(1),(j-1)*binsz(2)+1:j*binsz(2));
        J(i,j)=sum(c(:));
    end
end
