function [perm] = permute_corr(Orig,Est)
% generate permutation matrix using correlation between original one and
% estimated one
p=size(Orig,2);
corr = corrcoef([Orig Est]);
corr = abs(corr(p+1:2*p,1:p));
perm = zeros(p);
aux=zeros(p,1);
for i=1:p
    [ri, ci]=find(max(corr(:))==corr); 
    ri=ri(1); ci=ci(1); % in the case of more than one maximum
    perm(ri,ci)=1;
    corr(:,ci)=aux; corr(ri,:)=aux';
end