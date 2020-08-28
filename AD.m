
function [ad, rmsAD]=AD(m,mh)
% function [sad, rmsSAD]=SAD(estEM, refEM)
% compute spectral or abundance angle distance between
% original and estimated signatures and RMS value of sad
% p=size(m,2);
% sad=zeros(1,p);
% for i=1:p
%     mt=m(:,i);
%     mht=mh(:,i);
%     sad(i)=rad2deg(acos((mt'*mht)/(norm(mt)*norm(mht))));
% end

% RMS value
% rmsSAD=sqrt((sum(sad.^2)/p));

nm = diag(m*m'); 
nmh = diag(mh*mh');
dp = diag(m*mh')./sqrt(nm.*nmh);
% dp(dp>1) = 1;
ad = acos(dp); % 180/pi*
ad(isnan(ad)) = 0;
rmsAD = mean(ad.^2)^.5;