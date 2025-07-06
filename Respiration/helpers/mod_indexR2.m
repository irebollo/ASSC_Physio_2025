function mi = mod_indexR2(x,y);
% modulation index Tort
% now version for respiration


nbin=20;
%nbin=round(exp(0.626+0.4*log(length(x)-1)));
lnbin=log(nbin);
%work on whole array
edg=eqpop(x,nbin);
[~,bin] = histc(x,edg);
%[~,bin] = histc(x, linspace(-pi,pi,nbin+1));  % binned low frequency phase
binamp = zeros (size(y,2),nbin);      % binned amplitude

% for ip=1:nbin,
%     nip(ip)=length(find(bin==ip));
% end
nip=length(find(bin==1)); %sample per bin

for k = 1:nbin
    binamp(:,k) = mean(y(bin==k,:),1);
    %binamp(:,k) = trimmean(y(bin==k,:),20,'round',1);
end
binampn=binamp./repmat(sum(binamp,2),[1 nbin]);
MI=(lnbin-(-sum(binampn.*log(binampn),2)))/lnbin;

binl=5;
nseg=floor(nip/binl);
rind=repmat([1:binl],nseg,1);
nl=length(x);
MI2=zeros(200,60); % enter number of frequencies as second argument here!
for ib=1:200,
    for k = 1:nbin
        %create 5sample long continuous random vector
        rrind=rind'+repmat(randi([1 length(x)-10],nseg,1),1,binl)';     
        %rrind=rind+repmat(randperm(length(x)-5,nseg),binl,1)';
        rrind=mod(rrind(:),nl);
        rrind(rrind==0)=[];
        binamp(:,k) = mean(y(rrind(:),:),1);

        %binamp(:,k) = mean(y(randi(length(x),1,nip),:),1);
        %y2=circshift(y,randi([1000 nl]),1);
        %binamp(:,k) = mean(y2(bin==k,:),1);
        %binamp(:,k) = trimmean(y(randi(length(x),1,nip),:),20,'round',1);

    end
    binampn=binamp./repmat(sum(binamp,2),[1 nbin]);
    MI2(ib,:)=(lnbin-(-sum(binampn.*log(binampn),2)))/lnbin;
end
bm=mean(MI2,1);bs=std(MI2,0,1);
mi=(MI'-bm)./bs;

