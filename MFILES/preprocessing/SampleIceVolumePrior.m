function iced=SampleIceVolumePrior(ivtot,mu,cv,bins,t,timescale)

% iced=SampleIceVolumePrior(ivtot,mu,cv,bins,t,timescale)
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 12 June 2009

defval('t',1:length(ivtot));
defval('timescale',1);

for i=1:size(ivtot,2)
	sub=find((bins(1:end-1)<=(ivtot(i))).*(bins(2:end)>(ivtot(i))));
	if length(sub)==0
		if (ivtot(i))<bins(1)
			sub=1;
		elseif (ivtot(i))>bins(end)
			sub=(length(bins)-1);
		end
	end
	binned(i)=sub(1);
end

iced(1,:) = mvnrnd(mu(binned(1),:),cv(:,:,binned(1)),1);
i=2;
while i<=length(t)
	deltagiv=ivtot(i)-ivtot(i-1);
	basevelocity=sum(.05*diag(cv(:,:,binned(i))))/timescale + abs(deltagiv)/abs(t(i)-t(i-1));
	
	current=iced(i-1,:);
	acceptedtot=0; rejectedtot=0;
	accepted=0;
	
	while (~accepted)
		if rejectedtot>1000
			i=i-1;
			current=iced(i-1,:);
			deltagiv=ivtot(i)-ivtot(i-1);
			basevelocity=sum(.05*diag(cv(:,:,binned(i))))/timescale + abs(deltagiv)/abs(t(i)-t(i-1));
			rejectedtot=0;
		end			
	
		candidate=mvnrnd(mu(binned(i),:),cv(:,:,binned(i)),1);
		delta = candidate-current;
		deltasum=sum(abs(delta));
		ratTot = (1-2*(-.5+normcdf(deltasum/abs(t(i)-t(i-1)),0,basevelocity))) * (prod(1*(candidate>=0)));
		trialprob=rand;
		
		if trialprob<ratTot
			accepted=1;
			acceptedtot=acceptedtot+1;
%			disp(['   Accepted: ' num2str([ratA ratB ratTot trialprob])]);
		else
			rejectedtot=rejectedtot+1;
%			disp(['   Rejected: ' num2str([ratA ratB ratTot trialprob])]);
		end
	end
%	disp([num2str(acceptedtot) '/' num2str(acceptedtot+rejectedtot) ' accepted']);
	
	iced(i,:)=candidate;
	i=i+1;
end

iced=iced';

sub=find(iced<0); iced(sub)=0;

adj=ivtot./sum(iced);
iced=iced.*repmat(adj,[7 1]);