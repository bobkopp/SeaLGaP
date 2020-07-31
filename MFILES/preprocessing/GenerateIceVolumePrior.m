function [cv,mu,sd,bins,eslx,ivx]=GenerateIceVolumePrior(models)

% [cv,mu,sd,bins,eslx,ivx]=GenerateIceVolumePrior(models)
%
% Generates the prior for ice volume conditioned on
% GSL.
%
% INPUTS:
%	models	cell array of equally-weighted ice models
%
% OUTPUTS:
%	cv	covariances for each bin
%	mu	means for each bin
%	bins	bounds of each bin
%	eslx	synthetic equivalent sea level ice volumes
%	ivx		synthetic total ice volume
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 10 June 2009

defval('Nsamps',200000);
%defval('icevolmaxloss',[.16 7.0 0.02 5.0 58.0 0 0]);

for i=1:length(models)
	clear eslnew ivnew;
	
	esl0(:,i) = models{i}.sheetesl(:,end);
	iv0(i)=sum(esl0(:,i));
	esl{i}=models{i}.sheetesl(:,end:-1:1);
	iv{i}=sum(esl{i});
	Nsampseff(i)=ceil(Nsamps/length(iv{i}));
	
	desl{i}=diff(esl{i},[],2);
	div{i} =diff(iv{i},[],2);

	rnd=3.^(2*rand([size(desl{i}) Nsampseff(i)])-1);
	rndgsl=1.5.^(2*rand([size(div{i}) Nsampseff(i)])-1);

	eslnew(:,1,:) = repmat(esl0(:,i),[1 1 Nsampseff(i)]);
	eslnew(:,2:size(desl{i},2)+1,:) = repmat(desl{i},[1 1 Nsampseff(i)]).*rnd;
	eslnew = cumsum(eslnew,2);
	
	ivnew(1,1,:) = repmat(iv0(i),[1 1 Nsampseff(i)]);
	ivnew(1,2:size(div{i},2)+1,:) = repmat(div{i},[1 1 Nsampseff(i)]).*rndgsl;
	ivnew=cumsum(ivnew,2);
	
	corrections=ivnew./sum(eslnew,1);
	eslnew = eslnew .* repmat(corrections,[size(eslnew,1) 1 1]);

	% then add deviations to highstands
	sub=find(ivnew<(iv0(i)+10));
	rshift1=2*(rand(1,length(sub))-.5);
	rshift2=2*(rand(1,length(sub))-.5);
	rshift3=2*(rand(1,length(sub))-.5);
	rshift4=2*(rand(1,length(sub))-.5);
	rshift4b=2*randn(1,length(sub));
	
	shiftWAIS = esl0(4,i)*rshift1;
	NHmin=sum(esl0(1:3,i));
	shiftNHfrac = esl0(4,i)*rshift2/NHmin + ((NHmin-esl0(4,i))/NHmin)*rshift3;
	shiftNA=esl0(1,i)*shiftNHfrac;
	shiftEur=esl0(3,i)*shiftNHfrac;
	shiftGIS=esl0(2,i)*shiftNHfrac;
	shiftEAIS=(NHmin-esl0(4,i))*rshift4+rshift4b;
	
	eslnew(:,sub)=eslnew(:,sub)-[shiftNA ; shiftGIS ; shiftEur ; shiftWAIS ; shiftEAIS ; 0*shiftNA ; 0*shiftNA];
	
	sub=find(eslnew<0);
	eslnew(sub)=0;

	% then add some scenarios where only EAIS is left
	rshift5=.4*esl0(5,i)*rand(1,100);
	sub=size(eslnew,2)+[1:length(rshift5)];
	eslnew([1:4 6:7],sub)=0;
	eslnew(5,sub)=esl0(5,i)-rshift5;

	esl2{i} = eslnew;
	iv2{i} = sum(eslnew);
	
end



eslx=[]; ivx=[];
for i=1:length(esl2)
	eslx=[eslx reshape(esl2{i},size(esl2{1},1),[],1)];
	ivx=[ivx reshape(iv2{i},1,[],1)];
end



% chop up

gslrnd = round(ivx);
u=unique(gslrnd);
clear cnt;
for i=1:length(u)
 cnt(i)=sum(round(gslrnd)==u(i));
end
a=interp1(cumsum(cnt),u,[100 1000:1000:length(gslrnd)]);
a=floor(a)+.5;
a=[-0.001 a];
if max(a)<max(ivx)
 a(end+1)=ceil(max(ivx))+.5;
end
a=unique(a);
bins = a;
clear cv mu sd;
for i=1:(length(bins)-1)
 sub=find((ivx>=bins(i)).*(ivx<bins(i+1)));
 mu(i,:)=mean(eslx(:,sub)');
 sd(i,:)=std(eslx(:,sub)');
 sd=max(sd,.001);
 denormed=(eslx(:,sub)'-repmat(mu(i,:),[length(sub) 1])) ./ repmat(sd(i,:),[length(sub) 1]);

 cv(:,:,i)=cov(eslx(:,sub)');
%	cv(:,:,i)=cov(denormed);
end
