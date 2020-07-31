function [fIces,sdIces,VIces,tgood,maxmeanGSL,maxmeantimes,minmeanNH,minmeanSH,Vmean,maxGSL,maxtimes,minNH,minSH,minNHoverall,minSHoverall,minNHoveralltime,minSHoveralltime,randpaths]=IceRatios(data,fKeep,f_sdKeep,gKeep,f,sd,fIc,sdIc,VIc)

% [fIces,sdIces,VIces,tgood,maxGSL,maxtimes,minNH,minSH]=IceRatios(data,fKeep,f_sdKeep,gKeep,f,sd)
%
% Generate NH ice, SH ice, and GSL projections with full covariance matrix
% among them.
%
% This routine isn't quite a modular function -- more like a hybrid
% between a script and a function. You'll want to edit it, most likely.
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 3 June 2009

defval('Nsamples',100);

defval('DropGSLPoints',0);
defval('fIc',[]);
defval('sdIc',[]);
defval('VIc',[]);

defval('t',[140:-.5:110]');

% threshold gsl

smoothTemporalCV = [1 2.58 .1];

if ~exist('thresholdgsl','var')
	% prior variances
	
	icelat = [Inf Inf Inf Inf]';
	icelong = [1:3 Inf]';
	priorV = LookupCovariance([icelat icelong zeros(size(icelat))],[],'cvtable-jxm',smoothTemporalCV);
	priorsd = sqrt(diag(priorV));
	
	%testV = LookupCovariance([Inf*t Inf*t t],[],'cvtable-jxm',smoothTemporalCV);
	%prioravgsd = mean(sqrt(diag(avgop*testV*avgop')));
	%thresholdavggsl = prioravgsd * .15;
	thresholdgsl = priorsd(end) * .1;
end

% identifying time points with any sd below threshold

goodpts = (sd<thresholdgsl);
sgp = sum(goodpts,2);
goodtimes = sgp>0;
tgood=t(find(goodtimes));

% grid

icelat = [repmat(Inf,1,3)]';
icelong = [Inf 1 2 ]';
icelatR = repmat(icelat,[1 length(tgood)])';
icelongR = repmat(icelong,[1 length(tgood)])';
icetR = repmat(tgood(:)',[length(icelat) 1])';

icelatR=icelatR(:);
icelongR=icelongR(:);
icetR = icetR(:);


% invert to get vector of [NH SH GSL]' and associated covariance matrix

promptfreq = 20;
	
N=size(fKeep,2);
Nt = length(tgood);
Nices = length(icelat);

if DropGSLPoints
	sub=find(~isinf(data.long));
else
	sub=1:length(data.long);
end

if length(VIc)==0
	fIces = zeros(Nt*length(icelat),N);
	sdIces = zeros(Nt*length(icelat),N);
	VIces = zeros(Nt*length(icelat),Nt*length(icelat),N);
	
	tic
	for i=1:N
		if mod(i,promptfreq)==0
			disp(['    ' num2str(i) '/' num2str(size(fKeep,2)) '  [' num2str(toc) ' seconds]' ]);
		end

		[fIces(:,i),VIces(:,:,i)] = GPSpaceTimeRegression(fKeep(sub,i),data.lat(sub),data.long(sub),gKeep(sub,i),f_sdKeep(sub,i),icelatR(:),icelongR(:),icetR(:),[],smoothTemporalCV);

		sdIces(:,i)=sqrt(diag(VIces(:,:,i)));
		VIces(:,:,i) = tril(VIces(:,:,i)) + tril(VIces(:,:,i),-1)';
	end
else
	fIces = fIc;
	sdIces = sdIc;
	VIces = VIc;
end

if nargout > 4

	subGSL = find(icelongR==Inf);
	subNH = find(icelongR==1);
	subSH = find(icelongR==2);

	[maxmeanGSL,mi] = max(fIces(subGSL,:),[],1);
	maxmeantimes = tgood(mi);
	zeropos = ((1:size(fIces,2))-1)*size(fIces,1);
	
	miGSL = mi+zeropos;
	miNH = miGSL + length(subGSL);
	miSH = miNH + length(subNH);
	
	minmeanNH = fIces(miNH);
	minmeanSH = fIces(miSH);
	
	for i=1:size(fIces,2)
		
		Vmean(:,:,i) = VIces(mi(i) + [0 length(subGSL) (length(subGSL)+length(subNH))],mi(i) + [0 length(subGSL) (length(subGSL)+length(subNH))],i); 

	end
	
	mmean = [maxmeanGSL(:)' ; minmeanNH(:)' ; minmeanSH(:)'];
end

if nargout>9
	disp('Sampling...')

	randpaths=NaN * ones([size(fIces,1) Nsamples size(fIces,2)]);
	
	tic

% generate random paths from each sample

	for i=1:size(fIces,2)
		if mod(i,promptfreq)==0
			disp(['    ' num2str(i) '/' num2str(size(fKeep,2)) '  [' num2str(toc) ' seconds]' ]);
		end
		
		try
			randpaths(:,:,i) = mvnrnd(fIces(:,i)',VIces(:,:,i),Nsamples)';
		catch
			[m,n] = size(VIces(:,:,i));
			[U2,S2,V2] = svd(VIces(:,:,i),0);
	   		s2 = diag(S2); 		
	   		if sum(log(s2))<log(eps)
%					disp('    Determinant of covariance = 0!');
%					disp('    Searching for coefficient of correction...');
				guess = sqrt(exp(log(eps)/(1+m)));
				lambda = fzero( @(x) sum(log(s2+(x^2)))-log(1e-10),guess)^2;
				lambda = max(lambda);
%					disp(['    Adding ' num2str(lambda) '*I']); 
			end
			randpaths(:,:,i) = mvnrnd(fIces(:,i)',VIces(:,:,i)+eye([m n])*lambda,Nsamples)';
		end
	end

% identify time point of maximum GSL (minimum ice volume)
% and record NH, SH, and GSL at this time point

	subGSL = find(icelongR==Inf);
	subNH = find(icelongR==1);
	subSH = find(icelongR==2);
	sublogNH = find(icelongR==11);
	sublogSH = find(icelongR==12);
	
	randpaths = reshape(randpaths,Nt*length(icelong),[]);
	
	[maxGSL,mi] = max(randpaths(subGSL,:),[],1);
	maxtimes = tgood(mi);
	zeropos = ((1:size(randpaths,2))-1)*size(randpaths,1);
	
	miGSL = mi + zeropos;
	miNH = miGSL + length(subGSL);
	miSH = miNH + length(subNH);
	
	minNH = randpaths(miNH);
	minSH = randpaths(miSH);
	
	[minNHoverall,mi]=min(randpaths(subNH,:));
	minNHoveralltime = tgood(mi);
	[minSHoverall,mi]=min(randpaths(subSH,:));
	minSHoveralltime = tgood(mi);
	end


