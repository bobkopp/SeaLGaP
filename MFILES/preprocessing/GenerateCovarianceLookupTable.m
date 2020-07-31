function [mu,sd,coef,cvscor,muoverall,sdoverall,coefoverall,cvscoroverall,LAT,LONG,TIMS] = GenerateCovarianceLookupTable(hmaps,hmaptimes,hmapids,icedist,steric,lat,long,outputfile)

% [mu,sd,coef,cvscor,muoverall,sdoverall,coefoverall,cvscoroverall,LAT,LONG,tims] = 
% GenerateCovarianceLookupTable(hmaps,hmaptimes,hmapids,icedist,steric,lat,long,outputfile)
%
% INPUTS
%	hmaps:			array of sea level slices, with different positions on 1st dim and realizations/tims
%					on second
%	hmaptimes:		ages corresponding to each column of hmaps
%	hmapids:		realization id for each column of hmaps
%	icedist:		ice sheet equivalent sea levels for each column of hmaps
%	steric:			mean steric sea level change for each column of hmaps
%	lat
%	long
%	outputfile:		file in which to save output (default: 'cvtable');
%
% OUTPUTS
%	mu:				means of each space-time point
%	sd:				std of each space-time point
%	coef:			truncated principal components for space-time points
%	cvscor:			covariance of principal components
%	muoverall:		means of each space point
%	sdoverall:		std of each space point
%	coefoverall:	truncated principal components for space points
%	cvscoroverall:	covariance of principal components
%	LAT:			latitude of each row in coef
%	LONG:			longitude of each row in coef
%	tims:			tims of each row in coef
%
% EXAMPLE
%
%   load hmaps
%   load seeds
%   [gsl,icedist,steric,visco]=MatchHMAPSeeds(hmapids,hmaptimes,seeds,150:-1:100);
%   [mu,sd,coef,cvscor,muoverall,sdoverall,coefoverall,cvscoroverall,LAT,LONG,tims] = ... 
%       GenerateCovarianceLookupTable(hmaps,hmaptimes,hmapids,icedist,steric,lat,long);
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 24 July 2009
%

defval('outputfile','cvtable');
defval('latentcutoff',1-1e-2);
defval('icevol0',70.5);

uids = unique(hmapids);
tims = hmaptimes(find(hmapids==uids(1)));

% add ESL sums

sumesl(1,:) = sum(icedist([1 2 3 6],:)); % NH ice volume
sumesl(2,:) = sum(icedist([4 5 7],:)); % SH ice volume
sumesl(3,:) = sum(icedist); % total ice volume
sumesl(4,:) = steric;
sumesl(5,:) = icevol0-sumesl(3,:) + steric;

% add steric sea level to hmaps
hmaps = bsxfun(@plus,hmaps,steric);

% append longitudes and latitudes for ice volumes

lat = [lat(:) ; repmat(-Inf,[size(icedist,1) 1]) ; repmat(Inf,[size(sumesl,1) 1])];
long = [long(:) ; [1:size(icedist,1)]' ; [1 2 3 11 Inf]' ];

% reshape hmaps to distinguish each space-time point
hm = reshape([hmaps ; icedist ; sumesl],(size(hmaps,1)+size(icedist,1)+size(sumesl,1))*length(tims),[]);

LAT = repmat(lat,[length(tims) 1]);
LONG = repmat(long,[length(tims) 1]);
TIMS = repmat(tims(:)',[length(lat) 1])'; TIMS=TIMS(:);

mu = mean(hm,2);
sd = std(hm,[],2);

hmnorm = bsxfun(@minus,hm,mu);
hmnorm = bsxfun(@rdivide,hmnorm,sd+eps);

disp('Computing principal components...'); tic
[coeff,score,latent]=princomp(hmnorm','econ');
toc

cutoff=sum((cumsum(latent)/sum(latent))<latentcutoff);
disp(['   Reduced from ' num2str(size(hmnorm,1)) ' dimensions to ' num2str(cutoff) ' dimensions.']);

coef=coeff(:,1:cutoff); scor=score(:,1:cutoff);
cvscor=cov(scor);

% fall back prior for outside of range

muoverall = mean([hmaps ; icedist ; sumesl],2);
sdoverall =  std([hmaps ; icedist ; sumesl],[],2);

hmoverallnorm = bsxfun(@minus,[hmaps ; icedist ; sumesl],muoverall);
hmoverallnorm = bsxfun(@rdivide,hmoverallnorm,sdoverall+eps);

disp('Computing principal components geographically...'); tic
[coeffoverall,scoreoverall,latentoverall]=princomp(hmoverallnorm','econ');
toc

cutoffoverall=sum((cumsum(latentoverall)/sum(latentoverall))<latentcutoff);
disp(['   Reduced from ' num2str(size(hmoverallnorm,1)) ' dimensions to ' num2str(cutoffoverall) ' dimensions.']);

coefoverall=coeffoverall(:,1:cutoffoverall); scoroverall=scoreoverall(:,1:cutoffoverall);
cvscoroverall=cov(scoroverall);

save(outputfile,'mu','sd','coef','cvscor','muoverall','sdoverall','coefoverall','cvscoroverall','LAT','LONG','TIMS','lat','long','tims');
