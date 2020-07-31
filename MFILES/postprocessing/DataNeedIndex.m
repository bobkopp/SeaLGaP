function [dni,testlong,testlat]=DataNeedIndex(map,fname,gsl,gsl_sd,gsl_t,samplefreq)

% [dni,testlong,testlat]=DataNeedIndex(map,fname)
%
% Computes the Data Need Index.
%
% INPUTS:
%	map		sea level database structure
%	fname	file name of stored runs
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 15 August 2009

defval('timeslices',114:129);
defval('testlat',-90:10:90);
defval('testlong',0:15:360);
defval('samplefreq',20);
defval('dividinglongs',[-10 180 360]);
defval('savebackups',1);
defval('gsl',[]);
defval('sealevelrange',[-10 Inf]);

Ntimeslices = length(timeslices);

maplat=repmat(testlat(:),[1,length(testlong),Ntimeslices]);
maplong=repmat(testlong(:)',[length(testlat),1,Ntimeslices]);
maptimes=repmat(reshape(timeslices,1,1,[]),[length(testlat),length(testlong),1]);

load(fname,'fKeep','f_sdKeep','gKeep','datids');

[c,ia,ib]=intersect(map.dataids,datids);

disp('Calculating prior variances...')
priorvar = [];
for i=1:Ntimeslices
		for k=2:length(dividinglongs)
			sub=find(maptimes(:)==timeslices(i));
			sub = intersect(sub,find((maplong(:)>dividinglongs(k-1)).*(maplong(:)<=dividinglongs(k))));
			priorvar = [priorvar ; diag(LookupCovariance([maplat(sub) maplong(sub) maptimes(sub)]))];
	end
end

disp('Calculating posterior variances...')
meanfMap = []; varMap=[];
tic
col = 1;
for m=samplefreq:samplefreq:size(gKeep,2)
	meanfMap2 = []; varMap2=[];
	for i=1:Ntimeslices
		for k=2:length(dividinglongs)
			sub=find(maptimes(:)==timeslices(i));
			sub = intersect(sub,find((maplong(:)>dividinglongs(k-1)).*(maplong(:)<=dividinglongs(k))));
			[meanfMap1,V1,logpPrior1] = GPSpaceTimeRegression(fKeep(ib,m),map.lat(ia),map.long(ia),gKeep(ib,m),f_sdKeep(ib,m),maplat(sub),maplong(sub),maptimes(sub));
			
			varMap1 = diag(V1);
		
			meanfMap2 = [meanfMap2 ; meanfMap1];
			varMap2 = [varMap2 ; varMap1];
		end
	end

	disp(['   ' num2str(m) '/' num2str(size(gKeep,2)) ' [' num2str(toc) ' seconds]']);

	meanfMap(:,col) = meanfMap2;
	varMap(:,col) = varMap2;
	col=col+1;
end
toc

fMaps = meanfMap; varMaps = varMap;
if savebackups
	save timeslicemaps maplat maplong maptimes fMaps varMaps
end

varratio = bsxfun(@rdivide,varMaps,priorvar);
if length(gsl)==0
	goodslices = ones(size(fMaps));
else
	sub=samplefreq:samplefreq:size(gsl,2);
	gsl2 = interp1(gsl_t,gsl(:,sub),timeslices,'nearest','extrap');
	gsl2 = shiftdim(repmat(gsl2,[1 1 (length(testlat) * length(testlong))]),2);
	gsl2 = reshape(gsl2,size(fMaps));

	gsl2_sd = interp1(gsl_t,gsl_sd(:,sub),timeslices,'nearest','extrap');
	gsl2_sd = shiftdim(repmat(gsl2_sd,[1 1 (length(testlat) * length(testlong))]),2);
	gsl2_sd = reshape(gsl2_sd,size(fMaps));
	
	goodpts = normcdf(sealevelrange(2),gsl2,gsl2_sd) - normcdf(sealevelrange(1),gsl2,gsl2_sd);
end
meanvarratio = sum(goodpts.*varratio,2)./sum(goodpts,2);
dni = mean(reshape(meanvarratio,length(testlat),length(testlong),Ntimeslices),3);


if savebackups
	save dni testlong testlat dni 
end