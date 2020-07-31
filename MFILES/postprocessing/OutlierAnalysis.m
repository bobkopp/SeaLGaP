function [pLSLmean,pAgemean,pLSL,pAge,chisqLSL,chisqAge]=OutlierAnalysis(map,fKeep,f_sdKeep,gKeep,datids)

% [pLSLmean,pAgeman,pLSL,pAge,chisqLSL,chisqAge]=OutlierAnalysis(map,fKeep,f_sdKeep,gKeep,datids)
%
% Outlier analysis.
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 6 June 2009

defval('datids',map.dataids);

ntrials = size(fKeep,2);
[c,ia,ib]=intersect(map.dataids,datids);

Norig=length(map.sl_mode);
map=MapSubset(map,ia);
datids=datids(ib);
fKeep=fKeep(ib,:);
f_sdKeep=f_sdKeep(ib,:);
gKeep=gKeep(ib,:);


predecessor = zeros*ones(size(map.age_mode));
successor = zeros*ones(size(map.age_mode));
for i=1:length(ia)
	pred = find((map.stratseqid == (map.stratseqid(i)-1)) .* (map.siteid == map.siteid(i)));
	succ = find((map.stratseqid == (map.stratseqid(i)+1)) .* (map.siteid == map.siteid(i)));
	if length(pred) > 0
		predecessor(i) = pred(1);
	end
	if length(succ) > 0
		successor(i) = succ(1);
	end
end
sequenceControlled = find((predecessor>0).*(isfinite(map.strataftertime_mode)).*(isfinite(map.strataftertime_sd)).*(map.strataftertime_sd<map.age_sd));

gKeep(sequenceControlled,:) = -gKeep(sequenceControlled,:) + gKeep(predecessor(sequenceControlled),:);
map.age_mode(sequenceControlled) = map.strataftertime_mode(sequenceControlled);
map.age_sd(sequenceControlled) = map.strataftertime_sd(sequenceControlled);

matsl_mode = repmat(map.sl_mode,[1 ntrials]);
matsl_sd = repmat(map.sl_sd,[1 ntrials]);

matage_mode = repmat(map.age_mode,[1 ntrials]);
matage_sd = repmat(map.age_sd,[1 ntrials]);

clear chisqLSL chisqAge;
chisqLSL=ones(Norig,ntrials)*NaN;
chisqAge=ones(Norig,ntrials)*NaN;

chisqLSL(ia,:) = (matsl_mode - fKeep).^2 ./ ( (matsl_sd.^2 + f_sdKeep.^2) ); 
chisqAge(ia,:) = (matage_mode - gKeep).^2 ./ ( (matage_sd.^2) ); 

pLSL = 1-chi2cdf(chisqLSL,1);
pAge = 1-chi2cdf(chisqAge,1);

sub = find(map.sl_limiting == 1);
pLSL(ia(sub),:) = normcdf((matsl_mode(sub,:) - fKeep(sub,:)) ./ sqrt(matsl_sd(sub,:).^2 + f_sdKeep(sub,:).^2));
sub = find(map.sl_limiting == -1);
pLSL(ia(sub),:) = normcdf(-(matsl_mode(sub,:) - fKeep(sub,:)) ./ sqrt(matsl_sd(sub,:).^2 + f_sdKeep(sub,:).^2));

pLSLmean=mean(pLSL,2);
pAgemean=mean(pAge,2);
