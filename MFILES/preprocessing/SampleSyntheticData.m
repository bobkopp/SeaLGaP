function y=SampleSyntheticData(hmaps,hmaptimes,sheetesl,lat,long,dattemplate,maxvector)
	
% SampleSyntheticData(hmaps,hmaptimes,sheetesl,lat,long,dattemplate,
%					  [maxAgesd maxSLsd [maxUpliftRateSD] ]);
%
% Given real sea level values, sample at the location and with the errors given by
% dattemplate.
%
% EXAMPLE:
%
%    load hmaps-jxm
%    sldb=SealevelDBImport(fullfile(getenv('IFILES'),'sealevel.csv'));
%    sub=find(hmapids==5);
%    synsl=SampleSyntheticData(hmaps(:,sub),hmaptimes(sub),sheetesl(:,sub),lat,long,sldb);
%    synsl_lowerror = SampleSyntheticData(hmaps(:,sub),hmaptimes(sub),sheetesl(:,sub),...
%                     lat,long,sldb,[.1 .1 1e-6]);
%    save('syndata','synsl','synsl_lowerror');
% 
% Last updated by Bob Kopp rkopp-at-princeton.edu, 17 May 2009
	
	rand('state',sum(100*clock));
	randn('state',sum(100*clock));
	
	defval('maxvector',[Inf Inf Inf]);
	defval('maxUpliftRateSD',Inf);
	
	maxAgesd = maxvector(1); maxSLsd = maxvector(2);
	if length(maxvector)>2
		maxUpliftRateSD = maxvector(3);
	end
	
	[sorter,sortedi] = sort(length(dattemplate.siteid)*dattemplate.siteid + dattemplate.stratseqid); 

	Ndata = length(dattemplate.siteid);
	fields = fieldnames(dattemplate);
	for i=1:length(fields)
		if length(dattemplate.(fields{i})) == Ndata
			dattemplate.(fields{i}) = dattemplate.(fields{i})(sortedi);
		end
	end
	
	limitingDepth = 30;
	
	y=dattemplate;
	

	
	predecessor = zeros*ones(size(y.age_mode));
	successor = zeros*ones(size(y.age_mode));
	for i=1:Ndata
		pred = find((y.stratseqid == (y.stratseqid(i)-1)) .* (y.siteid == y.siteid(i)));
		succ = find((y.stratseqid == (y.stratseqid(i)+1)) .* (y.siteid == y.siteid(i)));
		if length(pred) > 0
			predecessor(i) = pred(1);
		end
		if length(succ) > 0
			successor(i) = succ(1);
		end
	end
	sequenceControlled = find((predecessor>0).*(isfinite(y.strataftertime_mode)).*(isfinite(y.strataftertime_sd)).*(y.strataftertime_sd<y.age_sd));
	
	for k=1:Ndata
		succset{k} = k;
		if successor(k)>0
			m=k;
			while ((successor(m))*ismember(successor(m),sequenceControlled)) > 0
				m=successor(m);
				succset{k}=[succset{k} m];
			end
		end
	end

	sub = find(isfinite(y.strataftertime_sd));
	y.strataftertime_sd(sub)=min(maxAgesd,y.strataftertime_sd(sub));

	y.age_sd = min(maxAgesd,y.age_sd(:));
	perturbablet_sd = y.age_sd;
	perturbablet_sd(sequenceControlled) = y.strataftertime_sd(sequenceControlled);

	y.age_real = y.age_mode;
	
	startrandomizer = randn(size(y.age_real)) .* perturbablet_sd;
	for k=1:length(y.age_real)
		minoffset = -Inf;
		maxoffset = Inf;
		offset = startrandomizer(k);
		if predecessor(k) > 0
			maxoffset = y.age_real(predecessor(k)) - y.age_real(k);
		end
		if successor(k) > 0
			minoffset = y.age_real(successor(k)) - y.age_real(k);
		end
		offset = max(minoffset,offset); offset = min(maxoffset,offset);
		y.age_real(succset{k}) = y.age_real(succset{k}) + offset;
	end
	
	y.upliftrate_sd = min(maxUpliftRateSD,y.upliftrate_sd(:));
	y.upliftrate_real = y.upliftrate_mode;
	
	imprecise = find(y.upliftrate_sd > 0);
	
	y.upliftrate_real(imprecise) = y.upliftrate_mode(imprecise) + randn(size(y.upliftrate_mode(imprecise))).*y.upliftrate_sd(imprecise);
	y.uplift_real = y.upliftrate_real .* y.age_real;
	
	y.uplift_mode = y.upliftrate_mode .* y.age_mode;
	y.uplift_sd = abs ( (y.age_sd ./ y.age_mode) .* y.uplift_mode );
	y.uplift_sd(imprecise) = abs ( sqrt ( (y.upliftrate_sd(imprecise)./(1e-9+y.upliftrate_mode(imprecise))).^2 + (y.age_sd(imprecise)./y.age_mode(imprecise)).^2) .* y.uplift_mode(imprecise) );
	
	subset1 = find((y.age_real<=max(hmaptimes(:))) .* (y.age_real>=min(hmaptimes(:))));
%	subset = intersect(subset1,find(isfinite(y.lat).*isfinite(y.long)));
	

% here we get the data from the inputs


	sumesl(1,:) = sum(sheetesl([1 2 3 6],:));
	sumesl(2,:) = sum(sheetesl([4 5 7],:));
	sumesl(3,:) = -sum(sheetesl);
	
	hmaps = [hmaps ; sheetesl ; sumesl];
	lat = [lat(:) ; repmat(-1e6,[size(sheetesl,1) 1]) ; repmat(1e6,[size(sumesl,1) 1])];
	long = [long(:) ; [1:size(sheetesl,1)]' ; [1 2 1e6]' ];

	subsl=find(isfinite(lat));
	subices = find(lat==-Inf); subsumices = find(lat==Inf);
	
	nearest = ones([1 length(y.lat)]);
	y.lat(find(y.lat==Inf))  =  1e6;
	y.lat(find(y.lat==-Inf)) = -1e6;
	y.long(find(y.long==Inf))  =  1e6;

	% identify nearest geographic location
	diff1x = repmat(y.long,[1 length(long)]) - repmat(long',[length(y.long) 1]);
	diff1y = repmat(y.lat,[1 length(lat)]) - repmat(lat', [length(y.lat) 1]);
	dist1 = sqrt(diff1x.^2 + diff1y.^2);
	[mindist1,mindist1i] = min(dist1,[],2);
	nearest = mindist1i;

	% okay, so now the coordinates of the geographic locations are in nearest
	% we need to interpolate to account for time
	
	y.sl_real = NaN(size(y.sl_mode));
	y.sl_real(subset1,1) = diag(interp1(hmaptimes,hmaps(nearest(subset1),:)',y.age_real(subset1)));

%%%
	y.rsl_real = y.sl_real + y.uplift_real;
	
%	y.sl_sd = min([repmat(maxSLsd,[length(y.sl_sd) 1]) y.sl_sd(:)],[],2);
	y.rsl_sd = min(maxSLsd,y.rsl_sd(:));
	y.rsl_sd = min(y.rsl_sd, y.sl_sd);
	y.sl_sd = sqrt( y.rsl_sd.^2 + y.uplift_sd.^2 );
	
	y.rsl_mode = dattemplate.rsl_mode * NaN;
	y.sl_mode = dattemplate.sl_mode * NaN;

	y.rsl_mode = y.rsl_real + randn(size(y.rsl_real)) .* y.rsl_sd;
	
	s1 = find((dattemplate.sl_limiting(:) == -1) .* isfinite(y.sl_real(:)))';
	y.rsl_mode(s1) = y.rsl_mode(s1) - abs(randn(size(y.rsl_mode(s1)))) * limitingDepth;		

	s1 = find((dattemplate.sl_limiting(:) == 1) .* isfinite(y.sl_real(:)))';
	y.rsl_mode(s1) = y.rsl_mode(s1) + abs(randn(size(y.rsl_mode(s1)))) * limitingDepth;		
	
	y.sl_mode = y.rsl_mode - y.uplift_mode;
	
	subset = find(isfinite(y.sl_mode));
	N = length(y.sl_mode);
	fields = fieldnames(y);
	for i=1:length(fields)
		if length(y.(fields{i})) == N
			y.(fields{i}) = y.(fields{i})(subset);
		end
	end
	
	y.lat(find(y.lat>1e3))=Inf;
	y.lat(find(y.lat<-1e3))=-Inf;
	y.long(find(y.long>1e3))=Inf;
	y.long(find(y.long<-1e3))=-Inf;
	
end