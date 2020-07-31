function y=AssimilateRun(files,Nburnin)

% y=AssimilateRun(files,Nburnin)
%
% Assimilates output data from multiple Monte Carlo runs into a single structure.
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 17 May 2009

	defval('Nburnin',0);

	load(files(1).name);
	anchor{1} = intersect(find(max(diff(gKeep,1,2)')==0),find(data.age_mode==0));
	y.acceptanceRateKeep=[acceptanceRateKeep(:,(Nburnin+1):end)];
	y.fKeep=[fKeep(:,(Nburnin+1):end)];
	y.f_sdKeep=[f_sdKeep(:,(Nburnin+1):end)];

%	if length(anchor{1}) > 0
%		gKeep = gKeep - gKeep(anchor{1}(1),1);
%	end
	y.ageOffset(1) = ageOffset;
	
	y.gKeep=[gKeep(:,(Nburnin+1):end)] + ageOffset;
	if exist('f_VKeep')
		y.f_VKeep=[f_VKeep(:,:,(Nburnin+1):end)];
		clear f_VKeep;
	end
	
	y.runid = ones(1,size(y.fKeep,2));
	y.gelman.runmean(:,1) = mean(y.gKeep,2);
	y.gelman.N(1) = size(y.gKeep,2);
	y.gelman.withinrunvar(:,1) = sum((y.gKeep - repmat(y.gelman.runmean(:,1),[1 y.gelman.N(1)])).^2,2) / (y.gelman.N(1)-1);
	
	y.data.dataids = data.dataids;
	y.data.long = data.long;
	y.data.lat = data.lat;

	for i=2:length(files)
		load(files(i).name);
		anchor{i} = intersect(find(max(diff(gKeep,1,2)')==0),find(data.age_mode==0));
%		if length(anchor{i}) > 0
%			gKeep = gKeep - gKeep(anchor{i}(1),1);
%		end
		y.ageOffset(i) = ageOffset;

		paired = []; matcher = [];
		for j=1:length(y.data.dataids)
			u = find(data.dataids==y.data.dataids(j),1);
			if length(u) > 0
				paired(end+1) = j;
				matcher(end+1) = u;
			end
		end
		
		y.acceptanceRateKeep = y.acceptanceRateKeep(paired,:);
		y.fKeep = y.fKeep(paired,:);
		y.f_sdKeep = y.f_sdKeep(paired,:);
		y.gKeep = y.gKeep(paired,:);
		
		y.gelman.runmean = y.gelman.runmean(paired,:);
		y.gelman.withinrunvar = y.gelman.withinrunvar(paired,:);
		
		if (exist('y.f_VKeep'))
			if size(y.f_VKeep,1)>0
				if exist('f_VKeep')
					y.f_VKeep = y.f_VKeep(paired,paired,:);
					y.f_VKeep(:,:,end+1:(end+size(f_VKeep,3)-Nburnin)) = f_VKeep(matcher,matcher,Nburnin+1:end);
				else
					y.f_VKeep = [];
				end
			end
		end

		y.data.dataids = y.data.dataids(paired,:);
		y.data.long = y.data.long(paired,:);
		y.data.lat = y.data.lat(paired,:);
		
		y.acceptanceRateKeep=[y.acceptanceRateKeep acceptanceRateKeep(matcher,(Nburnin+1):end)];
		y.fKeep=[y.fKeep fKeep(matcher,(Nburnin+1):end)];
		y.f_sdKeep=[y.f_sdKeep f_sdKeep(matcher,(Nburnin+1):end)];
		y.gKeep=[y.gKeep gKeep(matcher,(Nburnin+1):end)+ ageOffset];
		y.runid = [y.runid i*ones(1,size(gKeep(matcher,(Nburnin+1):end),2))];
		y.gelman.runmean(:,i) = mean(gKeep(matcher,(Nburnin+1):end)+ageOffset,2);
		y.gelman.N(i) = size(gKeep(matcher,(Nburnin+1):end),2);
		y.gelman.withinrunvar(:,i) = sum((gKeep(matcher,(Nburnin+1):end)+ageOffset - repmat(y.gelman.runmean(:,1),[1 y.gelman.N(i)])).^2,2) / (y.gelman.N(i)-1);


	end
	y.gelman.overallmean = mean(y.gelman.runmean,2);
	y.gelman.W = mean(y.gelman.withinrunvar,2);
	vars = (y.gelman.runmean - repmat(y.gelman.overallmean,[1 length(y.gelman.N)])).^2;
	vars = vars .* repmat(y.gelman.N,[length(y.gelman.runmean) 1]);
	y.gelman.B = sum(vars/(length(y.gelman.N)-1),2);
	Nmin = min(y.gelman.N);
	y.gelman.varhat = ((Nmin-1)/Nmin) * y.gelman.W + (length(y.gelman.N)+1)/(length(y.gelman.N) * Nmin) * y.gelman.B; 
	y.gelman.varratios= mean(y.gelman.withinrunvar./repmat(y.gelman.varhat,[1 length(y.gelman.N)]));
	
	y.anchor = anchor;
	
end