function  [fKeep,gKeep,yKeep,f_VKeep,f_sdKeep,acceptanceRateKeep] = SealevelAnalyze(data,N,thin,varargin)

% [fKeep,gKeep,yKeep,f_VKeep,f_sdKeep,acceptanceRateKeep] = SealevelAnalyze(data,N,thin,[varargin])
%
% ALGORITHM:
%	 s is our sea level observations, including censored data
%	 y is our picked substitute for s
%	 f is true sea level
%	
%	 t are our observed ages
%	 g are true ages
%	
%	 Process priors:
%		f ~ GP
%		g ~ constant
%	
%	 Data likelihoods:
%		t | g ~ N(g,E_t)
%		y | f, g ~ N(f(g), E_s)
%		s | f, g ~ N(f(g), E_s) * I(y in I)
%	
%	 We're searching for the posterior p(f,g | s,t)
%	
%	 having already determined the covariance function, 
%	 our algorithm is as follows:
%	
%	 (1) Pick y given s: y ~ p(y | f, g, s) ~ N(f(g), E_s) * I(y in I)
%	 (2) Krieg for f ~ p(f | g, t, y)
%	 (3) MCMC for g ~ p (g | t, f, y) ~ p (t | g, f, y) * p(g | f, y)
%	                                  ~ p (t | g) * p(f | g) * p(g)
%	
%	 This version separates treatment of RSL and uplift. Substitute "s - u" for s.
%
%
%
% INPUT:
%	data			sea level data structure
%	N				number of Monte Carlo steps
%	thin			how many MC steps to run per steps stored
%
% VALID OPTIONS FOR VARARGIN:
%	'minage'		specifies minimum age of points considered in next parameter
%	'maxage'		specifies maixmum age of points considered in next parameter
%	'label'			the next parameter is the label for run (Default: SeaLevelAnalyze)
%	'store'			the next parameter is a structure with additional information to save with output
%	'timeres'		the next parameter is the minimum age error to require a covariance reevaluation
%	'seed'			the next parameter is a structure containing starting values for
%					yKeep, fKeep, gKeep, f_sdKeep, and f_VKeep
%	'mcstepsize'		the next parameter is the initial step size for age adjustments
%	'grandomizestart'	randomize starting ages
%	'showplot'		produce plots as run progresses
%	'noplot'		do not produce plots (default)
%   'cvfile'		covariance file to use
%   'dataweight'    weight factor to give posterior in MC (vs. prior)
%
% OUTPUT:
%	fKeep			stored true sea level values
%	gKeep			stored true age values
%	yKeep			stored observed sea level values
%	f_VKeep			stored covariance matrices of f (3D array)
%	f_sdKeep		stored standard deviations of f
%	acceptanceRateKeep	stored acceptance rate of perturbations to g
%
% Results are also saved to the file [label '.' hostname '.' processid], e.g. 'SeaLevelAnalyze.bobsmac.1.mat'
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 7 August 2009

% first we set all the parameters

scriptid = '';
randn('state',sum(100*clock))
rand('state',sum(100*clock))

cvfile =[];
MCStepSize0 = 2;

adaptiveStep0 = 2; adaptiveDecayScale = 3;
adaptiveMin = .1; adaptiveMax = 8;
adaptationPeriod = 4;
idealAcceptanceProb = 0.24;
timeres = 0.01;

gRandomizeStart = 1;
showplot=0; debugmode = 0;
dataweight = 1;
tprsc = [];

minage = -Inf; maxage = Inf;
	
store=struct();

if length(varargin)>0
	for i=1:length(varargin)
		if strcmpi(varargin{i},'minage')
			minage = varargin{i+1};
			i=i+1;
		elseif strcmpi(varargin{i},'maxage')
			maxage = varargin{i+1};
			i=i+1;
		elseif strcmpi(varargin{i},'label')
			scriptid = [scriptid varargin{i+1}];
			i=i+1;
		elseif strcmpi(varargin{i},'store')
			store = varargin{i+1};
			i=i+1;
		elseif strcmpi(varargin{i},'timeres')
			timeres = varargin{i+1};
			i=i+1;
		elseif strcmpi(varargin{i},'seed')
			seed = varargin{i+1};
			i = i+1;
		elseif strcmpi(varargin{i},'mcstepsize')
			MCStepSize0 = varargin{i+1};
			i = i+1;
		elseif strcmpi(varargin{i},'grandomizestart')
			gRandomizeStart = 1;
		elseif strcmpi(varargin{i},'showplot')
			showplot = 1;
		elseif strcmpi(varargin{i},'noplot')
			showplot = 0;
		elseif strcmpi(varargin{i},'cvfile')
			cvfile = varargin{i+1};
			i = i+1;
		elseif strcmpi(varargin{i},'dataweight')
			dataweight = varargin{i+1};
			i = i+1;
		elseif strcmpi(varargin{i},'tprsc')
			if length(varargin{i+1})>0
				tprsc = varargin{i+1};
			end
			i = i+1;
		end
	end
end

if length(scriptid) == 0
	scriptid = 'SealevelAnalyze';
end

processid = getenv('PBS_JOBID'); mpiexecid = 0;
PBSrun = 0;

sub=find((data.age_mode>=minage).*(data.age_mode<=maxage));
data=MapSubset(data,sub);
  
MCStepSize = repmat(MCStepSize0,size(data.sl_mode));
                                                                              
disp(['Max age: ' num2str(maxage)]);                                            
disp(['Min age: ' num2str(minage)]);                                            
disp(['Tprsc: ' num2str(tprsc)]);                                               
                                                                                
disp(['Oldest sample: ' num2str(max(data.age_mode))]);                          
                           

if length(processid) == 0
	
	[stat,hostname] = system('hostname');
	dots = strfind(hostname,'.');
	hostname=hostname(1:dots-1);
	
	
	processid = 1;
	while exist([scriptid '.' hostname '.' num2str(processid) '.mat'],'file')
		processid = processid+1;
	end
	uniqueid = [scriptid '.' hostname '.' num2str(processid)]
else
	PBSrun = 1;
	hostname = 'PBS';
	
	processid = getenv('PBS_JOBID');
	mpiexecid = getenv('MPIEXEC_RANK');

	uniqueid = [scriptid '.' hostname '.' num2str(processid)]

	if length(mpiexecid) > 0
		uniqueid = [uniqueid '.' num2str(mpiexecid)];
	else
		mpiexecid = 0;
	end
end

% sort the data by siteid

[sorter,sortedi] = sort(length(data.siteid)*data.siteid + data.stratseqid); 

Ndata = length(data.siteid);
fields = fieldnames(data);
for i=1:length(fields)
	if length(data.(fields{i})) == Ndata
		data.(fields{i}) = data.(fields{i})(sortedi);
	end
end
	
% if there's an anchor sequence, anchor it
ageOffset  = 0;
anchor=[];
sub=find(data.age_sd==0);
if length(sub) == 0
	[m,i] = max(data.stratseqid);
%	if m > 12
		disp(['Anchoring sequence ' num2str(data.siteid(i))]);
		anchor = find((data.stratseqid == 1).*(data.siteid == data.siteid(i)));
%	end
else
	anchor=sub(1);
end
if length(anchor) > 0
	data.age_sd(anchor) = 0;
	ageOffset = data.age_mode(anchor);
	data.age_mode = data.age_mode - ageOffset;
	if isfield(data,'age_real')
		data.age_real = data.age_real - data.age_real(anchor);
	end
end

% initialize variables

runcycles = floor(N/thin); N=runcycles * thin;

s = data.rsl_mode;
s_sd = data.rsl_sd;
t = data.age_mode;
t_sd = data.age_sd;

du = data.upliftrate_mode;
du_sd = data.upliftrate_sd;

perturbablet = t;
perturbablet_sd = t_sd;

% set up data with constrained relative ages

perturbablet_predecessor = zeros*ones(size(t));
perturbablet_successor = zeros*ones(size(t));
for i=1:length(t)
	pred = find((data.stratseqid == (data.stratseqid(i)-1)) .* (data.siteid == data.siteid(i)));
	succ = find((data.stratseqid == (data.stratseqid(i)+1)) .* (data.siteid == data.siteid(i)));
	if length(pred) > 0
		perturbablet_predecessor(i) = pred(1);
	end
	if length(succ) > 0
		perturbablet_successor(i) = succ(1);
	end
end
sequenceControlled = find((perturbablet_predecessor>0).*(isfinite(data.strataftertime_mode)).*(isfinite(data.strataftertime_sd)).*(data.strataftertime_sd<data.age_sd));

for k=1:length(t)
	perturbablet_succset{k} = k;
	if perturbablet_successor(k)>0
		m=k;
		while ((perturbablet_successor(m))*ismember(perturbablet_successor(m),sequenceControlled)) > 0
			m=perturbablet_successor(m);
			perturbablet_succset{k}=[perturbablet_succset{k} m];
		end
	end
end
perturbablet(sequenceControlled) = data.strataftertime_mode(sequenceControlled);
perturbablet_sd(sequenceControlled) = data.strataftertime_sd(sequenceControlled);

perturbablet_sd(anchor) = 0; % just in case we goofed

counter = 1;
i=1;


% if this is a synthetic data and we know the real values, compute some variables for comparision
if isfield(data,'sl_real') && isfield(data,'age_real')
	[fReal,VReal,logpReal] = GPSpaceTimeRegression(data.sl_real,data.lat,data.long,data.age_real,ones(size(data.sl_real))*.01,cvfile,tprsc);
	sub=find(data.sl_limiting==0);
	[fInitial,VInitial,logpInitial] = GPSpaceTimeRegression(data.sl_mode(sub),data.lat(sub),data.long(sub),data.age_real(sub),data.sl_sd(sub),cvfile,tprsc);
	PosteriorSL = sum ( ((data.sl_real(sub)-data.sl_mode(sub))./data.sl_sd(sub)).^2 )/length(sub);
	realt = data.age_real;
	realt(sequenceControlled) = data.age_real(perturbablet_predecessor(sequenceControlled)) - data.age_real(sequenceControlled);
	PosteriorAge =  sum( ((realt-perturbablet)./(1e-4+perturbablet_sd)).^2)/length(data.age_real);
end



% First we need to initialize

if ~exist('seed','var')

	% set initial ages to measured ages
	g(:,i) = t;
	
	% if we going to perturb the starting ages, do that
	if gRandomizeStart
		startrandomizer = randn(length(t),1) .* MCStepSize .* perturbablet_sd;
		for k=1:length(t)
			minoffset = -Inf;
			maxoffset = Inf;
			offset = startrandomizer(k);
			if perturbablet_predecessor(k) > 0
				maxoffset = g(perturbablet_predecessor(k),i) - g(k,i);
			end
			if perturbablet_successor(k) > 0
				minoffset = g(perturbablet_successor(k),i) - g(k,i);
			end
			offset = max(minoffset,offset); offset = min(maxoffset,offset);
			g(perturbablet_succset{k},i) = g(perturbablet_succset{k},i) + offset;
		end
	end
	
	% calculate uplift based on ages
	u = (ageOffset + g(:,i)) .* du;
	u_sd = (ageOffset + g(:,i)) .* du_sd;
	
	%disp(['   Anchor initialized at ' num2str(u(anchor))]);
	
	f_V = zeros([length(s) length(s) 1]);


	% separately handle indicative points (data.sl_limiting==0) and
	% upper or lower bounds
	
	sub1 = find(data.sl_limiting==0);
	sub2 = find(data.sl_limiting ~= 0);
	
	y(sub1,i) = s(sub1) - u(sub1); % measured sea level = observed RSL - uplift
	
	% GPR for f given g and y
	[f(sub1,i),f_V(sub1,sub1,i),logp,errorflags,logpgrad,cvfunc,evalfunc] = GPSpaceTimeRegression(y(sub1,i),data.lat(sub1),data.long(sub1),g(sub1,i),sqrt(s_sd(sub1).^2 + u_sd(sub1).^2),cvfile,tprsc);
	f_sd(sub1,i) = sqrt(diag(f_V(sub1,sub1,i)));
	
	% handle upper and lower bounds
	if length(sub2) > 0
		[f(sub2,i),f_V(sub2,sub2,i)] = evalfunc([data.lat(sub2),data.long(sub2),g(sub2,i)]);
		f_sd(sub2,i) = sqrt(diag(f_V(sub2,sub2,i)));
	end
	
	sub = find(data.sl_limiting == -1); % lower bounds to sea level
	y(sub,1) = uncensorData(f(sub),sqrt(f_sd(sub).^2 + s_sd(sub).^2 + u_sd(sub,1).^2),s(sub)-u(sub),Inf);
	sub=find(data.sl_limiting == 1); % upper bounds to sea level
	y(sub,1) = uncensorData(f(sub),sqrt(f_sd(sub).^2 + s_sd(sub).^2 + u_sd(sub,1).^2),-Inf,s(sub)-u(sub));
	
	yKeep(:,counter) = y(:,i);
	fKeep(:,counter) = f(:,i);
	gKeep(:,counter) = g(:,i);
	f_sdKeep(:,counter) = f_sd(:,i);
	f_VKeep(:,:,counter) = f_V(:,:,i);
else
	yKeep(:,counter) = seed.yKeep(:,end);
	fKeep(:,counter) = seed.fKeep(:,end);
	gKeep(:,counter) = seed.gKeep(:,end);
	f_sdKeep(:,counter) = seed.f_sdKeep(:,end);
	f_VKeep(:,:,counter) = seed.f_VKeep(:,:,end);
end

% Then we start our algorithm

totalcount = 1;
% loop through the desired number of run cycles
for counter=2:(runcycles+1)

	fid = fopen(['status.' uniqueid],'w');
	fprintf(fid,['%0g\t' 'of ' num2str(runcycles+1) '\t' datestr(now) '\t' uniqueid '\n'],counter);
	fclose(fid);

	clear y f g f_sd f_V acceptanceRate;
	g(:,1) = gKeep(:,counter-1);
	
	% then within each run cycle, load through "thin" iterations before saving
	for i=1:thin

		disp(['Starting run ' num2str((counter-2)*thin+i) ' [cycle ' num2str(counter-1) ']...'])
		
		% first pick y
	
		% uplift calculation
		u = (ageOffset + g(:,i)) .* du;
		u_sd = (ageOffset + g(:,i)) .* du_sd;
	
		sub1 = find(data.sl_limiting==0);
		sub2 = find(data.sl_limiting ~= 0);
		
		impute_sd = 0*g(:,1);
		y(sub1,i) = s(sub1) - u(sub1);
		
		% GPR for f at indicative points given y
		[ftemp(sub1),f_Vtemp(sub1,sub1),logp,errorflags,logpgrad,cvfunc,evalfunc] = GPSpaceTimeRegression(y(sub1,i),data.lat(sub1),data.long(sub1),g(sub1,i),sqrt(s_sd(sub1).^2 + u_sd(sub1).^2),cvfile,tprsc);
		if length(sub2)>0
			% if limiting points, use the indicative points to calculate f there
			[ftemp(sub2),f_Vtemp(sub2,sub2)] = evalfunc([data.lat(sub2),data.long(sub2),g(sub2,i)]);
			impute_sd(sub2) = sqrt(diag(f_Vtemp(sub2,sub2)));
		end

		% then impute y values given these fs at the limiting points
		ftemp = ftemp(:);
		sub=find(data.sl_limiting == -1);
		y(sub,i) = uncensorData(ftemp(sub),sqrt(impute_sd(sub).^2+s_sd(sub).^2 + u_sd(sub).^2),s(sub)-u(sub),Inf);
		
		sub=find(data.sl_limiting == 1);
		y(sub,i) = uncensorData(ftemp(sub),sqrt(impute_sd(sub).^2+s_sd(sub).^2 + u_sd(sub).^2),-Inf,s(sub)-u(sub));

%	 then krieg approximation for f	given the full y

		[f(:,i),f_V(:,:,i),logp,errorflags,logpgrad,cvfunc,evalfunc] = GPSpaceTimeRegression(y(:,i),data.lat,data.long,g(:,i),sqrt(s_sd.^2 + u_sd.^2 + impute_sd.^2),cvfile);
		f_sd(:,i) = sqrt(diag(f_V(:,:,i)));

% and finally do a sequential Metropolis step for g

		gold = g(:,i);
		fold = f(:,i); f_sdold = f_sd(:,i);
		
		[fold2,Vold2,logp,errorflags,logpgrad,cvfunc,evalfunc] = GPSpaceTimeRegression(f(:,i),data.lat,data.long,gold,f_sd(:,i),cvfile);
		modelLogLikelihoodOld = logp;

	
		% begin monitoring and statistics routines for cases where plot is requested
		% and/or synthetic data is being used
		haveRealData = (isfield(data,'sl_real') && isfield(data,'age_real'));
		if haveRealData || showplot
				tempt = g(:,i);
				tempt(sequenceControlled) = gold(perturbablet_predecessor(sequenceControlled)) - gold(sequenceControlled);
				mintime = min(gold); maxtime = max(gold); evaltimes = linspace(min(gold)-5,max(gold)+5,50); evaltimes = evaltimes(:);
				infset = ones(size(evaltimes)) * Inf;
				[gsltemp,gsltemp_V] = evalfunc([infset infset  evaltimes]);
				gsltemp_sd = sqrt(diag(gsltemp_V));
		end
		if haveRealData			
				sub = find(isinf(data.lat+data.long));
				if length(sub)>0
					[sorted,sorti] = sort(data.age_real(sub));
					sub2 = find((evaltimes>=mintime).*(evaltimes<=maxtime));
					[sortedu,sortedui]=unique(sorted);
					sortedui = sorti(sortedui);

					realsl = interp1(sortedu,data.sl_real(sub(sortedui)),evaltimes(sub2),'linear','extrap');
					chiGSL(totalcount) = (1/length(sub2)) * sum( ((realsl(:)-gsltemp(sub2))./gsltemp_sd(sub2)).^2);
				else
					chiGSL(totalcount) = 0;
				end

			chiSL(totalcount) = (1/length(data.sl_real)) * sum( ((data.sl_real(:) - fold(:))./f_sdold(:)).^2);
			chiAge(totalcount) = (1/length(realt)) * sum( ((realt(:)-tempt(:))./(perturbablet_sd(:)+1e-4)).^2);


		end
		if showplot
			clf; 
			subplot(3,1,1)
			meanoff = [];
			if haveRealData
				errorxy([data.sl_real fold f_sdold]);
				l1 = xlim; l2 = ylim;
				m1 = min([data.sl_real(:) ; fold]); m2 = max([data.sl_real(:) ; fold]);
				hold on ; plot([m1 m2],[m1 m2],'--');
				xlim(l1); ylim(l2);

			title([num2str(counter) ' / ' num2str(i) ' -- \chi^2/N  = ' num2str(chiSL(totalcount)) '  -- logp/N=' num2str(logp/length(fold)) ' [' num2str(PosteriorSL) '] [' num2str(logpReal/length(fold)) ' / ' num2str(logpInitial/length(fold)) ']' ]);
			else
				errorxy([g(:,i) fold f_sdold]);
				title([num2str(counter) ' / ' num2str(i) '  -- logp/N=' num2str(logp/length(f(:,i)))]);
			end

			subplot(3,1,2)
			meanoff = [];
			if haveRealData
				errorxy([data.age_real(:) gold perturbablet_sd(:)]);
				l1 = xlim; l2 = ylim;
				m1 = min([data.age_real(:) ; gold]); m2 = max([data.age_real(:) ; gold]);
				hold on ;plot([m1 m2],[m1 m2],'--');
				xlim(l1); ylim(l2);
				title([num2str(counter) ' / ' num2str(i) ' -- \chi^2/N = ' num2str(chiAge(totalcount)) '  -- logp/N=' num2str(logp/length(fold)) ' [' num2str(PosteriorAge) '] [' num2str(logpReal/length(fold)) ']']);
			else
				errorxy([gold fold data.age_sd],'ColXe',3,'ColYe',NaN);
			end
			
			subplot(3,1,3)
				errorxy([evaltimes gsltemp gsltemp_sd]);
				if haveRealData
					sub = find(isinf(data.lat+data.long));
					[sorted,sorti] = sort(data.age_real(sub));
					hold on; plot(sorted,data.sl_real(sub(sorti)));
					mintime = min(data.age_real(sub)); maxtime = max(data.age_real(sub));
					title(['\chi^2/N = ' num2str(chiGSL(totalcount))]);
					ylim([min(-25,min(gsltemp)) max(10,max(gsltemp))]);
				end
			refresh
			pause(.1)
		end
		
		%%% end monitoring routines
		
		
		accepted(:,i)=zeros(length(g(:,i)),1);
		% we do this in a random order each time
		
		sub = find(perturbablet_sd>timeres); % we don't bother evaluating points with uncertainties below resolution
		
		if length(sub)>0
			for k=randsample(sub(:)',length(sub))			
				gtrial = gold;
				
				% calculate a trial step, taking into account sequence constraints
				minoffset = -Inf;
				maxoffset = Inf;
				if perturbablet_predecessor(k) > 0
					maxoffset = gold(perturbablet_predecessor(k)) - gold(k);
				end
				if perturbablet_successor(k) > 0
					minoffset = gold(perturbablet_successor(k)) - gold(k);
				end

				[perturbs(k),backtofwd] = proposalDraw(perturbablet_sd(k)*MCStepSize(k),minoffset,maxoffset);
				gtrial(perturbablet_succset{k}) = gold(perturbablet_succset{k}) + perturbs(k);

				% calculate likelihood of trial observation vs. old observation
				if ismember(k,sequenceControlled)
					obsLikelihoodTrial = pdf('norm', perturbablet(k),gtrial(perturbablet_predecessor(k)) - gtrial(k),perturbablet_sd(k)/dataweight);
					obsLikelihoodOld = pdf('norm',perturbablet(k),gold(perturbablet_predecessor(k)) - gold(k),perturbablet_sd(k)/dataweight);
				else
					obsLikelihoodTrial = pdf('norm',t(k),gtrial(k),t_sd(k)/dataweight);
					obsLikelihoodOld = pdf('norm',t(k),g(k,i),t_sd(k)/dataweight);
				end
				if perturbablet_predecessor(k) > 0
					obsLikelihoodOld = obsLikelihoodOld * (gold(k) < gold(perturbablet_predecessor(k))) ;
					obsLikelihoodTrial = obsLikelihoodTrial * (gtrial(k) < gtrial(perturbablet_predecessor(k)));
				end
				if perturbablet_successor(k) > 0
					obsLikelihoodOld = obsLikelihoodOld * (gold(perturbablet_successor(k)) < gold(k));
					obsLikelihoodTrial = obsLikelihoodTrial * (gtrial(perturbablet_successor(k)) < gtrial(k));;
				end
				
				% calculate the likelihood of the new assemblage of ages, given f
				[ftrial,Vtrial,modelLogLikelihoodTrial] = GPSpaceTimeRegression(f(:,i),data.lat,data.long,gtrial,f_sd(:,i),cvfile,tprsc);
%				f_sdtrial = sqrt(diag(Vtrial));

				% calculate the acceptance probability
				acceptanceProb = backtofwd * exp(modelLogLikelihoodTrial-modelLogLikelihoodOld) * obsLikelihoodTrial /  obsLikelihoodOld;
				
%				disp(['              ' num2str([obsLikelihoodTrial/obsLikelihoodOld exp(modelLogLikelihoodTrial-modelLogLikelihoodOld)]) ]);
				acceptanceProbRecord(k,totalcount) = min(max(1e-3,acceptanceProb),1);
				if ~isfinite(acceptanceProbRecord(k,totalcount))
					acceptanceProbRecord(k,totalcount) = 1;
				end
				
				% use a random number to determine whether or not to accept, and display comment
				if (rand < acceptanceProb) || (isnan(acceptanceProb))
					dispstring=sprintf(['     Accepted at position ' num2str(k) ' (p = ' num2str((acceptanceProb)) ', dt = ' num2str(perturbs(k)) ')\t[Step Size: ' num2str(MCStepSize(k)) ']']);
					if totalcount > 1
						dispstring = [dispstring ' [' num2str(round(100*mean([acceptedRecord(k,:) 1]))) '% / ' num2str(mean(acceptanceProbRecord(k,end-min(totalcount-1,adaptationPeriod):end))) ']'];
					end
					disp(dispstring);
					accepted(k,i) = 1;
					gold = gtrial; 
%						fold = ftrial; f_sdold = f_sdtrial; Vold = Vtrial;
					modelLogLikelihoodOld = modelLogLikelihoodTrial;
				else
					dispstring=(sprintf(['     Rejected at position ' num2str(k) ' (p = ' num2str((acceptanceProb)) ', dt = ' num2str(perturbs(k)) ')\t[Step Size: ' num2str(MCStepSize(k)) ']']));
					if totalcount > 1
						dispstring = [dispstring ' [' num2str(round(100*mean([acceptedRecord(k,:) 0]))) '% / ' num2str(mean(acceptanceProbRecord(k,end-min(totalcount-1,adaptationPeriod):end))) ']'];
					end
					disp(dispstring);
	
				end
				
				
				% early on in the Monte Carlo simulation, we adaptatively adjust the step size
				% for each measurement
				
				if mod(totalcount,adaptationPeriod)==0
					adaptiveStep = 1+ (adaptiveStep0-1) * exp(-((totalcount/adaptationPeriod)/(adaptiveDecayScale)));
					if adaptiveStep > 1.01
						meanaccept = mean(acceptanceProbRecord(k,end-(adaptationPeriod-1):end));
						if meanaccept/idealAcceptanceProb > 1.1
							MCStepSize(k) = MCStepSize(k) * adaptiveStep;
						elseif meanaccept/idealAcceptanceProb < 0.9
							MCStepSize(k) = MCStepSize(k) / adaptiveStep;
						end
						MCStepSize(k) = min(max(MCStepSize(k),adaptiveMin),adaptiveMax);
					end
				end
				
				
				
				
				
			end
		end
		disp(['Acceptance rate: ' num2str(mean(accepted(:,i)))]);
		acceptedRecord(:,totalcount) = accepted(:,i);
		totalcount = totalcount + 1;
		g(:,i) = gold;
		
		g(:,i+1) = g(:,i);

	end
	yKeep(:,counter) = y(:,i);
	fKeep(:,counter) = f(:,i);
	gKeep(:,counter) = g(:,i);
	f_sdKeep(:,counter) = f_sd(:,i);
	f_VKeep(:,:,counter) = f_V(:,:,i);
	acceptanceRateKeep(:,counter) = mean(accepted,2);
	
	if exist('chiSL')
		store.chiSL = chiSL;
	end
	if exist('chiAge')
		store.chiAge = chiAge;
	end
	if exist('chiGSL')
		store.chiGSL = chiGSL;
	end

	save([uniqueid '.mat'],'data','N','thin','fKeep','gKeep','f_sdKeep','acceptanceRateKeep','ageOffset','store');
	

end

fid = fopen(['status.' uniqueid],'w');
fprintf(fid,['Ended\t ' 'of ' num2str(runcycles+1) '\t' datestr(now) '\t' uniqueid '\n']);
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = uncensorData(means,sds,mins,maxs)

% y=uncensorData(means,sds,mins,maxs)
%
% Takes random samples from truncated normal distributions.
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 17 May 2009

if (length(mins) == 1) && (length(means)>1)
	mins = repmat(mins,size(means));
end
if (length(maxs) == 1) && (length(means)>1)
	maxs = repmat(maxs,size(means));
end

cdfmins = cdf('norm',mins,means,sds);
cdfmaxs = cdf('norm',maxs,means,sds);
cdfranges = cdfmaxs-cdfmins;

rands = rand(size(means));
rands = (rands .* cdfranges) + cdfmins;
y = icdf('norm',rands,means,sds);

% look for long tails, which will give rise to NaNs if we're not careful

exphandler = find(cdfranges==0);
if length(exphandler) > 0
	seta = intersect(exphandler, find(mins>means));
	setb = intersect(exphandler, find(maxs<means));
	
	cdfranges(seta) = cdf('norm',means(seta)-mins(seta),0,sds(seta));
	rands(seta) = rand(size(seta));
	rands(seta) = rands(seta) .* cdfranges(seta);
	y(seta) = means(seta) + (means(seta) - icdf('norm',rands(seta),means(seta),sds(seta)));
	seta2 = intersect(seta,find(isinf(y)));
	y(seta2) = mins(seta2);
	
	cdfranges(setb) = cdfmaxs(setb);
	rands(setb) = rand(size(setb));
	rands(setb) = rands(setb) .* cdfranges(setb);
	y(setb) = max(mins(setb),icdf('norm',rands(setb),means(setb),sds(setb)));
	setb2 = intersect(setb,find(isinf(y)));
	y(setb2) = maxs(setb2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y,backprobtofwdprob] = proposalDraw(stdev,mindraw,maxdraw,N)

% [y,backprobtofwdprob] = proposalDraw(stdev,mindraw,maxdraw,N)
%
% Generate a proposal step, subject to the bounds mindraw and maxdraw.
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 17 May 2009

defval('mindraw',-Inf); defval('maxdraw',Inf); defval('N',1);
probscalefwd = 1; probscaleback = 1;

r = rand(N,1);

if isfinite(maxdraw) || isfinite(mindraw)
	mincdf = cdf('norm',mindraw,0,stdev);
	maxcdf = cdf('norm',maxdraw,0,stdev);
	range = maxcdf-mincdf;
	y = icdf('norm',(r .* range) + mincdf,0,stdev);
	
	if ~isfinite(y)
		if isfinite(maxdraw)
			y = maxdraw;
		elseif isfinite(mindraw)
			y = mindraw;
		else
			y = icdf('norm',r,0,stdev);
		end
	end
	fwdprob = pdf('norm',y,0,stdev)./range;
	if isnan(fwdprob)
		fwdprob = Inf;
	end
			
	mindraw2 = mindraw - y; maxdraw2 = maxdraw - y;
	mincdf2 = cdf('norm',mindraw2,0,stdev);
	maxcdf2 = cdf('norm',maxdraw2,0,stdev);
	range2 = maxcdf2-mincdf2;
	backprob = pdf('norm',-y,0,stdev)./range2;
	if isnan(backprob)
		backprob = Inf;
	end
	
	backprobtofwdprob = backprob./fwdprob;
	if isnan(backprobtofwdprob)
		backprobtofwdprob = 1;
	end
else
	y = icdf('norm',r,0,stdev);
	backprobtofwdprob = 1;
end
