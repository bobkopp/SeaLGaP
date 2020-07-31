% wrapper MATLAB script for running SealevelAnalyze
% Last updated by Bob Kopp rkopp-at-princeton.edu, 8 July 2009

clear;

% get some parameters from the environment
anchorseq = getenv('ANCHORSEQ');
cvfile = getenv('CVFILE');
IFILES = getenv('SEALGAP_IFILES'); defval('IFILES',getenv('IFILES'));
maxage = str2num(getenv('MAXAGE'));
minage = str2num(getenv('MINAGE'));
N = str2num(getenv('MONTECARLO_N'));
thin = str2num(getenv('MONTECARLO_THIN'));
params = getenv('PARAMS');
scriptid = getenv('SCRIPTID');
addtopath = getenv('SEALGAP_PATH');
progpath = getenv('SEALGAP_ROOT');
slfile = getenv('SEALEVEL_FILE');
slvar = getenv('SEALEVEL_VAR');
workpath = getenv('WORKPATH');
dataweight = str2num(getenv('DATAWEIGHT'));
tprsc=str2num(getenv('TPRSC'));

% and some defaults if they're not specified
defval('scriptid','slanalysis');
defval('N',7500);
defval('thin',20);
defval('progroot',pwd);
defval('workpath',pwd);
defval('IFILES',fullfile(progroot,'IFILES'));
defval('slfile','sldb');
defval('slvar','sldb');
defval('maxage',160); defval('minage',80);
defval('anchorseq',1);
defval('cvfile','cvtable');
defval('dataweight',1);
%defval('tprsc',[]);

% set some paths
setenv('IFILES',IFILES);
addpath(IFILES);

addpath(progroot);

if length(addtopath)>0
	eval(['addpath ' addtopath]);
end

if ~exist(workpath,'dir')
	mkdir(workpath);
end
cd(workpath);


% split params based on ':' delimiter

nogsl = 0; nolongseqs=0; groupfilter=0; categoryfilter = 0; categorysubset = 0;

if length(params)>0
	delims = [0 strfind(params,':') (length(params)+1)];
	for i=1:length(delims)-1
		param{i} = params((delims(i)+1):(delims(i+1)-1));
	end
	
	for i=1:length(param)
		if strcmpi(param{i},'nogsl')
			nogsl = 1;
			disp('Using local sea level terms only');
			scriptid = [scriptid '_nogsl'];
		elseif strcmpi(param{i},'nolongseqs')
			nolongseqs = 1;
			disp('Eliminating all long sequences');
			scriptid = [scriptid '_nolong'];
		elseif strcmpi(params,'groupfilter')
			groupfilter = 1;
			disp('Group filter mode');
			scriptid = [scriptid '_groupfilter'];
			loadfile = [IFILES '/' 'groupsAnalyzed.mat'];
			touchfile = [IFILES '/' 'groupsAnalyzed.touch.mat'];
		elseif strcmpi(params,'categoryfilter')
			categoryfilter = 1;
			disp('Category filter mode');
			scriptid = [scriptid '_categoryfilter'];
			loadfile = [IFILES '/' 'catsAnalyzed.mat'];
			touchfile = [IFILES '/' 'catsAnalyzed.touch.mat'];
		elseif strcmpi(params,'categorysubset')
			categorysubset = 1;
			disp('Category subset mode');
			scriptid = [scriptid '_categorysubset'];
			loadfile = [IFILES '/' 'catsubsAnalyzed.mat'];
			touchfile = [IFILES '/' 'catsubsAnalyzed.touch.mat'];
		end
	end
end

processid = getenv('PBS_JOBID'); mpiexecid = 0;
PBSrun = 0;

strsub=find(slvar=='{');
if length(strsub)>0
	strsub=strsub(1);
else
	strsub=length(slvar)+1;
end
strsub=1:(strsub-1);
s = load(slfile,slvar(strsub));
%dat = s.(slvar); clear s;
dat = eval(['s.' slvar]); clear s;

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
	maxNumCompThreads(1);
	
	processid = getenv('PBS_JOBID');
	mpiexecid = getenv('MPIEXEC_RANK');

	uniqueid = [scriptid '.' hostname '.' num2str(processid)]

	if length(mpiexecid) > 0
		uniqueid = [uniqueid '.' num2str(mpiexecid)];
	else
		mpiexecid = 0;
	end
end

 
% load map
map=dat;

% find anchor sequence and toss out data way away from time interval of interest
subset=(find(map.siteid==anchorseq));
[m,i] = min(map.stratseqid(subset));
map.age_sd(subset(i)) = 0;
subset2 = intersect(find((map.age_mode-map.age_sd)<=maxage),find((map.age_mode+map.age_sd)>=minage));
subset=union(subset,subset2);
subset=intersect(subset,find(map.usable));

map=MapSubset(map,subset);

tossout = [];
if nogsl
	subset = find(isfinite(map.lat));
	map=MapSubset(map,subset);
end

if nolongseqs
	%%%
	% toss out long sequences
	u = unique(map.siteid);
	for i=1:length(u)
		maxsteps(i) = max(map.stratseqid(find(map.siteid==u(i))));
		if maxsteps(i) > 9
			tossout = [tossout u(i)];
		end
	end
	%%%
	subset = find(~ismember(map.siteid,tossout));
	map=MapSubset(map,subset);
end

% adjust all the sequences so that they start at 1
u = unique(map.siteid);
for i=1:length(u)
	[m,k] = min(map.stratseqid(find(map.siteid==u(i))));
	map.stratseqid(find(map.siteid==u(i))) = map.stratseqid(find(map.siteid==u(i))) + 1 - m;
end

% I forget why I did this, but I think it has to do with invertibility
correctingerror = 1e-3;
map.sl_sd = map.sl_sd+correctingerror;


% Group filtering mode
if groupfilter
	while exist(touchfile)
		pause(5)
	end
	
	touch = 1; save(touchfile,'touch');
	if ~exist(loadfile)
		possibleGroups = unique(map.ResearchGroup);
		groupsAnalyzedTimes = zeros(size(possibleGroups));
		save(loadfile,'possibleGroups','groupsAnalyzedTimes');
	else
		load(loadfile)
	end
	delete(touchfile)
	
	unanalyzed = possibleGroups( find(groupsAnalyzedTimes==min([groupsAnalyzedTimes(:)'])));
	pickedGroup = unanalyzed(randsample(length(unanalyzed),1));
	
	map = MapSubset(map,find(map.ResearchGroup~=pickedGroup));
	
	disp(['Excluding Research Group ' num2str(pickedGroup)]);
	
	scriptid = [scriptid '_without' num2str(pickedGroup)];
	
	while exist(touchfile)
		pause(5)
	end
	touch = 1;
	save(touchfile,'touch');
	load(loadfile);
	groupindex = find(possibleGroups == pickedGroup);
	groupsAnalyzedTimes(groupindex) = groupsAnalyzedTimes(groupindex) + 1;
	save(loadfile,'possibleGroups','groupsAnalyzedTimes');
	delete(touchfile)
end

% Category filtering/subset mode
if categoryfilter || categorysubset

	% classify
	isotopic = zeros(size(map.sl_mode));
	constructional = zeros(size(map.sl_mode));
	coral = zeros(size(map.sl_mode));
	erosional = zeros(size(map.sl_mode));
	facies = zeros(size(map.sl_mode));

	for i=1:length(map.sl_mode)
		if strfind(map.MaterialCategory{i},'isotop')
			isotopic(i) = 1;
		end
		if strfind(map.MaterialCategory{i},'constr')
			constructional(i) = 1;
		end
		if strfind(map.MaterialCategory{i},'coral')
			coral(i) = 1;
		end
		if strfind(map.MaterialCategory{i},'eros')
			erosional(i) = 1;
		end
		if strfind(map.MaterialCategory{i},'fac')
			facies(i) = 1;
		end
	end
	
	u = unique(map.siteid);
	for i=1:length(u)
		sub = find(map.siteid==u(i));
		siteIsotopic(i) = sum(isotopic(sub));
		siteConstructional(i) = sum(constructional(sub));
		siteCoral(i) = sum(coral(sub));
		siteErosional(i) = sum(erosional(sub));
		siteFacies(i) = sum(facies(sub));
	end
	
	% consolidate coral and constructional and form into matrix
	siteTypes = [siteIsotopic ; siteConstructional + siteCoral; siteErosional ; siteFacies];
	catnames = {'Iso','Cor','Ero','Fac'};
	
	while exist(touchfile)
		pause(5)
	end
	
	touch = 1; save(touchfile,'touch');
	if ~exist(loadfile)
		catAnalyzedTimes = zeros([length(catnames) 1]);
		save(loadfile,'catAnalyzedTimes');
	else
		load(loadfile)
	end
	delete(touchfile)
	
	unanalyzed = ( find(catAnalyzedTimes==min([catAnalyzedTimes(:)'])));
	pickedCat = unanalyzed(randsample(length(unanalyzed),1));

	if categoryfilter		
		sub = find(siteTypes(pickedCat,:)==0);
		disp(['Excluding ' catnames{pickedCat} ' Category']);
		scriptid = [scriptid '_without' catnames{pickedCat}];
	elseif categorysubset
		sub = find(siteTypes(pickedCat,:)==1);
		disp(['Examining only ' catnames{pickedCat} ' Category']);
		scriptid = [scriptid '_only' catnames{pickedCat}];
	end	
	

	subset = find(ismember(map.siteid,u(sub)));
	
	map = MapSubset(map,subset);

	
	while exist(touchfile)
		pause(5)
	end
	touch = 1;
	save(touchfile,'touch');
	load(loadfile);
	catAnalyzedTimes(pickedCat) = catAnalyzedTimes(pickedCat) + 1;
	save(loadfile,'catAnalyzedTimes');
	delete(touchfile)
end


SealevelAnalyze(map,N,thin,'label',scriptid,'cvfile',cvfile,'dataweight',dataweight,'tprsc',tprsc,'maxage',maxage,'minage',minage)