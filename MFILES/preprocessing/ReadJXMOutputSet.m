function [hmaps,hmapids,hmaptimes,sheetesl,steric,gsl,long,lat]=ReadJXMOutputSet(basemodels,slong,slat)

% EXAMPLE:
%      basemodels(1)=ImportJXMIceModelSet('ice-models/ice5g/ice.*',-26,0);
%      basemodels(2)=ImportJXMIceModelSet('ice-models/bassett/ice.*',-28,0);
%      sldb=SealevelDBImport('sealevel11.csv',[-Inf Inf],'keepall');
%      [ulatlong,ulli] = unique((sldb.lat)*1000+(sldb.long));		
%      sub = find(isfinite(ulatlong));
%      ulli = ulli(sub);
%      slat = (sldb.lat(ulli)); slong=(sldb.long(ulli));
%      [hmaps,hmapids,hmaptimes,sheetesl,long,lat]=ReadJXMOutputSet(basemodels,slong,slat)
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 31 July 2009

defval('synarcdiro','ARCHIVE-OUTPUT');
defval('synworkdiro','WORK');
defval('latbase',-90:5:90);
defval('longbase',0:10:360);
defval('bandwidth',256);
defval('bandwidthtrunc',256); 
defval('ageref',0); % the last column is actually t=+1 ky sea level
defval('minage',100);
defval('maxage',150);
defval('maxtoload',500);
defval('runstoskip',0);
defval('slong',[]);
defval('slat',[]);
defval('savefile','hmaps');
defval('seedtimes',150:-1:90);

startdir = pwd;
synarcdiro = fullfile(startdir,synarcdiro);
synworkdiro = fullfile(startdir,synworkdiro);
savefile = fullfile(startdir,savefile);

basemodelind = @(j) (1+mod(j-1,length(basemodels)));
for j=1:length(basemodels)
	basesheetesl(:,j) = basemodels(j).sheetesl(:,end);
	basesheeteslstart(:,j) = basemodels(j).sheetesl(:,1);
end

%load seeds

[LONG,LAT] = meshgrid(longbase,latbase);
long = LONG(:); lat = LAT(:);

[ulonglat,ulli]=unique(slat*1e3+slong);
sub=find(isfinite(ulonglat));
slat=slat(:); slong=slong(:);
long = [long ; slong(ulli(sub))]; lat=[lat ; slat(ulli(sub))];

% load data

if ~exist(synworkdiro,'dir')
	mkdir(synworkdiro);
end

hmaps = [];
hmapids = [];
hmaptimes = [];
sheetesl = [];
steric = [];
gsl = [];

files=dir(fullfile(synarcdiro,'out_*.tar.gz'));

cd(synworkdiro);

for i=1:min(maxtoload,length(files))

	% extract to work directory
	[stat,rslt]=system(['tar -xvzf ' fullfile(synarcdiro,files(i).name)]);
	slashes=strfind(rslt,'/');
	workingpath='';
	if length(slashes)>0
		workingpath = rslt(1:(slashes(1)-1));
	end
	
	system(['rm -f ' workingpath '/seagl_1.*']);
	
	s=regexp(files(i).name,'out_(....)\.tar.gz','tokens');
	if length(s{1})==1
		id(i) = str2num(s{1}{1});
		disp(['PROCESSING JOB ' num2str(id(i)) '........']);
		if id(i)>runstoskip
			if exist('Plm','var')
				[curhmap,lo,la,timesteps,Plm,hmapref]=readjxms(fullfile(workingpath,'seagl_*.mn'),minage,maxage,ageref,[long lat],bandwidthtrunc,Plm);
			else
				[curhmap,lo,la,timesteps,Plm,hmapref]=readjxms(fullfile(workingpath,'seagl_*.mn'),minage,maxage,ageref,[long lat],bandwidthtrunc);
			end
%			model = ConstructJXMIceModel(basemodels(basemodelind(id(i))),timesteps,seeds.iced(:,:,id(i)));
%			model.steric = seeds.steric(i,:);
%			model.visco = seeds.visco(i,:);
%			model.gsl = seeds.gsl(i,:);
%			
%			% here we adjust for different starting ice volumes
%			% deltagslstart = sum(model.sheetesl(:,1)) - sum(basesheeteslstart);
%			% curhmap = curhmap - deltagslstart;
%			
%			se = model.sheetesl - repmat(basesheetesl(:,basemodelind(id(i))),[1 size(model.sheetesl,2)]);
%			
			hmaps=[hmaps curhmap]; hmapids = [hmapids repmat(id(i),[1 length(timesteps)])];
			hmaptimes=[hmaptimes timesteps];
%			sheetesl = [sheetesl se];
%			steric = [steric model.steric];
%			gsl = [gsl model.gsl];

			model = ConstructJXMIceModel(basemodels(basemodelind(id(i))),[],seeds.iced(:,:,id(i)));
			
			% here we adjust for different starting ice volumes
			% deltagslstart = sum(model.sheetesl(:,1)) - sum(basesheeteslstart);
			% curhmap = curhmap - deltagslstart;
			
			se = interp1(-model.timesteps,model.sheetesl',timesteps)' - repmat(basesheetesl(:,basemodelind(id(i))),[1 length(timesteps)]);
			
			sheetesl = [sheetesl se];
			steric = [steric interp1(seedtimes,seeds.steric(id(i),:),timesteps,'linear')];
			gsl = [gsl interp1(seedtimes,seeds.gsl(id(i),:),timesteps,'linear')];
	


			
			if length(savefile)>0
				save(savefile,'hmaps','hmapids','hmaptimes','sheetesl','steric','gsl','lat','long');
			end
			
		end
	end
	
	% clean up	
	
	system(['rm -fr ' fullfile(synworkdiro,workingpath)]);
end

cd(startdir);