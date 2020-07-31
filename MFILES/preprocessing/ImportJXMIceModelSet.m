function y=ImportJXMIceModelSet(fnames,minage,maxage)
% ImportJXMIceModelSet(fnames,minage,maxage)
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 5 May 2009

defval('fnames','ice.*');
defval('minage',-28);
defval('maxage',0);

% constants
rad = 6.371e6; % radius of the Earth
RHOE=5511.57; RHOW=1000.0; RHOI=920.0; %densities
areaOcean = 361419000 * 1e6;

u=strfind(fnames,'/');
if length(u)>0
	fpath = fnames(1:u(end));
else
	fpath = '';
end

files = dir(fnames);
for i=1:length(files)
	age(i) = ImportJXMIceModel(fullfile(fpath,files(i).name));
end

[age,ind] = sort(-age);
files=files(ind);
sub=find((age<=maxage).*(age>=minage));

files=files(sub); age=age(sub);
y.timesteps = age;

if length(files)>0
	i=1;
	[ag,ic,lat,long] = ImportJXMIceModel(fullfile(fpath,files(i).name));
	y.ice(:,i) = ic(:);
	[mshlong,mshlat] = meshgrid(long,lat);
	y.long=mshlong(:); y.lat=mshlat(:); clear mshlong mshlat;
	
	for i=2:length(files)
		[ag,ic] = ImportJXMIceModel(fullfile(fpath,files(i).name));
		y.ice(:,i) = ic(:);
	end
end

midlat = [90 (lat(2:end)+lat(1:end-1))/2 -90];

longsize = long(2)-long(1);
y.volfactor = repmat(rad^2 * (longsize * pi/180).*(sind(midlat(1:end-1))-sind(midlat(2:end))),1,length(long));
y.volfactor = y.volfactor(:);
vol = repmat(y.volfactor,[1 size(y.ice,2)]) .* y.ice;
y.esl = (RHOI/RHOW) * vol / areaOcean;
	
		
icemasses = SeparateIceMasses(y.lat,y.long);
y.sheets = unique(icemasses);
for i=1:length(y.sheets)
	sub = find(icemasses == y.sheets(i));
	y.sheetesl(i,:) = (sum(y.esl(sub,:)));
end
y.sheetesl = full(y.sheetesl);


end
