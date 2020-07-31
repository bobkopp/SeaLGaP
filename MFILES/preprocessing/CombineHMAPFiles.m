function [hmaps,hmapids,hmaptimes,sheetesl,long,lat]=CombineHMAPFiles(fnames,savefile)

% [hmaps,hmapids,hmaptimes,sheetesl,long,lat]=CombineHMAPFiles(fnames,savefile)
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 16 May 2009

defval('savefile','hmaps-combined.mat');

hmaps = [];
hmapids = [];
hmaptimes = [];
sheetesl = [];
long = [];
lat = [];

files = dir(fnames);
for i=1:length(files)
	infile = load(files(i).name);
	hmaps = [hmaps infile.hmaps];
	hmapids = [hmapids infile.hmapids];
	hmaptimes = [hmaptimes infile.hmaptimes];
	sheetesl = [sheetesl infile.sheetesl];
	long = infile.long;
	lat = infile.lat;
end

save(savefile,'hmaps','hmapids','hmaptimes','sheetesl','lat','long');