function [gsl,icedist,steric,visco]=MatchHMAPSeeds(hmapids,hmaptimes,seeds,seedtimes)

% [gsl,icedist,steric,visco]=MatchHMAPSeeds(hmapids,hmaptimes,seeds,seedtimes)
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, July 3 2009

defval('seedtimes',150:-1:90);

st=hmapids(2:end)~=hmapids(1:end-1);
uids=[hmapids(find(st)) hmapids(end)];
utimes=hmaptimes(find(hmapids==uids(1)));

[c,ia,ib]=intersect(utimes,seedtimes);
subtimes=sort(ib);

gsl=seeds.gsl(uids,subtimes)';
gsl=gsl(:)';

icedist=reshape(seeds.iced(:,subtimes,uids),size(seeds.iced,1),[]);

steric=seeds.steric(uids,subtimes)';
steric=steric(:)';

visco=seeds.visco(uids);