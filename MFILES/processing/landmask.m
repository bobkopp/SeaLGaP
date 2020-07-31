function [y,inpoly]=landmask(testlat,testlong)

% [y,inpoly]=landmask(testlat,testlong)
%
% modified from PLOTCONT
%
%
% Saved matrix as space-saving unsigned integer 
% - but that translates the NaN's into some  high number - take that out.
% Note how A==NaN does not return the NaN's!
%
% You'll have to maintain your own databases of continental outlines!
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 24 July 2009

defval('ofs',0);

% Where are the data kept? Create a directory $IFILES/COASTS
% with $IFILES set as an environment variable...
defval('ddir',fullfile(getenv('IFILES'),'COASTS'))


testlat=testlat(:);
testlong=testlong(:);

% Load the data sets
fid=fopen(fullfile(ddir,'cont.mtl'),'r','b');
cont=fread(fid,[5217 2],'uint16');
fclose(fid);

% Recast data in good form
cont=cont/100-90;
cont(cont==max(max(cont)))=NaN;

% Collect the data
lon=cont(:,1); lat=cont(:,2);
% To use fill, eliminate all the NaNs
fdem=find(isnan(lon));
% Just need to know a lot about this data
want=fdem(22)+1:fdem(23)-1;
eant=fdem(95)+1:fdem(96)-1;
% Redo the data set by fixing England and adding Antarctica last, etc
lon=[lon([1:fdem(2)-1 fdem(2)+1:fdem(22) fdem(23)+1:fdem(95) ...
   fdem(96)+1:fdem(104)-1 fdem(104)+1:end]) ; ...
     NaN; lon(eant) ; lon(want) ; ...
     360  ; 360 ; 0 ; 0 ; NaN];
lat=[lat([1:fdem(2)-1 fdem(2)+1:fdem(22) fdem(23)+1:fdem(95) ...
   fdem(96)+1:fdem(104)-1 fdem(104)+1:end]) ; ...
     NaN; lat(eant) ; lat(want) ; ...
     lat(want(end)) ; -90 ; -90 ; lat(eant(1)) ; NaN];
% And partitioning again
fdem=find(isnan(lon));
% Take out last one
tri=101;
lon=lon([1:fdem(tri)-1 fdem(tri+1)+1:end]);
lat=lat([1:fdem(tri)-1 fdem(tri+1)+1:end]);
XYZ=[lon+ofs lat];
% And partitioning again
fdem=find(isnan(lon));
% Now start with the patches
beg=1; handl=nan(1,length(fdem));
hold on
for i=1:length(fdem)
  inpoly(:,i)=inpolygon(testlong,testlat,lon(beg:fdem(i)-1),lat(beg:fdem(i)-1));
%  handl(i)=fill(lon(beg:fdem(i)-1),lat(beg:fdem(i)-1),pcol);
  beg=fdem(i)+1; 
end
y=sum(inpoly,2)>0;

