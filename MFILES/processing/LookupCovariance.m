function [cv,mux1,mux2,sdx1,sdx2]=LookupCovariance(x1,x2,cvfile,taperscale,trunctaperscale,pctrunc,pcoveralltrunc)

% [cv,mux1,mux2,sdx1,sdx2] = LookupCovariance(x1,x2,[cvfile],[taperscale],[trunctaperscale])
%
% This version performs a look up for both space and time components.
%
% INPUTS
%	x1			[lat long age]
%	x2			[lat long age]
%	cvfile		lookup table (in directory given by IFILES, or current directory) (Default: 'cv.mat')
%	taperscale	taper time scale for values larger than the max in lookup table (Default: 3)
%
% OUTPUTS
%	cv			covariance matrix between specified points
%	mux1		mean values at x1
%	mux2 		mean values at x2
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 13 June 2009

defval('x2',x1);
defval('cvfile','cvtable.mat');
defval('taperscale',3);
%defval('trunctaperscale',1000);
defval('trunctaperscale',3);

defval('pctrunc',0);
defval('pcoveralltrunc',0);

needsloading = 0;

persistent mu sd coef cvscor muoverall sdoverall coefoverall cvscoroverall lat long tims cvfilelast;

if length(strfind(cvfile,'.'))==0
	cvfile = [cvfile '.mat'];
end

if exist('cvfilelast','var')
	if ~strcmpi(cvfile,cvfilelast)
		needsloading = 1;
	end
end
cvfilelast = cvfile;

IFILES = getenv('IFILES');
if length(IFILES) == 0
	IFILES = ['.'];
end

if strfind(cvfile,'/')
	loadfile=cvfile;
else
	loadfile = fullfile(IFILES,cvfile);
end

if ~exist(loadfile)
	disp(['ERROR: Lookup table ' loadfile ' not found!']);
	return;
end

if exist('cvscor')
	if prod(size(cvscor))==0
		needsloading = 1;
	end
else
	needsloading = 1;
end

if needsloading
	disp(['Loading covariance table ' loadfile]);
	load(loadfile);
	lat=lat(:); long=long(:);
	lat(find(lat==-Inf)) = -1e6;
	lat(find(lat==Inf)) = 1e6;
	long(find(long==Inf)) = 1e6;
	LAT(find(LAT==-Inf)) = -1e6;
	LAT(find(LAT==Inf)) = 1e6;
	LONG(find(LONG==Inf)) = 1e6;		
end

if pctrunc==0
	pctrunc=size(cvscor,1);
else
	pctrunc=min(size(cvscor,1),pctrunc);
end

if pcoveralltrunc==0
	pcoveralltrunc=size(cvscoroverall,1);
else
	pcoveralltrunc=min(size(cvscoroverall,1),pcoveralltrunc);
end

nearest1 = ones([1 size(x1,1)]);
nearest2 = ones([1 size(x2,1)]);

subFinite1 = find(isfinite(sum(x1,2)));
subFinite2 = find(isfinite(sum(x2,2)));

x1(subFinite1,1)=min(90,x1(subFinite1,1));
x2(subFinite2,1)=min(90,x2(subFinite2,1));
x1(subFinite1,1)=max(-90,x1(subFinite1,1));
x2(subFinite2,1)=max(-90,x2(subFinite2,1));

x1(subFinite1,2) = mod(x1(subFinite1,2),360);
x2(subFinite2,2) = mod(x2(subFinite2,2),360);


% reset infinities to 1e6
x1(find(x1(:,1)==Inf),1)  =  1e6;
x1(find(x1(:,1)==-Inf),1) = -1e6;
x1(find(x1(:,2)==Inf),2)  =  1e6;

x2(find(x2(:,1)==Inf),1)  =  1e6;
x2(find(x2(:,1)==-Inf),1) = -1e6;
x2(find(x2(:,2)==Inf),2)  =  1e6;

lat=lat(:); long=long(:); tims=tims(:);

% identify nearest geographic location

diff1x = repmat(x1(:,2),[1 length(long)]) - repmat(long',[size(x1,1) 1]);
diff1y = repmat(x1(:,1),[1 length(lat)])  - repmat(lat', [size(x1,1) 1]);
dist1 = sqrt(diff1x.^2 + diff1y.^2);
[mindist1,mindist1i] = min(dist1,[],2);
nearest1 = mindist1i;


diff2x = repmat(x2(:,2),[1 length(long)]) - repmat(long',[size(x2,1) 1]);
diff2y = repmat(x2(:,1),[1 length(lat)])  - repmat(lat', [size(x2,1) 1]);
dist2 = sqrt(diff2x.^2 + diff2y.^2);
[mindist2,mindist2i] = min(dist2,[],2);
nearest2 = mindist2i;

% okay, so now the coordinates of the geographic locations are in nearest1 and nearest2

timepenalty = ones(size(x1,1),size(x2,1));
cv=zeros(size(x1,1),size(x2,1));

% first, identify points that are out of temporal range
subOOR1 = find((x1(:,3)<min(tims(:)))+(x1(:,3)>max(tims(:))));
subOOR2 = find((x2(:,3)<min(tims(:)))+(x2(:,3)>max(tims(:))));

subNotOOR1 = setdiff(1:size(x1,1),subOOR1);
subNotOOR2 = setdiff(1:size(x2,1),subOOR2);

[T2,T1]=meshgrid(x2(:,3),x1(:,3));
TD=T1-T2;

timepenalty(subOOR1,subOOR2)=exp(-((TD(subOOR1,subOOR2)/taperscale).^2));
timepenalty(subOOR1,subNotOOR2)=exp(-((TD(subOOR1,subNotOOR2)/taperscale).^2));
timepenalty(subNotOOR1,subOOR2)=exp(-((TD(subNotOOR1,subOOR2)/taperscale).^2));

cv(subOOR1,subOOR2) = coefoverall(nearest1(subOOR1),1:pcoveralltrunc)*cvscoroverall(1:pcoveralltrunc,1:pcoveralltrunc)*coefoverall(nearest2(subOOR2),1:pcoveralltrunc)';
cv(subOOR1,subNotOOR2) = coefoverall(nearest1(subOOR1),1:pcoveralltrunc)*cvscoroverall(1:pcoveralltrunc,1:pcoveralltrunc)*coefoverall(nearest2(subNotOOR2),1:pcoveralltrunc)';
cv(subNotOOR1,subOOR2) = coefoverall(nearest1(subNotOOR1),1:pcoveralltrunc)*cvscoroverall(1:pcoveralltrunc,1:pcoveralltrunc)*coefoverall(nearest2(subOOR2),1:pcoveralltrunc)';

mux1(subOOR1) = muoverall(nearest1(subOOR1));
mux2(subOOR2) = muoverall(nearest2(subOOR2));

sdx1(subOOR1) = sdoverall(nearest1(subOOR1));
sdx2(subOOR2) = sdoverall(nearest2(subOOR2));

cv = cv.*timepenalty;

% ok, now we deal with the points that we can look up in the table
% we have to interpolate to account for time

difft1 = repmat(x1(subNotOOR1,3),[1 length(tims)]) - repmat(tims',[length(subNotOOR1) 1]);
[dtP1,mindtP1i] = min(abs(difft1)+1e6*(difft1<0),[],2);
tP1 = mindtP1i;
[dtN1,mindtP1i] = min(abs(difft1)+1e6*(difft1>0),[],2);
tN1 = mindtP1i;
dtP1 = dtP1+eps; dtN1 = dtN1+eps;
wtP1 = dtN1./(dtN1+dtP1); wtN1 = dtP1./(dtN1+dtP1);

if length(subNotOOR1)>0
	
	compcoef1=bsxfun(@times,coef(nearest1(subNotOOR1)+(tP1-1)*length(lat),:),wtP1) ...
	        + bsxfun(@times,coef(nearest1(subNotOOR1)+(tN1-1)*length(lat),:),wtN1);

	mux1(subNotOOR1) = mu(nearest1(subNotOOR1)+(tP1-1)*length(lat)) .* wtP1 + ...
	                   mu(nearest1(subNotOOR1)+(tN1-1)*length(lat)) .* wtN1;
	
	sdx1(subNotOOR1) = sd(nearest1(subNotOOR1)+(tP1-1)*length(lat)) .* wtP1 + ...
	                   sd(nearest1(subNotOOR1)+(tN1-1)*length(lat)) .* wtN1;

end

%%%

difft2 = repmat(x2(subNotOOR2,3),[1 length(tims)]) - repmat(tims',[length(subNotOOR2) 1]);
[dtP2,mindtP2i] = min(abs(difft2)+1e6*(difft2<0),[],2);
tP2 = mindtP2i;
[dtN2,mindtP2i] = min(abs(difft2)+1e6*(difft2>0),[],2);
tN2 = mindtP2i;
dtP2 = dtP2+eps; dtN2 = dtN2+eps;
wtP2 = dtN2./(dtN2+dtP2); wtN2 = dtP2./(dtN2+dtP2);

if length(subNotOOR2)>0
	compcoef2=bsxfun(@times,coef(nearest2(subNotOOR2)+(tP2-1)*length(lat),:),wtP2) ...
	        + bsxfun(@times,coef(nearest2(subNotOOR2)+(tN2-1)*length(lat),:),wtN2);
	
	mux2(subNotOOR2) = mu(nearest2(subNotOOR2)+(tP2-1)*length(lat)) .* wtP2 + ...
	                   mu(nearest2(subNotOOR2)+(tN2-1)*length(lat)) .* wtN2;
	
	sdx2(subNotOOR2) = sd(nearest2(subNotOOR2)+(tP2-1)*length(lat)) .* wtP2 + ...
	                   sd(nearest2(subNotOOR2)+(tN2-1)*length(lat)) .* wtN2;

	if length(subNotOOR1)>0
		timepenalty(subNotOOR1,subNotOOR2)=exp(-((TD(subNotOOR1,subNotOOR2)/(trunctaperscale)).^2));
		cv(subNotOOR1,subNotOOR2) = compcoef1(:,1:pctrunc)*cvscor(1:pctrunc,1:pctrunc)*compcoef2(:,1:pctrunc)';
		cv(subNotOOR1,subNotOOR2) = cv(subNotOOR1,subNotOOR2) .* timepenalty(subNotOOR1,subNotOOR2);
	end
end

% and finally turn into non-normalized form

[SDX2, SDX1]=meshgrid(sdx2,sdx1);
cv = cv .* SDX1 .* SDX2;
