function [ac,ac90,minac90,lastac90]=AutocorrelationCheck(datfiles,Nburnin)

% [ac,ac90,minac90,lastac90]=AutocorrelationCheck(datfiles,[Nburnin])
%
% Checks SeaLGaP output for autocorrelation.
%
% INPUTS:
%	datfiles	A structure outputted by dir() with data file names in
%				datfiles.name
%	Nburnin		The number of iterations to skip (Default = 0)
%
% OUTPUTS:
%	ac			Cell array of autocorrelation matrices.
%				Rows of ac correspond to different data points.
%				Columns correspond to increasing lags, with the first column
%				corresponding to a lag of 0.
%	ac90		Cell array of the 90th percentile value for each lag.
%	minac90		The minimum of ac90 for each datfile.
%	lastac90	The last value of ac90 for each datfile.
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 19 May 2009

defval('Nburnin',0);

doPlot = 1;
for i=1:length(datfiles)
	disp(datfiles(i).name);
	dat = AssimilateRun(datfiles(i),Nburnin);
	ac{i} = autocorr(dat.gKeep);
	len(i) = size(ac{i},2);
	sub = find(isfinite(sum(ac{i},2)));
	ac90{i} = quantile(ac{i}(sub,:),.9);
	lastac90(i) = ac90{i}(end);
	[minac90(i),mini(i)] = min(ac90{i});
	disp(['   Last(AC90): ' sprintf('%5.3f',lastac90(i))]);
	disp(['   Min(AC90): ' sprintf('%5.3f',minac90(i)) ' at lag of ' num2str(mini(i)-1)]);
	
end

%%%%%%%%%%%%%

function y=autocorr(X)

% y=autocorr(X)
%
% Returns an autocorrelation matrix for X, where different variables
% are in the rows of X and different equally spaced times in the columns
% of X.
%
% In y, rows correspond to the different variables and columns to different
% lags. Column 1 corresponds to a lag of 0, and should equal ones(size(X,1),1).
%
% Check proceeds up to a lag of floor(size(X,2)/2).
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 19 May 2009

if sum(size(X)==1)
	X=X(:)';
end

Nvars = size(X,1);
Ntimes = size(X,2);

y=zeros(size(X,1),floor(Ntimes/2));

% i is lag

for i=0:floor(Ntimes/2)

	for j=1:Nvars
		temp=[X(j,1:(end-i)) ; X(j,(i+1):end)]';
		c = corr(temp);
		csub = c(1,2);
		y(j,(i+1)) = (csub);
	end

	% display some sign of progress
	if mod(i,20) == 0
		ac90 = quantile(y(:,(i+1)),.9);
		disp(['   Lag ' num2str(i) '/' num2str(floor(Ntimes/2)) '     [AC90 = ' sprintf('%4.2f',ac90) ']']);
	end

end
