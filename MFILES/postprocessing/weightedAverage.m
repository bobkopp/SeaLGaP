function [avg,popsd,meansd,netsd] = weightedAverage(X,Xsd,ndim,addlweights)

% [avg,popsd,meansd,netsd] = weightedAverage(X,Xsd,[ndim],[addlweights])
%
% Computes the weighted average over dimension ndim.
% Weights default to Xsd.^-2.
%
% OUTPUT:
%	avg		weighted average
%	popsd	weighted population standard deviation
%	meansd	weighted average of standard deviations
%	netsd	quadratic sum of popsd and meansd
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 18 May 2009

	defval('ndim',length(size(X)));
	defval('addlweights',ones(size(Xsd)));
		
	weights = (Xsd.^-2 .* addlweights);
	wtsum = sum(X.*weights,ndim);
	sumweights = sum(weights,ndim);
	avg = wtsum./sumweights;
	wtsumsq = sum(X.^2 .* weights,ndim);
	
	if length(X)>1
		popsd = sqrt ( (wtsumsq .* sumweights - wtsum.^2) ./ (sumweights.^2 - sum(weights.^2,ndim)) );
	else
		popsd = 0;
	end
	%meanvar = sum(Xsd.^2 .* weights,ndim)./sumweights;
	%meanvar = sum(ones(size(Xsd)),ndim) ./ sumweights;
	meansd = sum(Xsd .* weights,ndim) ./ sumweights;
	netsd = sqrt(popsd.^2 + meansd.^2);
end