function y=PathExceedanceProbability(x,randpaths,Nboot,Nbootstrapsamples)

% y=PathExceedanceProbability(x,randpaths,[Nboot],[Nbootstrapsamples])
%
% Given a set of random samples, calculate exceedance probabilities.
%
% INPUTS:
%	x	values to calculate exceedance probabilities for
%	randpaths	matrix of random samples
%	Nboot	number of times to run calculation for bootstrap
%	Nbootstrapsamples	number of bootstrap samples to draw per calculation
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 24 May 2009

	defval('Nboot',1);
	defval('Nbootstrapsamples',1000);
	
	if size(randpaths,1)>1
		maxes=max(randpaths);
		maxes=maxes(:)';
	else
		maxes=randpaths(:)';
	end
	
	if Nboot == 1
		maxcheck = repmat(x(:),[1 length(maxes)]);
		checks = repmat(maxes,[length(x) 1]) <= maxcheck;
		y = 1 - sum(checks,2)/length(maxes);
	else
		disp('Bootstrapping confidence intervals...');
		Nbootstrapsamples = min(Nbootstrapsamples,length(maxes));
		for i=1:Nboot
			if mod(i,50)==0
				disp(['    ' num2str(i) '/' num2str(Nboot)]);
			end
			sample = randsample(length(maxes),Nbootstrapsamples,true);
			y(:,i) = PathExceedanceProbability(x,maxes(sample));
		end
	end

end