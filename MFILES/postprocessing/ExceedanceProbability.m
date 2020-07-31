function [y,randpaths,bootstraps] = ExceedanceProbability(x,f,V,goodpts,Nsamples,filtrsize,bootstrapquantiles)
	
% ExceedanceProbability(x,f,V,[goodpts],[Nsamples],[filtrsize],[bootstrapquantiles])
%
% Calculates exceedance probabilities.
%
% INPUTS:
%	x	vector of levels at which to calculate the probability of exceedance
%	f	data matrix, with different samples in each column
%	V	either 2D matrix of standard deviations or
%		3D matrix of variance matrices
%	goodpts	matrix with 1s corresponding to values to be included
%			and 0s indicating values to be excluded
%	Nsamples	number of samples to be drawn from distribution
%	filtrsize	width of Gaussian filter to apply
%	bootstrapquantiles	quantiles to calculate, if doing bootstrap calculation of uncertainties
%
% OUTPUTS:
%	y	exceedance probabilities of x
%	randpaths	path from sample 
%	bootstraps	bootstrap quantiles of y
%
% SEE ALSO: PathExceedanceProbability
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 13 August 2009
	
	defval('goodpts',[]);
	defval('Nsamples',100);
	defval('filtrsize',0);
	defval('bootstrapquantiles',[0.025 0.16 0.5 0.84 0.975]);
	defval('Nbootstraps',1000);
	defval('Nbootstrapsamples',1000)
		
	randpaths = [];
	disp('Sampling...')
	
	dogoodpts=0;
	randpaths=NaN * ones([size(f,1) Nsamples size(f,2)]);
	if length(goodpts)>0
		dogoodpts = prod(1*(size(goodpts) == size(f)));
		if ~dogoodpts
			disp('WARNING: goodpts provided but does is not the same size as f!');
		end
	end
	
	if dogoodpts
		sub1=find(sum(goodpts)>0);
	else
		sub1=1:size(f,2);
		sub2=1:size(f,1);
		goodpts = ones(size(f));
	end
	tic
	for i=sub1
		if dogoodpts
			sub2 = find(goodpts(:,i));
		end
		subbad = setdiff(1:size(f,1),sub2);
		if mod(i,50)==0
				disp(['    ' num2str(i) '/' num2str(size(f,2)) '  [' num2str(toc) ' seconds]' ]);
		end
		if ndims(V)<3
			randpaths(:,:,i) = mvnrnd(f(:,i)',diag(V(:,i).^2),Nsamples)';

		else
			try
				randpaths(:,:,i) = mvnrnd(f(:,i)',V(:,:,i),Nsamples)';
			catch
				V(:,:,i) = .5 * (V(:,:,i)+V(:,:,i)');
				lambda = 0;
				[m,n] = size(V(:,:,i));
				[U2,S2,V2] = svd(V(:,:,i),0);
		   		s2 = diag(S2); 		
		   		if sum(log(s2))<log(eps)
%					disp('    Determinant of covariance = 0!');
%					disp('    Searching for coefficient of correction...');
					guess = sqrt(exp(log(eps)/(1+m)));
					lambda = fzero( @(x) sum(log(s2+(x^2)))-log(1e-10),guess)^2;
					lambda = max(lambda);
%					disp(['    Adding ' num2str(lambda) '*I']); 
				end
				try
					randpaths(:,:,i) = mvnrnd(f(:,i)',V(:,:,i)+eye([m n])*lambda,Nsamples)';
				catch
					subbad=1:size(f,1);
					disp(['   Tossing out ' num2str(i)]);
				end
			end
		end
		randpaths(subbad,:,i) = NaN;
	end

	randpaths = reshape(randpaths,size(randpaths,1),[],1);
	maxes = max(randpaths);

	if filtrsize>0
		for i=1:size(f,1)
			filtr(i,:) = pdf('norm',1:size(f,1),i,filtrsize);
			filtr(i,:) = filtr(i,:)/sum(filtr(i,:));
		end
		randpaths = filtr*randpaths;
	end

	y=PathExceedanceProbability(x,randpaths);
%	maxes=max(randpaths);
%	maxcheck = repmat(x(:),[1 length(maxes)]);
%	checks = repmat(maxes,[length(x) 1]) <= maxcheck;
%	y = 1 - sum(checks,2)/length(maxes);
	
	if nargout>2
		bs = PathExceedanceProbability(x,randpaths,Nbootstraps,Nbootstrapsamples);
		bootstraps = quantile(bs,bootstrapquantiles,2);
	end
end