function [randpaths,discards]=RandomIceSamples(f,V,goodpts)

%
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 9 July 2009

defval('goodpts',[]);
defval('Nsamples',100);

discards=0;

dimsorig=size(f);
if (ndims(f)==3)
	% the last (third) dimension of f is iteration
	% and the penultimate (second) is ice sheet
	% so we need to reshape f into two dimensions
	
	f = reshape(f,[size(f,1)*size(f,2) size(f,3)]);
end	

randpaths=NaN * ones([size(f,1) Nsamples size(f,2)]);

dogoodpts=0;
if length(goodpts)>0
	dogoodpts = prod(1*(size(goodpts) == size(f)));
	if ~dogoodpts
		goodpts = repmat(goodpts,[dimsorig(2) 1]);
		dogoodpts = prod(1*(size(goodpts) == size(f)));
	end
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
				V(:,:,i) = (V(:,:,i)+V(:,:,i)')/2;
				randpaths(:,:,i) = mvnrnd(f(:,i)',V(:,:,i)+eye([m n])*lambda,Nsamples)';
			catch
				subbad=1:size(f,1);
				disp(['   Tossing out ' num2str(i)]); discards=discards+1;
			end
		end
	end
	randpaths(subbad,:,i) = NaN;
end

randpaths=reshape(randpaths,size(randpaths,1),[]);
keepers=find(sum(~isnan(randpaths),1));
randpaths=randpaths(:,keepers);
randpaths=reshape(randpaths,dimsorig(1),dimsorig(2),[]);