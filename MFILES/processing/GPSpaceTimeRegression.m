function [f,V,logp,errorflags,logpgrad,cvfunc,evalfunc,traincv,cvterms,invtraincv] = GPSpaceTimeRegression(y,lat,long,times,errors,varargin)

% [f,V,logp,errorflags,logpgrad,cvfunc,evalfunc,traincv,cvterms,invtraincv] = 
%	GPSpaceTimeRegression(y,lat,long,times,errors,[lat1],[long1],[times1],[cvfile],[tprsc])
%
% Does spatiotemporal Gaussian process regression using covariance lookup table.
%
% See also: LookupCovariance
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 24 July 2009

	doself = 1;
	autocorrect = 1;
	cvfile = 'cvtable';
	tprsc=[];
		
	if length(varargin)>2
		lat1=varargin{1}; long1=varargin{2}; times1=varargin{3};
		doself = 0;
	end
	
	if length(varargin)>3
		if length(varargin{4})>0
			cvfile=varargin{4};
		end
		if length(varargin)>4
			if length(varargin{5})>0
				tprsc = varargin{5};
			end
		end
		
	elseif length(varargin)<=2
		if length(varargin)>0
			if length(varargin{1})>0
				cvfile=varargin{1};
			end
		end
		if length(varargin)==2
			if length(varargin{2})>0
				tprsc = varargin{2};
			end
		end
	end
	
	%%%%
	
	function [cv,mu1,mu2] = cvfun(x1,x2)
		[cv,mu1,mu2] = LookupCovariance(x1,x2,cvfile,[],tprsc);
		mu1 = mu1(:); mu2=mu2(:);
	end
	
	cvfunc = @(x1,x2) cvfun(x1,x2);
	
	function mu1 = positionmean(x1)
		[cv,mu1] = cvfun(x1,x1);
	end
	
	function [nf,nV,nlogp,nalfa,nerrorflags,ninvtraincv] = GPmeanadjusted(nx0,ny0,nx,ntraincv,ncvfunc,ntestcv2)
		% pass y0 with means subtracted, and return f with means added back in
		if nargin == 5
			[nf,nV,nlogp,nalfa,nerrorflags,ninvtraincv] = GaussianProcessRegression(nx0,ny0,nx,ntraincv,ncvfunc);
		elseif nargin > 5
			[nf,nV,nlogp,nalfa,nerrorflags,ninvtraincv] = GaussianProcessRegression(nx0,ny0,nx,ntraincv,ncvfunc,ntestcv2);
		end
		nf = nf+positionmean(nx);
	end
	
	%%%%
	
	[cvclean,meany0] = cvfunc([lat long times],[lat long times]);
	cvterms = cvclean;

	if ~doself
		meany1 = positionmean([lat1 long1 times1]);
	end

	if ((size(errors,1)==1)||(size(errors,2)==1))
		traincv = cvclean + diag(errors.^2);
	else
		traincv = cvclean + errors;
	end
	goodresult=0;
	function q=PositiveVarianceTest(addlerrors)
		addlerrors = abs(addlerrors);
		if doself
			[f2,V2] = GaussianProcessRegression([],y-meany0,[],traincv+eye(size(traincv))*addlerrors,cvclean,cvclean);
			f2 = f2+meany0;
		else
			[f2,V2] = GaussianProcessRegression([lat long times],y-meany0,[lat1 long1 times1],traincv+eye(size(traincv))*addlerrors,cvfunc);
			f2 = f2+meany1;
		end
		d = diag(V2);
		q = eps-sum(abs(d(find(d<eps))));

	end
	
	
	while ~goodresult
		evalfunc = @(x) GPmeanadjusted([lat long times],y-meany0,x,traincv,cvfunc);
		if doself
			[f,V,logp,alfa,errorflags,invtraincv] = GaussianProcessRegression([],y-meany0,[],traincv,cvclean,cvclean);
			f=f+meany0;
		else
			[f,V,logp,alfa,errorflags,invtraincv] = evalfunc([lat1 long1 times1]);
		end		
		d=diag(V);

%%%%%
%		disp([size(V) doself]);
%		clf; plot(d)
%		drawnow
%%%%%
		
		if eps-sum(abs(d(find(d<eps)))) <=0
			disp('Non-positive variances!')
			if autocorrect
				disp(['    Searching for coefficient of correction...']);
				lambda = (fzero(@(x) PositiveVarianceTest(x^2),0,optimset('TolX',1e-4)))^2;
				disp(['    Correcting by adding error of ' num2str(lambda)]);
				traincv = traincv + eye(size(traincv))*lambda;
				autocorrect = 0;
			end
		elseif sum(~isfinite(f))
			disp('Non-finite results!');
			goodresult = 1;

		else
			goodresult = 1;
		end
	end

	logpgrad= NaN;
	
end

