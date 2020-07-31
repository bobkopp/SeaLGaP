function [f,V,logp,alfa,errorflags,invtraincv] = GaussianProcessRegression(x0,y0,x,traincv,cvfunc,varargin)

% [f,V,logp,alfa,errorflags,invtraincv] = GaussianProcessRegression(x0,y0,x,traincv,cvfunc,[testcv2])
%
% Modeled on Rasmussen & Williams (2006)
%
% INPUT
%
% 	x0			x at training points
%	y0			y at training points
%	x			x at which to approximate
%	traincv		covariance matrix among covariance points
%	cvfunc		EITHER a handle for a function that takes x and x0 as inputs
%				OR the covariance matrix between x and x0
%	testcv2		IF cvfunc is passed as a matrix, then the covariance matrix among x
% 
% Last updated by Bob Kopp rkopp-at-princeton.edu, 17 May 2009

%%%%%

	errorflags=0;

	if nargin == 5
		testcv = feval(cvfunc,x0,x);
		testcv2 = feval(cvfunc,x,x);
	else
		testcv = cvfunc;
		testcv2 = varargin{1};
	end
	
%	diag(traincv)
%	try
%		L=chol(traincv,'lower');
%		invtraincv = pinv(L') * pinv(L);
%		alfa=L'\(L\y0);
%		doChol=1;
%	catch
%		disp('Not positive definite!')
		errorflags=-1;
		doChol = 0;
		   [m,n] = size(traincv);
		   [U,S,V] = svd(traincv,0);
		   s = diag(S);
		   tol = max(m,n) * eps(max(s));
		   r = sum(s > tol);
		   invtraincv = V(:,1:r)*diag(s(1:r).^-1)*U(:,1:r)';      
	             alfa = invtraincv * y0;
%	            disp('SVD')
%	end
	try
	f=testcv'*alfa;
	catch
	keyboard
	end

	logpterms(1) = -.5*abs(y0'*alfa);
	logpterms(3) = -.5*length(y0)*log(2*pi);
	% logpterms(2) = -.5 * log(det(traincv));
	if doChol
		v = L\testcv;
		V=testcv2-v'*v;
		logpterms(2) = - sum(log(diag(L)));
	else
		V=testcv2-testcv'*invtraincv*testcv;
		logpterms(2) = -.5 * sum(log((s(1:r))));
	end
		
	logp=sum(logpterms);
%	if ~isreal(logp) || logp > 0
%		disp(-[logpterms(:)' logp])
%	end

%	logp = -.5*y0'*alfa - .5*log(prod(s(1:r)))- .5*length(y0)*log(2*pi);

end
