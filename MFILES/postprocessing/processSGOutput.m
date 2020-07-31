function [fKeep,f_sdKeep,gKeep,f,sd,datids,fIces,sdIces,fIcesDeriv,sdIcesDeriv,Vsheets,VIcesDeriv,data] = processSGOutput(datfiles,t,Nburnin,varargin)

% [fKeep,f_sdKeep,gKeep,f,sd,datids,fIces,sdIces,fIcesDeriv,sdIcesDeriv,Vsheets,VIcesDeriv,data] = processSGOutput(datfiles,t,Nburnin,varargin)
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 8 July 2009

	defval('Nburnin',50);

	runmean = 0;
	dt = diff(t);
	dt=dt(:);
%	diffop = (-[speye(length(t)-2) zeros(length(t)-2,2)] + [zeros(length(t)-2,2) speye(length(t)-2) ]) ./ repmat(-(dt(1:end-1)+dt(2:end)),[1 length(t)])  ;

	diffop = (-[speye(length(t)-1) zeros(length(t)-1,1)] + [zeros(length(t)-1,1) speye(length(t)-1) ]) ./ repmat(-(dt),[1 length(t)]);
	
	if length(varargin)>0
		for i=1:length(varargin)
			if strcmpi(varargin{i},'runmean')
				runmean = 1;
			end
		end
	end
	
	Nt = length(t);
	
	icelat = [repmat(-Inf,1,7) repmat(Inf,1,4)]';
	icelong = [1:7 1:3 Inf]';
	icelatR = repmat(icelat,[1 length(t)])';
	icelongR = repmat(icelong,[1 length(t)])';
	icetR = repmat(t(:)',[length(icelat) 1])';
	
	
	% we end up with a matrix with rows representing points in time
	% and columns representing different ice masses

	y = AssimilateRun(datfiles,Nburnin);
	sub = find(isfinite(sum(y.fKeep)));
	fKeep = y.fKeep(:,sub); gKeep = y.gKeep(:,sub); f_sdKeep = y.f_sdKeep(:,sub); data = y.data;
	if runmean
		sub = find(isfinite(sum(y.fKeep)));
		[fKeep,meanfsdkeep,meansdkeep,f_sdKeep] = weightedAverage(y.fKeep(:,sub),y.f_sdKeep(:,sub));
		gKeep = mean(y.gKeep(:,sub),2);
	end
	
	if isfield(y,'f_VKeep')
		if size(y.f_VKeep,1)>0
			f_VKeep=y.f_VKeep(:,:,sub);
		end
	end
	
	promptfreq = 20;
	
	N=size(fKeep,2);
	Nt = length(t);
	Nices = length(icelat);
	Nderiv = size(diffop,1);
	
	logpPriorIces=zeros(1,N);
	f = zeros(Nt,N);
	sd = zeros(Nt,N);
	fIces = zeros(Nt,Nices,N);
	sdIces = zeros(Nt,Nices,N);
	Vsheets = zeros(Nt,Nt,Nices,N);
	fIcesDeriv = zeros(Nderiv,Nices,N);
	VIcesDeriv = zeros(Nderiv,Nderiv,Nices,N);
	sdIcesDeriv = zeros(Nderiv,Nices,N);
	
	tic
	if nargout > 3
		for i=1:N
			if mod(i,promptfreq)==0
				disp(['    ' num2str(i) '/' num2str(size(fKeep,2)) '  [' num2str(toc) ' seconds]' ]);
			end
		
	%		[f(:,i),V,logpPrior] = GPRSpaceTimeRegressionExplicit(fKeep(:,i+Nburnin),data.lat(:),data.long(:),gKeep(:,i+Nburnin),f_sdKeep(:,i+Nburnin),[  Inf*ones(size(gsl1))],[ Inf*ones(size(gsl1))],[t]);
	%		sd(:,i)=sqrt(diag(V));
	
			if exist('f_VKeep')

				[fice,V,logpPriorIces(i)] = GPSpaceTimeRegression(fKeep(:,i),data.lat(:),data.long(:),gKeep(:,i),f_VKeep(:,:,i),icelatR(:),icelongR(:),icetR(:));
			else
				[fice,V,logpPriorIces(i)] = GPSpaceTimeRegression(fKeep(:,i),data.lat(:),data.long(:),gKeep(:,i),f_sdKeep(:,i),icelatR(:),icelongR(:),icetR(:));
			end

			fIces(:,:,i) = reshape(fice,(size(icelatR)));
			sdIces(:,:,i)=reshape(sqrt(diag(V)),size(icelatR));
	
			f(:,i) = fIces(:,end,i);
			sd(:,i) = sdIces(:,end,i);
	
			if nargout>8
		
				for j=1:length(icelong)
					chunk = [1:Nt] + (j-1)*Nt;
					Vsheets(:,:,j,i) = tril(V(chunk,chunk)) + tril(V(chunk,chunk),-1)';
					fIcesDeriv(:,j,i) = diffop*fIces(:,j,i);
					VIcesDeriv(:,:,j,i) = diffop*Vsheets(:,:,j,i)*diffop';
					sdIcesDeriv(:,j,i) = sqrt(diag(VIcesDeriv(:,:,j,i)));
				end
			end		
		end
	end
	
	datids = data.dataids;	
	

end