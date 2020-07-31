function [evalpts,M,quant,hquants]=cumdensity2d(x,y,sd,quants,Npts,weights,doplot,linecolor,linesyms,linewidths)

% [evalpts,M,quant]=cumdensity2d(x,y,sd,quants,Npts,weights,doplot,linecolor,linesyms,linewidths)
%
% INPUTS:
%	x		x coordinates of the rows of y (default: 1:size(y,2))
%	y		data (different x values in each row, different samples in each col.)
%	sd		standard deviations corresponding to y
%	quants	quantiles to return (default: [.025 .167 .5 .833 .975])
%	Npts	number of y values at which evaluate density function
%			OR y values at which to evaluate
%	weights	weights to assign each column of y (default: all 1)
%	doplot	produce plot (default = 0)
%					0: no plot
%					1: density plot
%					2: density plot with quantiles marked
%					3: quantiles marked only
%	linecolor	line color in plot type 2
%	linesyms	line symbols for quantile lines in plot type 2
%	linewidths	corresponding linewidths (default: [1 1 2 1 1]);
%
% OUTPUTS:
%	evalpts		y values at which M is evaluated
%	M			density at each evalpts and x value
%	quant		quantiles
%
% Computes a cumulative density function for data from across MC runs.
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 2 June 2009

	evalpts=[];
	
	defval('x',1:size(y,2));
	defval('quants',[.025 .167 .5 .833 .975]);
	defval('Npts',200);
	defval('linewidths',[1 1 2 1 1]);
	defval('weights',ones([1 size(y,2)]));
	defval('doplot',0);
	defval('linecolor','g');
	
	if length(Npts)>2
		evalpts = Npts;
		Npts = length(evalpts);
	end
	
	if length(evalpts)==0
		miny=min(y(:)-2*sd(:));
		maxy=max(y(:)+2*sd(:));
		evalpts=linspace(miny,maxy,Npts);
	end
	
	if ~exist('linesyms','var') && (doplot>=2)
		Nquants = length(quants);
		for i=1:Nquants
			linesyms{i} = '--';
		end
		if mod(Nquants,2)==1
			linesyms{ceil(Nquants/2)} = '-';
		else
			linesyms{(Nquants/2)} = '-';
			linesyms{(Nquants/2) + 1} = '-';
		end
		if Nquants>4
			linesyms{1} = ':';
			linesyms{Nquants} = ':';
		end
	end
	
	evalptsM = repmat(evalpts,[size(y,1) 1]);

	M = zeros(size(evalptsM));
	%disp('Calculating density...');
	for i=1:size(y,2)
		M = M + weights(i) * cdf('norm',evalptsM,repmat(y(:,i),[1 length(evalpts)]),repmat(sd(:,i),[1 length(evalpts)]));
	end
	M = M/sum(weights);
	
	for i=1:size(y,1)
		sub = find((M(i,:)>eps).*(M(i,:)<(1-eps)));
		[u,ind]=unique(M(i,sub),'first');
		sub=sub(ind);
		if length(sub)<2
			quant(i,:)=NaN*ones(size(quants));
		else
			quant(i,:) = interp1(M(i,sub),evalpts(sub),quants);
		end
	end
	
	if ((doplot>0)&&(doplot<3))
		gmap = gray;
		gmap = gmap(end:-1:1,:);
		imagesc(x,.5*(evalpts(1:end-1)+evalpts(2:end)),diff(M',1,1)) ; set(gca,'YDir','Normal') ; colormap(gmap);
	end
	
	if doplot >= 2
		hold on;
		for i=1:length(quants)
				hquants(i)=plot(x,quant(:,i),[linecolor linesyms{i}],'LineWidth',linewidths(i));
		end
	end

end