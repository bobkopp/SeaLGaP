function [M,xpts,ypts]=histogram2D(x,y,Nptsx,Nptsy,weights,doplot)

% [M,xpts,ypts]=histogram2D(x,y,Nptsx,Nptsy,weights,doplot)
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 1 June 2009

defval('Nptsx',200);
defval('Nptsy',200);
defval('weights',ones(size(y)));
defval('doplot',2);

if length(Nptsx)>1
	xpts = Nptsx;
	Nptsx = length(xpts);
else
	xpts = linspace(min(x),max(x),Nptsx);
end

if length(Nptsy)>1
	ypts = Nptsy;
	Nptsy = length(ypts);
else
	ypts = linspace(min(y),max(y),Nptsy);
end

M = zeros(Nptsx-1,Nptsy-1);
[LBX,LBY] = meshgrid(xpts(1:end-1),ypts(1:end-1));
[UBX,UBY] = meshgrid(xpts(2:end),ypts(2:end));

for i=1:length(y)
	M = M + weights(i) * ( (x(i)>=LBX) .* (x(i)<UBX) .* (y(i)>=LBY) .* (y(i)<LBY));
end

if doplot>0
	gmap = gray;
	gmap = gmap(end:-1:1,:);
	
	imagesc(.5*(xpts(1:end-1)+xpts(2:end)),.5*(ypts(1:end-1)+ypts(2:end)),M') ; set(gca,'YDir','Normal') ; colormap(gmap);
		
end

