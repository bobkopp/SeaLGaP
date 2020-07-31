function model=ConstructJXMIceModel(model0,t,iced,esl0)

% model=ConstructJXMIceModel(model0,t,iced,esl0)
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 12 June 2009

defval('viewoutput',0);

defval('t',150:-1:90);

defval('esl0',[0.18 7.3 0.02 5.0 58.0 0 0]);
icevol0 = sum(esl0); esl0=esl0(:);

% constants
rad = 6.371e6; % radius of the Earth
RHOE=5511.57; RHOW=1000.0; RHOI=920.0; %densities
areaOcean = 361419000 * 1e6;

esladj=model0.sheetesl(:,end)-esl0;
iced2=iced+repmat(esladj,[1 size(iced,2) size(iced,3)]);
icevol0new=icevol0+sum(esladj);

sub=find(iced2<0);
[i,j,k]=ind2sub(size(iced2),sub);
% put all the diff on to EAIS
negice=iced2(sub);
eaisshift=5-i;
for i=1:length(sub)
	iced2(sub(i))=iced2(sub(i))-negice(i);
	iced2(sub(i)+eaisshift(i))=iced2(sub(i)+eaisshift(i))+negice(i);
end

terminalconfig = model0.sheetesl(:,1);
deltaterminalconfig = repmat(terminalconfig,[1 1 size(iced2,3)]) - iced2(:,end,:);
transition = repmat(.1:.1:1,[size(deltaterminalconfig,1) 1 size(deltaterminalconfig,3)]) .* repmat(deltaterminalconfig,[1 10 1]) + repmat(iced2(:,end,:),[1 10 1]);
%icevoltrans = squeeze(sum(transition,1)-icevol0new)';

iced2(:,end+1:end+10,:) = transition;
t=[t (t(end)-[1:10])];
% find the closest ice configuration
clear iceconfscale iceconfind;
for i=1:size(iced2,1)
	icelist=squeeze(iced2(i,:,:));
	[X,Y]=meshgrid(icelist(:),model0.sheetesl(i,:));
	[m,iceconfind(i,:)]=min(abs(X-Y)+100*(X>Y),[],1);
	iceconfscale(i,:) = icelist(:)./model0.sheetesl(i,iceconfind(i,:))';
	clear X Y;
end
iceconfscale(find(isnan(iceconfscale)))=0;
iceconfind=reshape(iceconfind,size(iced2));
iceconfscale=reshape(iceconfscale,size(iced2));

icemasses=SeparateIceMasses(model0.lat,model0.long);
C = zeros(length(model0.long),size(iced,2));
k=1;
for j=1:size(iced2,2)
for i=1:size(iced2,1)
	sub=find(icemasses==i);
	C(sub,j)=model0.ice(sub,iceconfind(i,j,k))*iceconfscale(i,j,k);
end
end


model.timesteps=-t;
model.ice=C;
model.long = model0.long;
model.lat = model0.lat;
model.volfactor = model0.volfactor;

vol = repmat(model.volfactor,[1 size(model.ice,2)]) .* model.ice;
model.esl = (RHOI/RHOW) * vol / areaOcean;

%model.sheets=[1:7]';
%model.sheetesl=iced;

icemasses = SeparateIceMasses(model.lat,model.long);
model.sheets = unique(icemasses);
for i=1:length(model.sheets)
	sub = find(icemasses == model.sheets(i));
	model.sheetesl(i,:) = (sum(model.esl(sub,:)));
end
model.referenceesl = model0.sheetesl(:,end);
model.sheetesl = full(model.sheetesl);
difrence=sum(sum(abs(model.sheetesl-iced2)));
if difrence>1e-3
	warning(['model.sheetesl differs from target ice distribution -- sum(abs(model.sheetesl-iced2)) = ' num2str(difrence)]);
end

% look at it

if viewoutput
	disp('Displaying output...');
	long=unique(model.long);
	lat=unique(model.lat); lat=lat(end:-1:1);
	for j=1:size(model.ice,2)
	clf
	imagesc(long,lat,reshape(model.ice(:,j),[length(lat) length(long)])) ; colorbar ; caxis([0 5000]); set(gca,'ydir','norm');
	[axlim,hcont]=plotcont; set(hcont','LineW',2,'Color',[.5 .5 .5]);
	axis image; 
	titlelabel{1}=[num2str(t(j)) ' ka'];
	titlelabel{2}=['\DeltaGlobal Ice Volume: ' num2str(sum(model.referenceesl)-sum(model.sheetesl(:,j))) ' m'];
	titlelabel{3}=['Laurentide: ' num2str(model.sheetesl(1,j)) ' m'];
	title(titlelabel);
	pause
	end
end