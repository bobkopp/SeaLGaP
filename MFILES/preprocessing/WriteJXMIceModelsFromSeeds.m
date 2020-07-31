function WriteJXMIceModelsFromSeeds(seeds,N,modeldirs,targetdir)

% WriteJXMIceModelsFromSeeds(seeds,N,modeldirs,targetdir)
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 12 June 2009

defval('modeldirs',{'ice-models/ice5g','ice-models/bassett'});
defval('targetdir','new-ice-models2');
defval('modelages',[-26 0 ; -28 0]);
defval('t',150:-1:90);
defval('esl0',[0.18 7.3 0.02 5.0 58.0 0 0]);
defval('dnameprefix',['kopp2009d-']);

startdir = pwd;

for i=1:length(modeldirs)
	cd(startdir);
	cd(modeldirs{i});
	basemodel{i}=ImportJXMIceModelSet([],modelages(i,1),modelages(i,2));
end

cd(startdir);
cd(targetdir);

for i=N
	curbase=mod((i-1),length(modeldirs))+1;
	disp(['Generating model ' num2str(i) ' based on base model ' modeldirs{curbase} '...']);
	dname=[dnameprefix sprintf('%04i',i)];
	model = ConstructJXMIceModel(basemodel{curbase},t,seeds.iced(:,:,i),esl0);
	WriteJXMIceModelSet(model,dname,seeds.visco(i,:),modeldirs{curbase},1);
end
