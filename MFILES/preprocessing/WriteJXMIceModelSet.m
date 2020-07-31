function WriteJXMIceModelSet(model,dname,visc,basemodel,compress)

% WriteJXMIceModelSet(model,dname,,visc,basemodel,compress)
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 12 June 2009

	defval('compress',1);
	defval('visc','');
	defval('basmodel','');
	
	mkdir(dname);
	frmt = [repmat('%8.3f ',1,9) '%8.3f\n'];
	
	i=1;
	fid=fopen([dname '/' IceFileName(i)],'w+');
	fprintf(fid,'%6.2f\n',model.timesteps(i)*sign(model.timesteps(1)));
	fprintf(fid,frmt,full(model.ice(:,i)));
	fclose(fid);

	for i=2:length(model.timesteps)
		fid=fopen([dname '/' IceFileName(i)],'w+');
		fprintf(fid,'%6.2f\n',model.timesteps(i)*sign(model.timesteps(1)));
		fprintf(fid,frmt,full(model.ice(:,i)));
		fclose(fid);
	end
	
	if length(basemodel)>0
		fid=fopen([dname '/basemodel.dat'],'w+');
		fprintf(fid,basemodel);
		fclose(fid);
	end
	
	if length(visc)>0
		fid=fopen([dname '/visc.dat'],'w+');
		fprintf(fid,'%f\n',visc(1:3));
		fclose(fid);
	end
	
	if compress
		system(['zip -q -r ' dname '.zip ' dname]);
		system(['rm -fr ' dname]);	
	end
		
function ifn=IceFileName(ind)
	if ind==1
		ifn='ice.i';
	else
		ifn=['ice.' num2str(ind-1)];
	end