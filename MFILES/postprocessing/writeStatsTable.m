function writeStatsTable(fname,ss)

% writeStatsTable(fname,ss)
%
% INPUT:
%	fname	name of file to write
%	ss		data structure
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 31 May 2009

	fields = fieldnames(ss);
	headers = '';
	for j=1:length(ss.(fields{1}))
		datarows{j}='';
	end
	
	for i=1:length(fields)
		headers = [headers fields{i} '\t'];
		for j=1:length(ss.(fields{i}))
			if iscell(ss.(fields{i}))
				datarows{j} = [datarows{j} ss.(fields{i}){j} '\t'];
			else
				datarows{j} = [datarows{j} num2str(ss.(fields{i})(j)) '\t'];
			end
		end
	end

	fid=fopen(fname,'w+');
	fprintf(fid,[headers '\n']);
	for j=1:length(datarows)
		fprintf(fid,[datarows{j} '\n']);
	end
	fclose(fid);
end