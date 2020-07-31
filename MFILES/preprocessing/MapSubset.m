function filtered=MapSubset(y,subset)

% filtered=MapSubset(y,subset)
%
% Returns a subset of the sea level database structure y, with indices specified by subset.
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 17 May 2009

	fields = fieldnames(y);
	N =length(y.sl_mode);
	for i=1:length(fields)
		if length(y.(fields{i})) == N
			y.(fields{i}) = y.(fields{i})(subset);
		end
	end
	filtered = y;


end