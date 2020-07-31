function y=GenerateIceModelSeeds(N)

% GenerateIceModelSeeds(N)
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 12 June 2009

	defval('litho_opts',[95 70 120]);
	defval('umantle_opts',[0.5 0.3 0.8 1.0]);
	defval('lmantle_opts',[10 2 3 5 8 20]);

	RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

	disp('Generating ice distributions...');
	[y.gsl,y.iced,y.steric] = SampleLR04BasedPrior(N);
	
	for i=1:N
		y.visco(i,:) = [randsample(litho_opts,1) randsample(umantle_opts,1) randsample(lmantle_opts,1)  ];
	end
	
end