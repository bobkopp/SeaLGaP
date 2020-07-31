function y=SealevelDBImport(filenam,varargin)
	minage = -Inf;
	maxage = Inf;
	keepall = 0;

	if length(varargin)>0
		minage=varargin{1}(1);
		maxage=varargin{1}(2);
	end
	if length(varargin) > 1
		for i=2:length(varargin)
			if strcmpi(varargin{i},'keepall')
				keepall = 1;
			end
		end
	end

% y = SealevelDBImport(filename)

	columns={ 'ID',
	 'SiteID',
	 'Usable',
	 'Region',
	 'Site',
	 'SequenceOrder',
	 'StratigraphicAfterTime',
	 'StratigraphicAfterTimeError_SD',
	 'ResearchGroup',
	 'BibRef',
	 'LatDD',
	 'LongDD',
	 'DatingTech',
	 'DatingMat',
	 'AgeNote',
	 'AgeComment',
	 'AgeMsmtsN',
	 'AgeReported',
	 'AgeReportedError2S',
	 'AgeMode',
	 'AgeMin',
	 'AgeMax',
	 'AgeSD',
	 'Altitude',
	 'AltitudeErr2S',
	 'Zdatum',
	 'MaterialCategory',
	 'Material',
	 'Species',
	 'LimitingDir',
	 'RefLevel',
	 'TidalStation',
	 'TidalMeanRange',
	 'TidalSpringRange',
	 'IndicMeaningwrtMTL',
	 'IndicMeaningErr_2s',
	 'PaleoRSL',
	 'VertError_M_2s',
	 'TectonicUpliftCalib',
	 'TectonicUpliftCalibCircularity',
	 'TectonicUpliftRateEst',
	 'TectonicUpliftRateEstMin',
	 'TectonicUpliftRateEstMax',
	 'TectonicUpliftRateEst_Mode',
	 'TectonicUpliftRateEst_SD',
	 'TectonicUpliftEst',
	 'TectonicUpliftEstMin',
	 'TectonicUpliftEstMax',
	 'PaleoSL',
	 'PaleoSLMin',
	 'PaleoSLMax',
 	 'PaleoSLSD'
	};

	fid=fopen(filenam);
	datarows=textscan(fid,'%s','delimiter','\n');
	fclose(fid);
	y.datrows=datarows{1};
	
	y.headerrow=y.datrows{1};
	y.datarows={y.datrows{2:length(y.datrows)}};
	for i=1:length(y.datarows)
		y.datarows{i} = regexprep(y.datarows{i},',,',', ,');
		y.datarows{i} = regexprep(y.datarows{i},',,',', ,');
		y.datarows{i} = regexprep(y.datarows{i},'#VALUE!','NaN');
	end
	
	k=1;
	done=0;
	rem=y.headerrow;
	while done==0
	   [y.headerdat{k},rem]=strtok(rem,',');
	   if isempty(y.headerdat{k})
	       done=1;
	   end
	   
	   matches = find(strcmpi(y.headerdat{k},columns));
	   if length(matches) > 0
	   	y.col.(columns{matches(1)}) = k;
	   end
   
   	   k=k+1;
	end
	
	rem = y.datarows;	
	for k=1:length(y.headerdat)
	   [y.datstrings{k},rem]=strtok(rem,',');
	end
	
	if isfield(y.col,'Site') y.site = y.datstrings{y.col.Site}; end;
	if isfield(y.col,'BibRef') y.cite = y.datstrings{y.col.BibRef}; end;
	
	if isfield(y.col,'Site') && isfield(y.col,'BibRef')
		prune = setdiff(1:length(y.site),strmatch(' ',y.site));
		y.site=y.site(prune); y.cite=y.cite(prune);
	end
		
	if isfield(y.col,'ID') y.dataids = str2num(str2mat(y.datstrings{y.col.ID})); end;
	if isfield(y.col,'Usable') y.usable = str2num(str2mat(y.datstrings{y.col.Usable})); end;
	if isfield(y.col,'SiteID') y.siteid = str2num(str2mat(y.datstrings{y.col.SiteID})); end;
	if isfield(y.col,'SequenceOrder') y.stratseqid = str2num(str2mat(y.datstrings{y.col.SequenceOrder})); end;
	if isfield(y.col,'AgeMode') y.age_mode= str2num(str2mat(y.datstrings{y.col.AgeMode})); end;
	if isfield(y.col,'AgeSD') y.age_sd= str2num(str2mat(y.datstrings{y.col.AgeSD})); end;
	if isfield(y.col,'LatDD') y.lat= str2num(str2mat(y.datstrings{y.col.LatDD})); end;
	if isfield(y.col,'LongDD')
		y.long= str2num(str2mat(y.datstrings{y.col.LongDD}));
		sub = find(isfinite(y.long));
		y.long(sub) = mod(y.long(sub),360);
	end;
	if isfield(y.col,'PaleoSL') y.sl_mode = str2num(str2mat(y.datstrings{y.col.PaleoSL})); end;
	if isfield(y.col,'PaleoSLSD') y.sl_sd=str2num(str2mat(y.datstrings{y.col.PaleoSLSD})); end;
	if isfield(y.col,'LimitingDir') y.sl_limiting = str2num(str2mat(y.datstrings{y.col.LimitingDir})); end;
	if isfield(y.col,'AgeMsmtsN') y.agemsmtsn = str2num(str2mat(y.datstrings{y.col.AgeMsmtsN})); end;
	if isfield(y.col,'TectonicUpliftCalibCircularity') y.tectonicupliftcalibcircularity = str2num(str2mat(y.datstrings{y.col.TectonicUpliftCalibCircularity})); end;
	
	if isfield(y.col,'PaleoRSL') y.rsl_mode = str2num(str2mat(y.datstrings{y.col.PaleoRSL})); end;
	if isfield(y.col,'VertError_M_2s') y.rsl_sd = str2num(str2mat(y.datstrings{y.col.VertError_M_2s}))/2; end;
	
	if isfield(y.col,'TectonicUpliftRateEst_Mode') y.upliftrate_mode =  str2num(str2mat(y.datstrings{y.col.TectonicUpliftRateEst_Mode})); end;
	if isfield(y.col,'TectonicUpliftRateEst_SD')y.upliftrate_sd =  str2num(str2mat(y.datstrings{y.col.TectonicUpliftRateEst_SD})); end;
	
	if isfield(y.col,'StratigraphicAfterTime') y.strataftertime_mode = str2num(str2mat(y.datstrings{y.col.StratigraphicAfterTime})); end
	if isfield(y.col,'StratigraphicAfterTimeError_SD') y.strataftertime_sd = str2num(str2mat(y.datstrings{y.col.StratigraphicAfterTimeError_SD})); end;
	
	if isfield(y.col,'ResearchGroup') y.ResearchGroup=str2num(str2mat(y.datstrings{y.col.ResearchGroup})); end
	if isfield(y.col,'DatingTech') y.DatingTech = y.datstrings{y.col.DatingTech}; end
	if isfield(y.col,'MaterialCategory') y.MaterialCategory = y.datstrings{y.col.MaterialCategory}; end

	if ~keepall
		subset=find((y.usable).*((y.age_mode+2*y.age_sd)>minage).*((y.age_mode-2*y.age_sd)<maxage));
	else
		subset = 1:length(y.dataids);
	end
	
	N = length(y.dataids);
	fields = fieldnames(y);
	for i=1:length(fields)
		if length(y.(fields{i})) == N
			y.(fields{i}) = y.(fields{i})(subset);
		end
	end
	
	
	if isfield(y,'DatingTech')
		clear t;
		t.U_Th=[];
		t.stratigraphic = [];
		t.d18O = [];
		t.palynology = [];
		t.AAR = [];
		t.OSL = [];
		t.IRSL = [];
		t.C14 = [];
		t.biostratigraphy = [];
		t.archaeology = [];
		t.soil = [];
		for i=1:length(y.DatingTech)
			q=y.DatingTech{i};
			if regexp(q,'U-Th')
				t.U_Th = [t.U_Th i];
			end
			if regexp(q,'stratigraphic')
				t.stratigraphic = [t.stratigraphic i];
			end
			if regexp(q,'oxygen isotopes')
				t.d18O = [t.d18O i];
			end
			if regexp(q,'palynology')
				t.palynology = [t.palynology i];
			end
			if regexp(q,'AAR')
				t.AAR = [t.AAR i];
			end
			if regexp(q,'OSL')
				t.OSL = [t.OSL i];
			end
			if regexp(q,'IRSL')
				t.IRSL = [t.IRSL i];
			end
			if regexp(q,'C14')
				t.C14 = [t.C14 i];
			end
			if regexp(q,'biostrat')
				t.biostatigraphy = [t.biostratigraphy i];
			end
			if regexp(q,'archaeol')
				t.archaeology = [t.archaeology i];
			end
			if regexp(q,'soil devel')
				t.soil = [t.soil i];
			end
		end
		y.Index.DatingTech = t;
	end
	
	if isfield(y,'MaterialCategory')
		clear t;
		t.facies = [];
		t.constructional = [];
		t.erosional = [];
		t.corals = [];
		t.isotopes = [];
		for i=1:length(y.MaterialCategory)
			q=y.MaterialCategory{i};
			if regexp(q,'facies')
				t.facies = [t.facies i];
			end
			if regexp(q,'constructional')
				t.constructional = [t.constructional i];
			end
			if regexp(q,'erosional')
				t.erosional = [t.erosional i];
			end
			if regexp(q,'corals')
				t.corals = [t.corals i];
			end
			if regexp(q,'isotopes')
				t.isotopes = [t.isotopes i];
			end
		end
		y.Index.MaterialCategory = t;
	end	
end