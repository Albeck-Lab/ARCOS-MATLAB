classdef arcos_utils
    methods(Static)
		function out = formatr(filename,XCoord,YCoord,bin)
            array = zeros(size(bin,1)*size(bin,2),5); %preallocate array
            sz = size(bin,1);
            for col = 1:size(bin,2)
                for row = 1:size(bin,1)
                    if col >1
                        ind = row+sz*(col-1); 
                    else
                        ind = row;
                    end
                    array(ind,1) = col; %Time index
                    if isnan(XCoord(row,col)) %If NaN set to 0
                        array(ind,2) = 0;
                    else
                        array(ind,2) = XCoord(row,col); %XCoord
                    end
                    if isnan(YCoord(row,col)) %If NaN set to 0
                        array(ind,3) = 0;
                    else
                        array(ind,3) = YCoord(row,col); %YCoord
                    end
                    array(ind,4) = bin(row,col); %Inactive = 0, active = 1
                    array(ind,5) = row; %Cell ID
                end
            end
            table = array2table(array); %Convert array to table
            table.Properties.VariableNames(1:5) = ["t","x","y","m","id"]; %Set table field names
			out = table;
            writetable(table,append(filename,'.csv')); %Write to file
        end 
        function [eps,minpts] = prep_dbscan(XCoord, YCoord)
            minpts = ndims(XCoord)*2;
            [~,d] = knnsearch([XCoord,YCoord], [XCoord, YCoord],'K', minpts+1); %k-nearest neighbors search
            d = max(d,[],2); %Biased toward greater distances as opposed to average of k-nearest
            max_d = sort(d);
            scaled = max_d * length(max_d)/max(max_d); %Scale the data
            smoothed = smoothdata(scaled,'gaussian'); %Smooth it
            slopes = gradient(smoothed); %Take first derivative
            [~,ix]=min(abs(slopes-1)); %Get the index of avg_d for ideal eps (where slope of line tangent to that point is 1);
            eps = max_d(ix);
        end 
        function [eps,minpts] = prep_dbscan2(XCoord, YCoord, time)
            minpts = ndims(XCoord)*2; %Minpts defined by dimensionality of data
            vals = zeros(1,length(time)); %Preallocation of data
            for it = 1:numel(time)
                t = time(it);
                [~,d] = knnsearch([XCoord(:,t),YCoord(:,t)], [XCoord(:,t), YCoord(:,t)],'K', minpts+1); %k-nearest neighbors search
                d = max(d,[],2); %Biased toward greater distances as opposed to average of k-nearest
                max_d = sort(d);
                scaled = max_d * length(max_d)/max(max_d); %Scale the data
                smoothed = smoothdata(scaled,'gaussian'); %Smooth it
                slopes = gradient(smoothed); %Take first derivative
                [~,ix]=min(abs(slopes-1)); %Get the index of avg_d for ideal eps (where slope of line tangent to that point is 1);
                eps = max_d(ix); %Set epsilon to the value of max_d at index ix
                vals(it) = eps; %Store eps in vals
            end
            real = ~isnan(vals); %Logical map of values ~= NaN
            eps = mean(vals(real)); %Mean of non-NaN epsilon values
        end
        function [eps,minpts] = prep_dbscan3(XCoord,YCoord,bin)
            xy = [XCoord,YCoord];
            P = 0;
            n = 11;
            [~,d] = knnsearch(xy,xy,'K', n);
            d_real = sort(d(~isnan(d(:,n)),n));
            eps = median(d_real); %Median distance to 11th neighbor of all points
            p = size(XCoord(bin),1)/size(XCoord,1); %Probability of cell being active
            for k = n:-1:1
                newP = nchoosek(n,k) .* p.^k .* (1-p).^(n-k);
                P = P + newP;
                if P > 0.01
                    %k = k-1;
                    break;
                end
            end
            if k < 3
                minpts = 3;
            else
                minpts = k;
            end
		end
		function [bin,threshold] = binarize(raw_data,xy,ch,perc)
			bin = cell(1,numel(xy));
			MegaDataHolder = [];
			%% Loop through wells, append channel data to MegaDataHolder
			for i = 1:numel(xy)
				well = xy(i);
				channel = raw_data{well}.data.(ch{1});
				MegaDataHolder = vertcat(MegaDataHolder, channel);
			end
			%% Get mean percentile of all well's channel data
			threshold = mean(prctile(MegaDataHolder,perc));
			%% Loop through wells again and binarize
			for i = 1:numel(xy)
				well = xy(i);
				channel = raw_data{well}.data.(ch{1});
				
				bin{well} = channel>threshold;
			end
		end 
		function out = pulse2bin(pulseAnalysisData,raw_data,channel)
    		%% Check for XYs with data and get their indicies
    		xy = 1:numel(pulseAnalysisData);
    		goodxys = ~cellfun(@isempty,pulseAnalysisData);% check to see if the input xys are good
    		xy = xy(goodxys);
    		
    		%% Loop through XYs
    		numXYs = numel(xy);
    		
    		out = {};
    		for iwell = 1:numXYs
        		well = xy(iwell);
        		%% Store pulse data and coords
        		pulseData = pulseAnalysisData{well}.data.(channel);
    		
        		% Make a data container
        		dataCont = false(size(raw_data{well}.data.XCoord,1),size(raw_data{well}.data.XCoord,2));
        		
                % Assign tracked points as zeros FIX - data is NaN unless real?
        		%dataCont(~isnan(raw_data{well}.data.XCoord)) = 0;
    		
        		%% Loop through pulse data
        		for i = 1:size(pulseData,1)
            		if ~isempty(pulseData(i).pkpos)
                		for ii = 1:size(pulseData(i).pkpos,1)
                    		%% Get start position and end (from dur of each pulse
                    		tStart = pulseData(i).pkpos(ii);
                    		tEnd = tStart + pulseData(i).dur(ii);
                    		%% Fill in the these points with 1s!
                    		dataCont(i,tStart:tEnd) = 1;
                		end
            		end
        		end
        		out{well} = dataCont;
    		end %XY loop
		end %convertPulseToBin
        function bin_synth = gensynth(raw_data,xy,varargin)
            %% Default parameters
            p.numspreads = 7; %How many spreads occur per cycle
            p.freq = 10; %Spread frequency
            p.bin = []; %Optional binarized data - useful for adding "noise"
            p.dist = 1; %Distance the spread grows per frame (epsilon from prep_dbscan is a useful metric)
            p.lifetime = 3; %How long a cell within a spread remains active before switching off
            p.seed = []; %Optional seed value for random number generator
            p.maxsize = 2.8; %Max spread size
            %% Setup
            nin = length(varargin);     %Check for even number of add'l inputs
            if rem(nin,2) ~= 0; warning('Additional inputs must be provided as option-value pairs'); end%Splits pairs to a structure
            for s = 1:2:nin; p.(lower(varargin{s})) = varargin{s+1}; end
            %% Random Number Generator Seed
			if ~isempty(p.seed)
				seed = str2double(evalc('disp(p.seed)'));
				rng(seed)
			end%Use user-provided seed
            %% Pre-allocate output
            bin_synth = cell(1,length(xy));
            %% Loop through XYs (wells)
            for iwell = 1:numel(xy)
                well = xy(iwell);
                XCoord = raw_data{well}.data.XCoord; % Assign X and Y Coords from data
                YCoord = raw_data{well}.data.YCoord;
                binSynth = boolean(zeros(size(XCoord,1), size(XCoord,2)));  % Pre-allocate logical array
                %% Generate random start indices
				for r = 1:round(size(binSynth,2)/p.freq+1)
                    col = r*p.freq-(p.freq-1);
					for q = 1:p.numspreads %Generate random indices for starting pts
                        istartPts(q,r) = randi(size(XCoord,1)); %#ok<AGROW>
                        binSynth(istartPts(q,r),col) = 1;
					end
				end
				%% Starting parameters
                cnt = 1;
                sz = 1;
                %% Loop through time series
                for time = 1:size(XCoord,2)
					movingframe = (time-p.freq:time); %A sliding frame of time indices from the current time back to time-p.lifetime
					movingframe = movingframe(movingframe>0); %Filter by values > 0
					if ismember(time,(1:p.freq:size(XCoord,2)))
                        startPts = [XCoord(istartPts(:,cnt),time),YCoord(istartPts(:,cnt),time)]; %Get coords for start pts
                        cnt = cnt+1;
                        sz = 1;
					else
						%% Set points within radius active
						pts = [XCoord(:,time),YCoord(:,time)]; %Get all xy coords at current time
						[idx,d] = rangesearch(pts,startPts,p.dist*sz); %Get points within dist*time of start pts
						sz = sz+1;
						for i = 1:size(idx)
							binSynth(idx{i}(d{i}<= p.dist*p.maxsize),time) = 1; %Set points within dist*time of start pts to active
							for row = 1:size(binSynth,1)
								if sum(binSynth(row,movingframe))>=p.lifetime+1
									binSynth(row,time)=0; %Set cells to inactive if they've been active > lifetime
								end
							end
						end
						%% Set points inactive if they've been alive for p.lifetime (Donut)
						
					end
                end
                %% Layer real data on top of synthetic
                if ~isempty(p.bin) %Use user-provided binary data
                    bin_synth{well} = logical(binSynth+p.bin{well});
                else
                    bin_synth{well} = logical(binSynth);
                end
            end %end well loop
		end 
		function out = reformat(cdata)
			%% Load data into struct
			timerange = size(cdata,2);
			max_id = cdata(end).newmax;
			clusters = repmat(struct('cid',[],'data',struct('time',{},'XYCoord',{},'id',{},'numpts',{},'bounds',{},'inactive',{},'area',{},'compl',{},'rocarea',{},'roccount',{}),'t_start',[],'t_end',[],'dur',[], 'maxarea',[],'maxcount',[]),max_id,1); %set up struct of structs
			for time = 1:timerange
				eps = cdata(time).eps;
				minpts = cdata(time).minpts;
				dtp = cdata(time).tracked; %dtp = data at timepoint
				for cluster = 1:size(dtp,2)
					id = mode(dtp(cluster).id(:,2));
					clusters(id).cid = id;
					clusters(id).data(time).time = time;
					clusters(id).data(time).XYCoord = dtp(cluster).XYCoord;
					clusters(id).data(time).id = dtp(cluster).id;
					clusters(id).data(time).eps = eps;
					clusters(id).data(time).minpts = minpts;
				end
			end
			%% Loop through substructs and remove empty entries
			for i = 1:size(clusters,1)  
				map = false(1,size(clusters(i).data,2));
				for ii = 1:size(clusters(i).data,2)
					if isempty(clusters(i).data(ii).XYCoord) || isempty(clusters(i).data(ii).id)
						map(ii) = 1;
					end
				end
				clusters(i).data(map) = [];
			end
			out = clusters;
		end
		function out = binarize_robust(raw_data,xy,ctrl,chan,type,perc)
			MegaDataHolder = []; %Array of channel data, vertically concatenated
			if perc > 200 || perc < 0; error('Perc should be a percentage <200, ex 85'); end
			perc = perc/100;
			%% Loop over all the xys and combine them into a mega data
			% This loop defined by xys or control wells you specify
			if isempty(ctrl); xy2bin = xy; else; xy2bin = ctrl; end
			numXY2bin = numel(xy2bin);
			for iXY = 1:numXY2bin
    			tXY = xy2bin(iXY);
    			if ~isempty(raw_data{tXY}) %check the XY isn't empty
    			if isfield(raw_data{tXY}, 'data') %check the XY has data
    			if isfield(raw_data{tXY}.data, chan) %check if that XY has the channel in it
        			MegaDataHolder = vertcat(MegaDataHolder, raw_data{tXY}.data.(chan)); %#ok<AGROW> %concat to megadataholder
    			end %if the chan exists
    			end %if data exists
    			end %if raw_data is not empty
			end %end normdata collection loop
    		%% Calculate the method of normalization requested
			switch type
    			case 'znorm' %Z-Score normalize (subtract mean of all data, divide by std dev)
        			MegaDataHolder = MegaDataHolder - min(MegaDataHolder,[],2); %subtract the min first
        			Subtract = mean(MegaDataHolder, 'omitnan'); %get the mean
        			DivideBy = std(MegaDataHolder, 0, 'all', 'omitnan'); %get std dev
    			case {'robust','selfrobust'} %Robust Scalar (scales to median and quantiles)
        			Subtract = median(MegaDataHolder, 'omitnan'); %get the median
        			Quantilz = quantile(MegaDataHolder, (0:0.25:1),'all');
        			Q25 = Quantilz(2); Q75 = Quantilz(4); %get the 25th and 75th quantiles
        			DivideBy = Q75 - Q25; %get the interquantile range
    			case {'specialrobust'} %Robust Scalar (scales to median and quantiles)
        			Subtract = median(MegaDataHolder, 2,'omitnan'); %get the median for each cell, to self subtract
        			Quantilz = quantile(MegaDataHolder, (0:0.25:1),2);
        			Q25 = Quantilz(:,2); Q75 = Quantilz(:,4); %get the 25th and 75th quantiles
        			DivideBy = median(Q75 - Q25,'omitnan'); %get the median interquantile range to divide all data by
    			otherwise %use the default of zscore, but warn them about it
        			%fprintf('Normalization method %s given as input not recognized, using Z-Score normalization instead. \n', p.normalize)
        			type = 'znorm';
        			Subtract = mean(MegaDataHolder, 'omitnan'); %get the mean
        			DivideBy = std(MegaDataHolder, 0, 'all', 'omitnan'); %get std dev
			end       
            %% Loop over all the xys again to normalize them
			numXYs = numel(xy);
			out = cell(1,numXYs);
            for iXY = 1:numXYs
                tXY = xy(iXY);
                if ~isempty(raw_data{tXY}) %check the XY isn't empty
                if isfield(raw_data{tXY}, 'data') %check the XY has data
                if isfield(raw_data{tXY}.data, chan) %check if that XY has the channel in it
                    switch type
                        case 'specialrobust'
                            Subtract = median(raw_data{tXY}.data.(chan), 2,'omitnan'); %get the median for each cell, to self subtract
                            out{tXY} = (((raw_data{tXY}.data.(chan)) - Subtract) ./ DivideBy) > DivideBy*perc; %make normalized data
                        otherwise
                            out{tXY} = (((raw_data{tXY}.data.(chan)) - Subtract) ./ DivideBy) > DivideBy*perc; %make normalized data
                    end %end switch
					%mean or median of DivideBy value (divide it by 2),
					% use that number to threshold.
                end %if the chan exists
                end %if data exists
                end %if raw_data is not empty
            end %end second iXY loop
		end
    end
end