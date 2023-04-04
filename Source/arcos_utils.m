%% ARCOS Utils (Utilities)
% Class with various static utility functions
%
%% formatR
% Method for formatting data such that it can be read by the R version of
% ARCOS. Outputs a csv (comma separated value) file.
%
% *Inputs*
%
% * *filename* - |String|, |Char| - Desired output csv file name. Must
% contain legal characters only. Do not include the file extension in the
% name.
% * *XCoord* - |Array| - An array of doubles containing the X coordinates
% for tracked cells. Must be organized such that rows are individual
% tracked cells and columns are timepoints. 
% * *YCoord* - |Array| - An array of doubles containing the Y coordinates
% for tracked cells. Must be organized such that rows are individual
% tracked cells and columns are timepoints. 
% * *bin* - |Array| - A logical array indicating active and inactive cells.
% Must be organized such that rows are individual
% tracked cells and columns are timepoints.
%
% *Output* 
% 
% * *out* - |Table| - A table formatted to be read by the R version of
% ARCOS. This will be saved to the workspace, but the csv file will be
% saved to the current working directoy unless a full path is given for the
% filename.
%% prep_dbscan
% Method for determining optimal epsilon and minpts values for DBSCAN.
%
% This method assumes that the optimal epsilon value is at the "elbow" or
% "knee" of the graph of k-nearest-neighbors distances for all datapoints.
%
% This method supplies epsilon and minpts values for a single microscope
% field at a single timepoint. 
%
% This method first asserts minpts = 2*dim where dim is the dimensionality of the
% data.
% 
% Next, it calculates the distances to the k-nearest neighbors where k =
% minpts, takes the maximum distances, scales them, applies gaussian
% smoothing and takes the first derivative. 
% 
% The ideal epsilon is the first point whose tangent line has a slope of 1.
% 
% *Inputs*
%
% * *XCoord* - |Array| - An array of doubles containing the X coordinates
% for tracked cells. Must be organized such that rows are individual
% tracked cells and columns are timepoints.
% * *YCoord* - |Array| - An array of doubles containing the Y coordinates
% for tracked cells. Must be organized such that rows are individual
% tracked cells and columns are timepoints. 
%
% *Outputs*
%
% * *eps* - |Double| - The epsilon parameter of DBSCAN.
% * *minpts* - |Integer| - The minpts parameter of DBSCAN.
%% prep_dbscan2
% A method for determining optimal epsilon and minpts values for DBSCAN.
%
% This method increases minpts when data points are densely distributed or noisey and
% lowers minpts when they are sparsely distributed. 
%
% To ensure that DBSCAN runs as expected, minpts has a lower cap of 3.
% Minpts values < 3 
%
% Epsilon is calculated as the median of k-nearest neighbors distances
% where k = 11. 
%
% *Inputs*
%
% * *XCoord* - |Array| - An array of doubles containing the X coordinates
% for tracked cells. Must be organized such that rows are individual
% tracked cells and columns are timepoints.
% * *YCoord* - |Array| - An array of doubles containing the Y coordinates
% for tracked cells. Must be organized such that rows are individual
% tracked cells and columns are timepoints. 
% * *bin* - |Logical array| - A logical array representing the binarized
% channel data for the currently indexed well/xy (1 = active, 0 =
% inactive). 
% * _varargin_ = |Option-value pair| - Accepts optional inputs as
% option-value pairs.
%
% *Optional Inputs*
%
% * *n* - |Integer| - The nth nearest neighbor to use for knnsearch.
%
% *Ouputs*
%
% * *eps* - |Double| - The epsilon parameter of DBSCAN.
% * *minpts* - |Integer| - The minpts parameter of DBSCAN.
%% binarize
% Method for binarizing channel data using a simple percentile threshold.
% Channel data for all xys/wells are loaded into memory and the specified
% percentile of the data is used to threshold and binarize.
%
% *Inputs*
%
% * *raw_data* - |Cell| - Cell array containing processed data. Ex
% dataloc.d 
% * *xy* - |Array| - An array of integers representing indices of wells/xys
% to process. Supports discontinuous ranges - Ex (1:5, 10:20)
% * *ch* - |String|, |Char| - The name of the channel to process - Ex
% 'nEKAR'.
% * *perc* - |Double|, |Integer| - The percentile of the data by which it will be
% thresholded. 
%
% *Outputs*
%
% * *bin* - |Cell| - A cell array where each cell is binary data for the
% indexed well/xy. 
% * *threshold* - |Double| - The value by which the channel data are
% thresholded. Signals above this value are considered "on", signals below
% are "off". 
%% pulse2bin
% Method of binarizing pulse analysis data.
%
% *Inputs*
%
% * *pulseAnalysisData* - |Unknown| - Pulse analysis data. Ask Nick DeCuzzi for
% clarification.
% * *raw_data* - |Cell| - Cell array containing processed data. Ex
% dataloc.d 
% * *channel* - |String|, |Char| - The channel to binarize
%
% *Outputs*
%
% *out* - |Cell| - Cell array where each cell is binarized data for the
% indexed well/xy. 
%% gensynth
% Method to create artificial collective activity for process validation
% and demonstration.
%
% This method uses X and Y coordinates from real data and layers synthetic activity
% onto those points.
%
% *Inputs*
%
% * *raw_data* - |Cell| - Cell array containing processed data. Ex
% dataloc.d 
% * *xy* - |Array| - An array of integers representing indices of wells/xys
% to process. Supports discontinuous ranges - Ex (1:5, 10:20)
% * _varargin_ - |Option-value pair| - Supports additional inputs as
% option-value pairs. Ex 'numspreads', 7
%
% *Optional Inputs*
%
% * *numspreads* - |Integer| - The number of spreads you'd like to generate
% per cycle. *Default value: 7*
% * *freq* - |Integer| - How frequently cycles occur. *Default value: 10*
% * *bin* - |Array| - Logical array with binarized channel data. This is
% layered over synthetic spreads to add noise. 
% * *dist* - |Double| - Distance the spread grows per frame (typically a
% the epsilon value given by prep_dbscan works well for dist) *Default
% value: 65
% * *lifetime* - |Integer| - How many timepoints cells remain in the "on"
% or "active" state. *Default value: 3*
% * *seed* - |Integer|, |Double|, |String|, |Char| - a seed value to
% initialize the random number generator. Using the same seed value with
% the same data will yield the same results each time, making the number
% generator pseudo-random. Ex: 'potato', 6.022, 4815162342, "Orange"
% *Default value: []*
% * *maxsize* - |Double|, |Integer| - Maximum spread size, as multiples of
% the dist value. If dist is 60 and maxsize is 3 then the maximum size will be 180. *Default value: 2.8*
%
% *Output*
%
% * *bin_synth* - |Cell| - Cell array containing synthetic binarized data for the desired
% wells/xys
%
%% reformat
% Helper function to take ARCOS data organized clusters by time and
% reformat it to clusters by ID.
%
% *Input*
%
% * *cdata* - |Cell| - Cluster data. If using arcos_core, the main output of that
% method. If using arcos, the clust_by_time output. 
%
% *Output*
%
% * *out* - |Cell| - Clusters by ID. Analogous to the clust_by_id output of
% arcos.
%
%% binarize_robust
% Method with three different ways to robustly binarize channel data.
%
% # *znorm*
% # *robust* aka *selfrobust*
% # *specialrobust*
%
% *Inputs*
%
% * *raw_data* - |Cell| - Cell array containing processed data. Ex
% dataloc.d 
% * *xy* - |Array| - An array of integers representing indices of wells/xys
% to process. Supports discontinuous ranges - Ex (1:5, 10:20)
% * *ctrl* - |Array| - An array of integers representing indices of control
% wells/xys. 
%Supports discontinuous ranges - Ex (1:5, 10:20)
% * *chan* - |Char| - The channel to binarize.
% * *type* - |Char| - One of the 3 methods listed above. Ex 'znorm'. Be
% sure to exactly match the way it's typed above.
% * *perc* - |Double|, |Integer| - The percentage of the data by which it
% will be thresholded/binarized. Give as a full percentage and not a decimal. 80% should be
% given as 80, not 0.8. Accepts values as high as 200 and as low as 0. 
%
% *Ouput*
%
% *out* - |Cell| - Cell array where each cell is binarized data for the
% indexed well/xy.
%% Convert ARCOS Analysis Data
% arcosDF = Convert_Arcos_Analysis_Data(dataloc, channel, varargin)
% Converts ARCOS in given dataloc (version > 3 required) structure in a table (that can be used like a dataframe)
% 
% dataloc - dataloc struct or cell array of (can be mixed) dataloc structs and paths to datalocs
% channel - ARCOS channel you want to convert (ex: EKAR), can be cell array for each dataloc or single channel
% 
% Optional Inputs:
% tstartaftertx - how many HOURS before starting to consider data after the given treatment (default is 0)
% tmaxaftertx - how many HOURS after the tx are allowed to be considered? (default is whole movie after tx)
% aftertx - after which tx should data be considered? (default is first tx for that xy)
% exclude - treatments to exclude (default = none)
% txs -txs to consider (default = {'pTx','Tx1','Tx2','Tx3','Tx4'})
% 
% boxsize - how big should the square for "hotspot" detection be? (in uM) (will be boxsize uM by boxsize uM) (default = 30)
% smlthresh - time thresholds for short, medium, and long spread durations (default = [0,0.5,1]) given in hours
% 
% meanreplicates - mean technical replicates? (default = false)
% standardize - standardize the data? (default = false)
%
% Copyright. Nicholaus DeCuzzi. 2023
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
            [~,ix]=min(abs(slopes-1)); %Get the index of max_d for ideal eps (where slope of line tangent to that point is 1);
            eps = max_d(ix);
        end 
    function [eps,minpts] = prep_dbscan2(XCoord,YCoord,bin,varargin)
			inp.n = 11;
			nin = length(varargin);
			if rem(nin,2) ~= 0; warning('Additional inputs must be provided as option, value pairs'); end  %#ok<WNTAG>
			for s = 1:2:nin; inp.(lower(varargin{s})) = varargin{s+1}; end
			n = inp.n;
            xy = [XCoord,YCoord];
            P = 0;
            [~,d] = knnsearch(xy,xy,'K', n);
            d_real = sort(d(~isnan(d(:,n)),n));
            eps = mean(median(d_real)); %Median distance to 11th neighbor of all points
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
			if ischar(ch); ch = {ch}; end % assert ch is cell
			bin = cell(1,numel(xy));
			MegaDataHolder = [];
			%%Loop through wells, append channel data to MegaDataHolder
			for i = 1:size(raw_data,2)
				if ~isempty(raw_data{i})
					channel = raw_data{i}.data.(ch{1});
					MegaDataHolder = vertcat(MegaDataHolder, channel); %#ok<AGROW> 
				end
			end
			%%Get mean percentile of all well's channel data
			threshold = mean(prctile(MegaDataHolder,perc));
			%%Loop through wells again and binarize
			for i = 1:numel(xy)
				well = xy(i);
				channel = raw_data{well}.data.(ch{1});
				bin{well} = channel>threshold;
			end
		end 
		function out = pulse2bin(pulseAnalysisData,raw_data,channel)
    		%%Check for XYs with data and get their indicies
    		xy = 1:numel(pulseAnalysisData);
    		goodxys = ~cellfun(@isempty,pulseAnalysisData);% check to see if the input xys are good
    		xy = xy(goodxys);
    		
    		%%Loop through XYs
    		numXYs = numel(xy);
    		
    		out = {};
    		for iwell = 1:numXYs
        		well = xy(iwell);
        		%%Store pulse data and coords
        		pulseData = pulseAnalysisData{well}.data.(channel);
    		
        		% Make a data container
        		dataCont = false(size(raw_data{well}.data.XCoord,1),size(raw_data{well}.data.XCoord,2));
        		
                % Assign tracked points as zeros FIX - data is NaN unless real?
        		%dataCont(~isnan(raw_data{well}.data.XCoord)) = 0;
    		
        		%%Loop through pulse data
        		for i = 1:size(pulseData,1)
            		if ~isempty([pulseData(i).mpos])
                		for ii = 1:size(pulseData(i).mpos,1)
                    		%% Get start position and end (from dur of each pulse
                		    tStart = pulseData(i).mpos(ii) - ((pulseData(i).dur(ii)-1)/2);
                		    tEnd = pulseData(i).mpos(ii) + ((pulseData(i).dur(ii)-1)/2); % Duration includes the initial tp, so you have to subtract 1
                    		%% Fill in the these points with 1s!
                    		dataCont(i,tStart:tEnd) = 1;
                		end
            		end
        		end
        		out{well} = dataCont; %#ok<AGROW> 
    		end %XY loop
		end %convertPulseToBin
        function bin_synth = gensynth(raw_data,xy,varargin)
            %%Default parameters
            p.numspreads = 7; %How many spreads occur per cycle
            p.freq = 10; %Spread frequency
            p.bin = []; %Optional binarized data - useful for adding "noise"
            p.dist = 65; %Distance the spread grows per frame (epsilon from prep_dbscan is a useful metric)
            p.lifetime = 3; %How long a cell within a spread remains active before switching off
            p.seed = []; %Optional seed value for random number generator
            p.maxsize = 2.8; %Max spread size
            %%Setup
            nin = length(varargin);     %Check for even number of add'l inputs
            if rem(nin,2) ~= 0; warning('Additional inputs must be provided as option-value pairs'); end%Splits pairs to a structure
            for s = 1:2:nin; p.(lower(varargin{s})) = varargin{s+1}; end
            %%Random Number Generator Seed
			if ~isempty(p.seed)
				seed = mean(double(evalc('disp(p.seed)')));
				rng(seed)
			end%Use user-provided seed
            %%Pre-allocate output
            bin_synth = cell(1,length(xy));
            %%Loop through XYs (wells)
            for iwell = 1:numel(xy)
                well = xy(iwell);
                XCoord = raw_data{well}.data.XCoord; % Assign X and Y Coords from data
                YCoord = raw_data{well}.data.YCoord;
                binSynth = boolean(zeros(size(XCoord,1), size(XCoord,2)));  % Pre-allocate logical array
                %%Generate random start indices
				for r = 1:round(size(binSynth,2)/p.freq+1)
                    col = r*p.freq-(p.freq-1);
					for q = 1:p.numspreads %Generate random indices for starting pts
                        istartPts(q,r) = randi(size(XCoord,1)); %#ok<AGROW>
                        binSynth(istartPts(q,r),col) = 1;
					end
				end
				%%Starting parameters
                cnt = 1;
                sz = 1;
                %%Loop through time series
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
						%%Set points inactive if they've been alive for p.lifetime (Donut)
						
					end
                end
                %%Layer real data on top of synthetic
                if ~isempty(p.bin) %Use user-provided binary data
                    bin_synth{well} = logical(binSynth+p.bin{well});
                else
                    bin_synth{well} = logical(binSynth);
                end
            end %end well loop
		end 
		function out = reformat(cdata)
			%%Load data into struct
			timerange = size(cdata,2);
			max_id = cdata(end).newmax;
			clusters = repmat(struct('cid',[],'data',struct('time',{},'XYCoord',{},'id',{},'numpts',{},'bounds',{},'inactive',{},'area',{},'compl',{},'rocarea',{},'roccount',{}),'t_start',[],'t_end',[],'dur',[], 'maxarea',[],'maxcount',[]),max_id,1); %set up struct of structs
			for time = 1:timerange
				eps = cdata(time).eps;
				minpts = cdata(time).minpts;
				dtp = cdata(time).tracked; %dtp = data at timepoint
                if max([dtp(end).id]) > 0 % EDIT FIX HI ITS NICK BREAKING STUFF
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
			end
			%%Loop through substructs and remove empty entries
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
        end % end of reformat(data)
		
        function out = binarize_robust(raw_data,xy,ctrl,chan,type,perc)
			MegaDataHolder = []; %Array of channel data, vertically concatenated
			if perc > 200 || perc < 0; error('Perc should be a percentage <200, ex 85'); end
			perc = perc/100;
			%%Loop over all the xys and combine them into a mega data

			% This loop defined by xys or control wells you specify
			% Control wells should be conditions where the cells "do
            % nothing". EX: ERK inhibitor treated EKAR cells
            if isempty(ctrl); xy2bin = xy; else; xy2bin = ctrl; end

			numXY2bin = numel(xy2bin);
			for iXY = 1:numXY2bin
    			tXY = xy2bin(iXY);
    			if ~isempty(raw_data{tXY}) %check the XY isn't empty
    			if isfield(raw_data{tXY}, 'data') %check the XY has data
    			if isfield(raw_data{tXY}.data, chan) %check if that XY has the channel in it
        			MegaDataHolder = vertcat(MegaDataHolder, raw_data{tXY}.data.(chan)); %#ok<AGROW> %concat all data into one place
    			end %if the chan exists
    			end %if data exists
    			end %if raw_data is not empty
			end %end normdata collection loop

    		%% Calculate the method of normalization requested
			switch type
                case 'znorm' % Z-Score normalize to the whole dataset (or controls)
                    % Method: Subtract mean of all (or control) data then divide that by the std dev of all the data
                    % Use that std dev as threshold to determin on state
        			Subtract = mean(MegaDataHolder,'all','omitnan'); %get the mean of all the data
                    DivideBy = std(MegaDataHolder, 0, 'all', 'omitnan'); % get the std
                    DiffSigThresh = DivideBy; %Determine Value that shows enough change to be postive
                
                case 'meansub znorm' % Z-Score normalize to after self mean subtraction
                    % Method: Subtract each data's mean from itself, then find what is considered a std dev from that mean 
                    MegaDataHolder = MegaDataHolder - mean(MegaDataHolder,2,'omitnan'); %subtract the mean first
                    DiffSigThresh = std(MegaDataHolder, 0, 2, 'omitnan'); %get std dev of each data point
                    DiffSigThresh = mean(DiffSigThresh); % get the average std

                case 'minsub znorm' % Z-Score normalize using means of the datas 
                    % Method: Subtract each data's min from itself, subtract the mean, then divide that individual cell by its own std dev 
                    MegaDataHolder = MegaDataHolder - min(MegaDataHolder,[],2,'omitnan'); %subtract the min first
        			DiffSigThresh = std(MegaDataHolder, 0, 2, 'omitnan'); % get std dev
                    DiffSigThresh = mean(DiffSigThresh,"all",'omitnan');

                case 'robust' % Robust Scalar (scales to median and quantiles)
                    % Method: Get the median and quantiles of the data over time
                    % Subtract the median and divide by the IQR
                    Subtract = median(MegaDataHolder,2,'omitnan');
                    MegaDataHolder = MegaDataHolder - Subtract; % subtract the median of each data point
        			Quantilz = quantile(MegaDataHolder, (0:0.25:1),'all');
        			Q25 = Quantilz(2); Q75 = Quantilz(4); %get the 25th and 75th quantiles
        			DivideBy = Q75 - Q25; %get the interquantile range
                    DiffSigThresh = DivideBy; 

                case 'minmediansub robust' % Robust Scalar (scales to median and quantiles, after subtracting the mins)
                    % Method: Get the median and quantiles of the data over time
                    % Subtract the median and divide by the IQR
                    MegaDataHolder = MegaDataHolder - min(MegaDataHolder,[],2,'omitnan'); %subtract the min first
        			MegaDataHolder = MegaDataHolder - median(MegaDataHolder,2,'omitnan'); % then sub the median
        			Quantilz = quantile(MegaDataHolder, (0:0.25:1),2);
        			Q25 = Quantilz(:,2); Q75 = Quantilz(:,4); %get the 25th and 75th quantiles
        			DiffSigThresh = Q75 - Q25; %get the interquantile range
                    DiffSigThresh = median(DiffSigThresh,'omitnan');

                case 'mediansub robust' % Robust Scalar (scales to median and quantiles, after subtracting the mins)
                    % Method: Get the median and quantiles of the data over time
                    % Subtract the median and divide by the IQR
        			MegaDataHolder = MegaDataHolder - median(MegaDataHolder,2,'omitnan'); % sub the median
        			Quantilz = quantile(MegaDataHolder, (0:0.25:1),2);
        			Q25 = Quantilz(:,2); Q75 = Quantilz(:,4); %get the 25th and 75th quantiles
        			DiffSigThresh = Q75 - Q25; %get the interquantile range
                    DiffSigThresh = median(DiffSigThresh,'omitnan');

    			case 'specialrobust' % Robust Scalar (scales to median and quantiles)
        			Subtract = median(MegaDataHolder, 2,'omitnan'); %get the median for each cell, to self subtract
        			Quantilz = quantile(MegaDataHolder, (0:0.25:1),2);
        			Q25 = Quantilz(:,2); Q75 = Quantilz(:,4); %get the 25th and 75th quantiles
        			DiffSigThresh = median(Q75 - Q25,'omitnan'); %get the median interquantile range to divide all data by
                   
                otherwise
        			error('Normalization method %s given as input not recognized. \n', type)
            end

            %% Loop over all the xys again to normalize them
			numXYs = numel(xy);
			out = cell(1,numXYs);
            for iXY = 1:numXYs
                tXY = xy(iXY);
                if ~isempty(raw_data{tXY}) %check the XY isn't empty
                if isfield(raw_data{tXY}, 'data') %check the XY has data
                if isfield(raw_data{tXY}.data, chan) %check if that XY has the channel in it
                    tDat = []; %ensure data isn't used twice by accident
                    switch type
                        case 'znorm'
                            tDat = (((raw_data{tXY}.data.(chan)) - Subtract) ./ DivideBy) > ((raw_data{tXY}.data.(chan)) - Subtract) + DiffSigThresh*perc;  % is the data above the thresh?
                            tDat = double(tDat); % convert the data back to a double
                            tDat(isnan(raw_data{tXY}.data.(chan))) = nan; % find where the nan was before make the data nan again
                            out{tXY} = tDat; % give the data to the output structure
                        case {'minsub znorm','meansub znorm'}
                            tDat = raw_data{tXY}.data.(chan) > (mean(raw_data{tXY}.data.(chan),2,'omitnan') + (DiffSigThresh*perc));
                            tDat = double(tDat); % convert the data back to a double
                            tDat(isnan(raw_data{tXY}.data.(chan))) = nan; % find where the nan was before make the data nan again
                            out{tXY} = tDat; % give the data to the output structure
                        case {'minmediansub robust','mediansub robust'}
                            tDat = raw_data{tXY}.data.(chan) > (median(raw_data{tXY}.data.(chan),2,'omitnan') + (DiffSigThresh*perc));
                            tDat = double(tDat); % convert the data back to a double
                            tDat(isnan(raw_data{tXY}.data.(chan))) = nan; % find where the nan was before make the data nan again
                            out{tXY} = tDat; % give the data to the output structure
                        case 'specialrobust'
                            Subtract = median(raw_data{tXY}.data.(chan), 2,'omitnan'); %get the median for each cell, to self subtract
                            tDat = (((raw_data{tXY}.data.(chan)) - Subtract) ./ DivideBy) > (DivideBy*perc); %make normalized data
                            tDat = double(tDat); % convert the data back to a double
                            tDat(isnan(raw_data{tXY}.data.(chan))) = nan; % find where the nan was before make the data nan again
                            out{tXY} = tDat; % give the data to the output structure
                        otherwise
                            tDat = (((raw_data{tXY}.data.(chan)) - Subtract) ./ DivideBy) > (DivideBy*perc); %make normalized data
                            tDat = double(tDat); % convert the data back to a double
                            tDat(isnan(raw_data{tXY}.data.(chan))) = nan; % find where the nan was before make the data nan again
                            out{tXY} = tDat; % give the data to the output structure
                    end %end switch
					%mean or median of DivideBy value (divide it by 2),
					% use that number to threshold.
                end %if the chan exists
                end %if data exists
                end %if raw_data is not empty
            end %end second iXY loop
        end % end of binarize_robust()
        function arcosDF = Convert_Arcos_Analysis_Data(dataloc, channel, varargin)
        
        %treatments existing
        inp.tstartaftertx = [];
        inp.tmaxaftertx = [];
        inp.aftertx = [];
        inp.exclude = [];
        inp.txs = {'pTx','Tx1','Tx2','Tx3','Tx4'};
        
        inp.boxsize = 30; % how big should the square for "hotspot" detection be? (in uM) (will be boxsize uM by boxsize uM)
        inp.lmhthresh = [0,3,5]; % how many spreads need to start from a given box to be considered low, medium, or high? (numbers here are the min for that group)
        inp.smlthresh = [0,0.5,1]; % time thresholds for short, medium, and long spread durations
        
        inp.meanreplicates = false; % mean technical replicates?
        inp.standardize = false; %stdize the data?
        
        inp.alreadysizenormalized = true; % is the given data already normalized to uM (rather than pixels, which is the default output of arcos data)
        
        inp.varNames = {'exp',     'xy',  'full',  'txinfo', 'rep',   'freq',  'dur',  'maxarea','freq_by_region','region_sizexy','smlthresh','smlfractions'}; % variable names for the table
        inp.varTypes = {'string','string','string','string','string','double','double','double',      'cell',        'string',      'string',    'double'}; % variable types for the varNames in the table
        
        
        %% Check input varargin parameters
        nin = length(varargin); if rem(nin,2) ~= 0; warning('Additional inputs must be provided as option, value pairs'); end %Splits pairs to a structure
        for s = 1:2:nin; inp.(lower(varargin{s})) = varargin{s+1}; end; clear nin s varargin;
        
        if ~isempty(inp.exclude) && ~iscell(inp.exclude); inp.exclude = {inp.exclude}; end
        if isstruct(dataloc); dataloc = {dataloc}; end % put dataloc in a cell
        nData = numel(dataloc); % how many datalocs are there
        if ~iscell(channel); channel = {channel}; end % put channel in a cell
        if numel(channel) < nData; channel = repmat(channel,1,nData); end % make a channel for every dataloc if you didn't already
        
        % Map the treatments for each dataloc
        arcosDF = table;
        for iData = 1:nData
            dHold = dataloc{iData};
            % Check if dataloc is a dataloc struct or is a path to a dataloc file
            if isstring(dHold) || ischar(dHold) 
                dHold = load(dHold);
                dHold = dHold.dataloc;
            end
            dHold = makeArcosDF(dHold, channel{iData}, inp); % make the live cell dataframe
            arcosDF = [arcosDF; dHold]; % append the dHold dataframe to datalocDF
            clear dHold;
        end 
        
        end % ARCOS Dataframe constructor
        
        function arcosDataFrame = makeArcosDF(dataloc, channel, inp)
        
        %% Set up everything for analysis
        
        inp.tktm = 60/dataloc.movieinfo.tsamp; %set up how many tps per hour
        if ~isempty(inp.tstartaftertx);  inp.tstartaftertx = floor(inp.tstartaftertx * inp.tktm); end 
        if ~isempty(inp.tmaxaftertx); inp.tmaxaftertx = floor(inp.tmaxaftertx * inp.tktm); end %convert respond time from hours to tps
        
        txData = ct_maptx(dataloc.platemapd.idx,'time',false,'reduce',false,'expar',inp.txs);
        txNames = fieldnames(txData);
        
        if ~isempty(inp.exclude)
            for iExclude = 1:numel(inp.exclude)
                tossThese = contains(txNames,inp.exclude{1});
                tossThese = txNames(tossThese);
                txData = rmfield(txData,tossThese);
                clear tossThese;
            end
            txNames = fieldnames(txData);
        end
        
        arcosDataFrame = table('Size',[0,numel(inp.varNames)],'VariableNames',inp.varNames,'VariableTypes', inp.varTypes); % set up the arcos dataframe
        
        %% Set up squares for hotspot detection
        
        if inp.alreadysizenormalized % if the arcos data is already normalized to uM for coordinates 
            xTrueBox = inp.boxsize; % how many uM should a box be (in the x)? uM already
            yTrueBox = inp.boxsize; % how many uM should a box be (in the y)? uM already
            xuMsize = dataloc.movieinfo.PixNumX * dataloc.movieinfo.PixSizeX; % how big is the x (in uM)
            yuMsize = dataloc.movieinfo.PixNumY * dataloc.movieinfo.PixSizeY; % how big is the y (in uM)
            xBoxes = (0:xTrueBox:xuMsize); % box range for x in uM
            yBoxes = (0:yTrueBox:yuMsize); % box range for y in uM
        
        else % if the arcos data is NOT already normalized to uM for coordinates 
            xTrueBox = round(inp.boxsize / dataloc.movieinfo.PixSizeX); % how many pixels are in a p.boxsize uMs (in the x)? must be converted from uM to pixels 
            yTrueBox = round(inp.boxsize / dataloc.movieinfo.PixSizeY); % how many pixels are in a p.boxsize uMs (in the y)? must be converted from uM to pixels
            xBoxes = (0:xTrueBox:dataloc.movieinfo.PixNumX); % box range for x in pixels
            yBoxes = (0:yTrueBox:dataloc.movieinfo.PixNumY); % box range for y in pixels
        end % make the boxes for the hotspot measurements 
        
        edgesForYou = {xBoxes,yBoxes}; % put em together
        
        
        %% Loop through the treatments
        for iTx = 1:size(txNames,1)
            tTx = txNames{iTx};
            repCounter = 0;
            tempDF = table('Size',[0,numel(inp.varNames)],'VariableNames',inp.varNames,'VariableTypes', inp.varTypes); % set up temp arcos dataframe
        
            for iXY = txData.(tTx).xy'
                if ~isempty(dataloc.(['arcos_',channel])(iXY)) && ~isempty(dataloc.(['arcos_',channel]){iXY})
                
                repCounter = repCounter + 1; % go up in rep
        
                tempDF2 = table('Size',[1,numel(inp.varNames)],'VariableNames',inp.varNames,'VariableTypes', inp.varTypes); % set up temp arcos dataframe
                
                tempDF2.rep = repCounter; % assign rep number
                tempDF2.exp = dataloc.file.base; % assign the experiment name
                tempDF2.xy = num2str(iXY); % the xy
                tempDF2.full = replace(tTx,'_',' + '); % add the full sorted name
                % tempDataFrame.cell = mappedTxs.(use2Map{iMap}).tx(1).name; % I currently don't care about density stuff FIX?
                % thisTx = structfun(@(x) ismember(iXY,x.xy),txInfo); % get a logical for the tx
        
                tempDF2.region_sizexy = [num2str(xTrueBox),', ',num2str(yTrueBox)]; % add the box sizes for hot spot info
                %tempDF2.lmhthresh = [num2str(inp.lmhthresh(1)),', ',num2str(inp.lmhthresh(2)),', ',num2str(inp.lmhthresh(3))]; % add the thresholds for low, medium, and hot spot info
                tempDF2.smlthresh = [num2str(inp.smlthresh(1)),', ',num2str(inp.smlthresh(2)),', ',num2str(inp.smlthresh(3))]; % add the thresholds for short, medium, and long spread info
        
                allTxs = [txData.(tTx).tx.time]'; % get the treatment times
                allTxs = allTxs(allTxs > 0); % ignore pretreatments
                % get the tp of the last treatment (or otherwise)
                if isempty(allTxs); tStart = 1;
                else
                    if ~isempty(inp.aftertx)
                        if numel(allTxs) < inp.aftertx; tStart = allTxs(end);
                        else; tStart = allTxs(inp.aftertx);
                        end
                    else; tStart = allTxs(end);
                    end
                end
        
                if ~isempty(inp.tmaxaftertx) % adjust max time considered
                    if (tStart + inp.tmaxaftertx) > size(dataloc.d{iXY}.data.XCoord,2)
                        tEnd = size(dataloc.d{iXY}.data.XCoord,2);
                    else; tEnd = tStart + inp.tmaxaftertx;
                    end
                else; tEnd = size(dataloc.d{iXY}.data.XCoord,2);
                end
        
                if ~isempty(inp.tstartaftertx) %adjust t start 
                    if (inp.tstartaftertx + tStart) > tEnd; warning('Your tstart after tx is after the end of your tmax after tx, ignoring tstart after tx');
                    else; tStart = tStart + inp.tstartaftertx;
                    end
                end
        
                %% Collect all of the treatment info and determine which data are allowed given the times provided
                txInfo = [txData.(tTx).tx]; %pull the one the xy is a part of and assign it
                txInfo2 = '';
                for iTTx = 1:length(txInfo)
                    txInfo2 = [txInfo2, ' ', num2str(txInfo(iTTx).dose),txInfo(iTTx).dunit, ' ', txInfo(iTTx).name];
                end
                txInfo2 = txInfo2(2:end);
                tempDF2.txinfo = {txInfo2}; % NOW put it into the dataframe
                
                aData = dataloc.(['arcos_',channel]){iXY}; % pull that xys arcos data
        
                sStart = [aData.t_start]'; % get when spreads start
                sEnd = [aData.t_end]'; %when spreads end
        
                sGoodRange(:,1) = (tStart < sStart); % only keep data that falls within the allowed time span
                sGoodRange(:,2) = (sEnd < tEnd); % only keep data that falls within the allowed time span
                sGoodRange = all(sGoodRange,2); % only keep data that falls within the allowed time span
        
                aData = aData(sGoodRange); % keep the arcos data that falls in the given range
        
                %% calculate SPREAD freq
                a = sum(sGoodRange); % take how many spreads occur in the given time 
                a = a / ((tEnd - tStart) / inp.tktm); % and divide it all by the allowed time span
                a = a / (xuMsize * yuMsize);  % x um/px * numXpix * y um/px * numYpix normalize the data for the image size
                tempDF2.freq = a * 1000000; %  * 10^6 (convert um^2 to mm^2) and append the data
                clear a;
        
                %% Divide the image into inp.boxsize um x inp.boxsize um regions and get the spread "rate" (spreads per hr) per square for that image
                a = {aData.data}'; % get the spread data
                a = cellfun(@(x){nanmean(x(1).XYCoord,1)},a); 
                a = cell2mat(a); % get the first xy of the spread
                if isempty(a); a = [1,1]; end
                h = histogram2(a(:,1),a(:,2),'XBinEdges',edgesForYou{1},'YBinEdges',edgesForYou{2});
                countz = h.Values(:); % pull the counts per square from the data
                close(gcf)
                tempDF2.freq_by_region = {(countz / ((tEnd - tStart) / inp.tktm))'}; % get the fraction of counts per group (out of all squares)
                clear a h countz totalSegs;
        
                %% mean and distribution of spread durations
                a = [aData.dur]' / inp.tktm; % get the durs and make them into hours
                totalSegs = [((inp.smlthresh(1) < a) & (a <= inp.smlthresh(2))), ((inp.smlthresh(2) < a) & (a <= inp.smlthresh(3))),inp.smlthresh(3) < a];
                tempDF2.smlfractions = sum(totalSegs,1)/length(a); % what fraction of spreads are short vs medium vs long
                tempDF2.dur = mean(a,'all','omitnan'); % get the mean duration of spreads and append the data
                clear totalSegs a;
        
                %% mean max spread size
                a = [aData.maxarea]'; %get the spread data
                if ~inp.alreadysizenormalized
                    a = a * (dataloc.movieinfo.PixSizeX * dataloc.movieinfo.PixSizeY);  % x um/px  * y um/px 
                end
                tempDF2.maxarea = mean(a,'all','omitnan'); % get the mean max spread size and append the data
                clear a;
                
                clear sGoodRange;
        
                if inp.meanreplicates
                    tempDF = [tempDF; tempDF2];
                else; arcosDataFrame = [arcosDataFrame; tempDF2];
                end
        
                end % empty d{} check
            end % xy loop
        
            %% Mean the replicates if desired
            if inp.meanreplicates && ~isempty(tempDF)
                arcosDataFrame = [arcosDataFrame; tempDF(1,:)]; %transfer the info from the first xy in the list (exp, full, txinfo, rep, lmhthresh, boxsizexy, smlthresh)
                arcosDataFrame.xy(end) = join([tempDF.xy]', ', ');
                arcosDataFrame.freq(end) = mean(tempDF.freq,1,'omitnan');
                fbr = cell2mat(tempDF.freq_by_region);
                arcosDataFrame.freq_by_region{end} = {fbr(:)'}; clear fbr;
                arcosDataFrame.dur(end) = mean(tempDF.dur,1,'omitnan');
                arcosDataFrame.smlfractions(end,:) = mean(tempDF.smlfractions,1,'omitnan');
                arcosDataFrame.maxarea(end) = mean(tempDF.maxarea,1,'omitnan');
            end
        
        end %tx loop
        
        if exist('arcosDataFrame','var')
            if inp.standardize
                arcosDataFrame = (arcosDataFrame-nanmean(arcosDataFrame))./std(arcosDataFrame,'omitnan'); %standardize vals
            end
        else; arcosDataFrame =[];
        end
        end % makeArcosDF function% of methods
    end
    
end