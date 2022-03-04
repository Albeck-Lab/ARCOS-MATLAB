%% ARCOS Utils
% A collection of utilities for ARCOS
%%% Format R
% Formats data for ARCOS R package processing
%%% Prep DBSCAN
% Provides optimal epsilon and minpts values for DBSCAN
%%% Binarize
% Simple binarization function
%%% Gen Synth
% Generates synthetic spread events
classdef arcos_utils
    methods(Static)
        function formatr(filename,XCoord,YCoord,bin)
            array = zeros(size(bin,1)*size(bin,2),6);
            sz = size(bin,1);
            for col = 1:size(bin,2)
                for row = 1:size(bin,1)
                    if col >1
                        ind = row+sz*(col-1); 
                    else
                        ind = row;
                    end
                    array(ind,1) = col;
                    array(ind,2) = XCoord(row,col);
                    array(ind,3) = YCoord(row,col);
                    array(ind,4) = bin(row,col);
                    array(ind,5) = row;
                end
            end
            table = array2table(array);
            table.Properties.VariableNames(1:5) = ["t","x","y","m","id"];
            writetable(table,append(filename,'.csv'));
        end %EOF
        function [minpts, eps] = prep_dbscan(XCoord, YCoord)
            minpts = ndims(XCoord)*2;
            [~,d] = knnsearch([XCoord,YCoord], [XCoord, YCoord],'K', minpts+1); %k-nearest neighbors search
            d = max(d,[],2); %Biased toward greater distances as opposed to average of k-nearest
            max_d = sort(d);
            scaled = max_d * length(max_d)/max(max_d); %Scale the data
            smoothed = smoothdata(scaled,'gaussian'); %Smooth it
            slopes = gradient(smoothed); %Take first derivative
            [~,ix]=min(abs(slopes-1)); %Get the index of avg_d for ideal eps (where slope of line tangent to that point is 1);
            eps = max_d(ix);
        end %EOF
        function [bin,threshold] = binarize(ch,perc)
            threshold = mean(prctile(ch,perc));
            bin = ch>threshold;
        end %EOF
        function bin = gensynth(XCoord,YCoord,varargin)
            p.numspreads = 7; %How many spreads occur per cycle
            p.freq = 10; %Spread frequency
            p.bin = []; %Optional binarized data - useful for adding "noise"
            p.dist = 1; %Distance the spread grows per frame (epsilon from prep_dbscan is a useful metric)
            p.lifetime = 3; %How long a cell within a spread remains active before switching off
            p.seed = []; %Optional seed value for random number generator
            p.maxsize = 2.8; %Max spread size
            %%Setup
            nin = length(varargin);     %Check for even number of add'l inputs
            if rem(nin,2) ~= 0
                warning('Additional inputs must be provided as option-value pairs');  
            end%Splits pairs to a structure
            for s = 1:2:nin
                p.(lower(varargin{s})) = varargin{s+1};   
            end
            
            assert(sum(size(XCoord)==size(YCoord),'all')==2,'XCoord and YCoord must be equal size.');
            binSynth = boolean(zeros(size(XCoord,1), size(XCoord,2)));
            if ~isempty(p.seed) %Use user-provided seed
            rng(p.seed);
            end
            for r = 1:round(size(binSynth,2)/p.freq)
                col = r*p.freq-(p.freq-1);
                for q = 1:p.numspreads %Generate random indices for starting pts
                    istartPts(q,r) = randi(size(XCoord,1)); %#ok<AGROW>
                    binSynth(istartPts(q,r),col) = 1;
                end
            end
            %binSynth(istartPts(:,),1:round(size(binSynth,2)/p.freq):end) = 1; %Set start pts to active
           
            cnt = 1;
            sz = 1;
            for time = 1:size(XCoord,2)
                if time==1 || mod(time,p.freq) == 0 && cnt <=size(istartPts,2)
                    startPts = [XCoord(istartPts(:,cnt),time),YCoord(istartPts(:,cnt),time)]; %Get coords for start pts
                    cnt = cnt+1;
                    sz = 1;
                end
                pts = [XCoord(:,time),YCoord(:,time)]; %Get all xy coords at current time
                [idx,d] = rangesearch(pts,startPts,p.dist*sz); %Get points within dist*time of start pts
                sz = sz+1;
                for i = 1:size(idx)
                    binSynth(idx{i}(d{i}<= p.dist*p.maxsize),time) = 1; %Set points within dist*time of start pts to active
                    for row = 1:size(binSynth,1)
                        if sum(binSynth(row,:))>=p.lifetime
                        binSynth(row,time)=0; %Set cells to inactive if they've been active > lifetime
                        end
                    end
                end
            end
            if ~isempty(p.bin) %Use user-provided binary data
                bin = binSynth+p.bin;
            else
                bin = binSynth;
            end
        end %EOF
    end
end