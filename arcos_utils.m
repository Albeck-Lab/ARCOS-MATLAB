%% ARCOS Utils
%%% Format R
%%% Prep DBSCAN
%%% Binarize
%%% Gen Synth
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
        function bin = gensynth(XCoord,YCoord,interval,varargin)
            p.numspreads = 7;
            p.bin = [];
            p.dist = 1;
            p.lifetime = 2;
            p.seed = [];
            p.maxsize = 2.8;
            %%Setup
            nin = length(varargin);     %Check for even number of add'l inputs
            if rem(nin,2) ~= 0
                warning('Additional inputs must be provided as option-value pairs');  
            end%Splits pairs to a structure
            for s = 1:2:nin
                p.(lower(varargin{s})) = varargin{s+1};   
            end
            assert(sum(size(XCoord)==size(YCoord),'all')==2,'XCoord and YCoord must be equal size.');
            binIn = boolean(zeros(size(XCoord,1), size(XCoord,2)));
            if ~isempty(p.seed) %Use user-provided seed
            rng(p.seed);
            end
            if isempty(interval)
                interval(1)=2;
                interval(2)=size(XCoord,2);
            end
            for q = 1:p.numspreads %Generate random indices for starting pts
                istartPts(q,1) = randi(size(XCoord,1)); %#ok<AGROW>
            end
            binIn(istartPts,1) = 1; %Set start pts to active
            startPts = [XCoord(istartPts,1),YCoord(istartPts,1)]; %Get coords for start pts
            for time = interval(1):interval(2)
                pts = [XCoord(:,time),YCoord(:,time)]; %Get all xy coords at current time
                [idx,d] = rangesearch(pts,startPts,p.dist*time); %Get points within dist*time of start pts
                for i = 1:size(idx)
                    binIn(idx{i}(d{i}<= p.dist*p.maxsize),time) = 1; %Set points within dist*time of start pts to active
                    for row = 1:size(binIn,1)
                        if sum(binIn(row,:))>=p.lifetime
                        binIn(row,time)=0; %Set cells to inactive if they've been active > lifetime
                        end
                    end
                end
            end
            if ~isempty(p.bin) %Use user-provided binary data
                bin = binIn+p.bin;
            else
                bin = binIn;
            end
        end %EOF
    end
end