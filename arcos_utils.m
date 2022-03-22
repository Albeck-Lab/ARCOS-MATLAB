classdef arcos_utils
    methods(Static)
        function formatr(filename,XCoord,YCoord,bin)
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
            minpts = ndims(XCoord)*2;
            vals = zeros(1,time(2)+1-time(1));
            for t = time(1):time(2)
                [~,d] = knnsearch([XCoord(:,t),YCoord(:,t)], [XCoord(:,t), YCoord(:,t)],'K', minpts+1); %k-nearest neighbors search
                d = max(d,[],2); %Biased toward greater distances as opposed to average of k-nearest
                max_d = sort(d);
                scaled = max_d * length(max_d)/max(max_d); %Scale the data
                smoothed = smoothdata(scaled,'gaussian'); %Smooth it
                slopes = gradient(smoothed); %Take first derivative
                [~,ix]=min(abs(slopes-1)); %Get the index of avg_d for ideal eps (where slope of line tangent to that point is 1);
                eps = max_d(ix);
                vals(t) = eps;
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
            for k = n : -1 : 1
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
        function [bin,threshold] = binarize(ch,perc)
            threshold = mean(prctile(ch,perc));
            bin = ch>threshold;
        end 
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
                bin = logical(binSynth+p.bin);
            else
                bin = logical(binSynth);
            end
        end 
    end
end