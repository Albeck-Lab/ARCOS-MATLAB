%% ARCOS
% A limited feature set port of the R package "ARCOS"
%
% "Automated Recognition of Collective Signalling (ARCOS) is an R package
% to identify collective spatial events in time series data."
%
% <https://github.com/dmattek/ARCOS Original ARCOS GitHub>
%
% <<https://github.com/dmattek/ARCOS/raw/main/man/figures/README-ex1plotTS-1.png>>
%
% <<https://github.com/dmattek/ARCOS/raw/main/man/figures/README-ex1plotColl-1.png>> 
% 
%% Inputs
% * *XCoord* - |2D array| - x coordinates of cells. Each row is a cell,
% each column is a timepoint
% * *YCoord* - |2D array| - y coordinates of cells. Each row is a cell,
% each column is a timepoint
% * *bin* - |2D array| - binarized data indicating 'active' (1) and
% 'inactive' (0) cells
% * _varargin_ - |option value pairs| - accepts optional inputs as option-value pairs.
%%% Optional Inputs
% * *eps* - |single| , |double| - Search radius used to find points in the
% dbscan clustering algorithm. Must be greater than zero. *Default: Calculated for every timepoint*
% * *minpts* - |integer| - Minimum number of neighbors required for a
% core point, specified as a positive integer. Must be a positive integer greater than zero. 
% *Default: ndims(XCoord)x2*
% * *time* - |1x2 array of ints| - Desired start and end time points. First
% element is start time, second element is end time. *Default: Earliest
% time to latest time index* 
% * *genplot* - |boolean| - Specify whether to plot the results. *Default:
% false*
%% Outputs
% *cdata* - outputs a cell array with columns representing timepoints. Each
% cell contains a cell array of structs representing spread events for that
% timepoint.
%
% Each struct contains the following fields
%
% * *pts* - the xy coordinates of cells in the spreading event
% * *hull* - indices of points that make up the spread's convex hull
% * *area* - the computer area of the convex hull
%
%% Examples
% *Using default parameters*
%
%   cdata = arcos(XCoord,YCoord,bin);
%
% *Using optional parameters, specifying epsilon of 55, minpts 5*
%
%   cdata = arcos(XCoord,YCoord,bin,'eps',55,'minpts',5);
%% See Also
% * <https://www.mathworks.com/help/stats/dbscan.html#mw_a35e3831-bd7a-437f-b128-889a3a444aa2
% Density-based spatial clustering of applications with noise (DBSCAN)>
% * <https://www.mathworks.com/help/matlab/ref/convhull.html Convex Hull>
% * <https://www.aaai.org/Papers/KDD/1996/KDD96-037.pdf A Density-Based
% Algorithm for Discovering Clusters>
%% To Do
% * Test for false-positive rate
% * Test for divergence events
% * Test for convergence events

%%
function cdata = arcos(XCoord, YCoord, bin, varargin)
    p.eps = [];
    p.minpts = [];
    p.sm = 1; %Tracking speed multiplier. Searches a wider radius for neighbor clusters. *Set at your own risk - may affect accuracy of tracking algorithm*
    %%Prepare additional inputs
    nin = length(varargin);     %Check for even number of add'l inputs
    if rem(nin,2) ~= 0
        warning('Additional inputs must be provided as option, value pairs');  
    end%Splits pairs to a structure
    for s = 1:2:nin
        p.(lower(varargin{s})) = varargin{s+1};   
    end

    %%Perform checks
    assert(~isempty(XCoord), 'No x coordinate data detected');
    assert(~isempty(YCoord), 'No y coordinate data detected');
    assert(numel(XCoord) == numel(YCoord), 'XCoords and YCoords not equal')
    if isempty(bin)
       warning('No binarized data supplied');
    end


    %%Loop through time
    cdata = cell(1,size(XCoord,2)); %Preallocate cdata for speed
    for time = 1:size(XCoord,2)
        if isempty(p.eps) || isempty(p.minpts)
            [minpts, eps] = arcos_prep_dbscan(XCoord(:,time),YCoord(:,time));
        else
            minpts = p.minpts;
            eps = p.eps;
        end
        activeXY = [XCoord(bin(:,time)==1,time), YCoord(bin(:,time)==1,time)];
        cdata{1,time} = arcos_core(activeXY, eps, minpts);
        cdata{2,time} = eps;
        cdata{3,time} = arcos_track(cdata,time, p.sm);
    end
end %wrapper function end

%%Get minpts and epsilon
function [minpts, eps] = arcos_prep_dbscan(XCoord, YCoord)
    %%Calculate minpts
    % Default is 4 for 2D data according to Ester et al., 1996)
    minpts = ndims(XCoord)*2;
    %%Calculate eps
    [~,d] = knnsearch([XCoord,YCoord], [XCoord, YCoord],'K', minpts+1); %k-nearest neighbors search
    d = max(d,[],2); %Biased toward greater distances as opposed to average of k-nearest
    max_d = sort(d);
    scaled = max_d * length(max_d)/max(max_d); %Scale the data
    smoothed = smoothdata(scaled,'gaussian'); %Smooth it
    slopes = gradient(smoothed); %Take first derivative
    [~,ix]=min(abs(slopes-1)); %Get the index of avg_d for ideal eps (where slope of line tangent to that point is 1);
    eps = max_d(ix);
end %prep dbscan function end

%%ARCOS Core
% 
function events = arcos_core(activeXY, eps, minpts)
    assert(eps>0, 'eps must be greater than 0') %These checks might need to happen earlier
    assert(minpts>1, 'minpts must be greater than 1')
    if (isempty(activeXY))
        events = [];
        return
    end
    if minpts <=2 
        warning('minpts less than 3 may yield inaccurate results'); 
    end
    clusters = dbscan(activeXY, eps, minpts);
    events = cell(max(clusters),4); %Preallocate cell array for speed
    for cl = 1:max(clusters)
        pts = activeXY(clusters==cl,:);
        if size(pts,1)>2
            [hull, area] = convhull(pts);
            events{cl,1} = num2cell(pts); % points in that cluster
            events{cl,2} = cl; %cluster identity
            events{cl,3} = num2cell(hull);
            events{cl,4} = area;
        else
            continue
        end 
    end
end %arcos_core end

%%getPos - a utility for arcos_track for returning xy coords
function pos = getPos(data)
    lpos = [0,0,0];
    for cl = 1:size(data,1) %Loop through clusters and get their coords
        clCurr = cell2mat(data{cl,1});
        for p = 1:size(clCurr,1) %Append all the cluster points into one matrix
            lpos(end+1,1) = clCurr(p,1);
            lpos(end,2) = clCurr(p,2);
            lpos(end,3) = data{cl,2};
        end
    end
    lpos(1,:)=[]; %remove the empty first line
    pos = lpos;
end
function tdata = arcos_track(cdata,t,sm)
    if (t-1==0) %If t-1 is out of bounds
        tdata = cdata{1,t};
    else
        eps = cdata{2,t}; %epsilon for current timepoint
        posCurr = getPos(cdata{1, t}); %Get all xy for clusters at time current
        posPrev = getPos(cdata{3, t-1}); %Get xy for previous clusters
        [idx,d] = knnsearch(posPrev(:,1:2),posCurr(:,1:2)); %the closest neighbor to point in posCurr is posPrev(idx) with distance d
        if size(idx)>0
                for pt = 1:size(idx)
                   if d(pt) <= eps*sm %If the point is within epsilon distance, set cid to prev cid
                       posCurr(pt,4) = posPrev(idx(pt),3);
                   else
                       posCurr(pt,4) = 0;
                   end
                end

            clusters = cell(max(posCurr(:,4)),4); %generate empty cell array
            for cl = 1:max(posCurr(:,4)) %format and output data
                cls = posCurr(:,4)==cl; %logical map for rows of current cl
                clustXY = posCurr(cls==1,1:2); %xys for those rows
                clusters{cl,1} = num2cell(clustXY);
                clusters{cl,2} = cl;
                if size(clustXY,1)>2
                    [hull, area] = convhull(clustXY);
                    clusters{cl,3} = num2cell(hull);
                    clusters{cl,4} = area;
                end
            end
            tdata = clusters; %Returning an empty array until I can get this working.
        else
            tdata = [];
        end
    end
end
