function out = arcos(data,xy,ch,varargin)
    %% Optional Parameters
    p.upbin = []; %user-provided binarized data %%check if it's the same size as the X and Y coord data
    p.bin_perc = []; %Percentile for threshold binarization. 
    p.eps = [];
    p.minpts = [];
    %% Prep varargin struct
    nin = length(varargin);
    if rem(nin,2) ~= 0
        warning('Additional inputs must be provided as option, value pairs');  
    end
    for s = 1:2:nin
        p.(lower(varargin{s})) = varargin{s+1};   
    end
    %% Preallocate out
    out = cell(1,size(data,2));
    %% Loop through XY
    for well = xy(1):xy(2)
        %% Progress bar setup
        if well==1
            bar = waitbar(well/(1+xy(2)-xy(1)),append('Processing well ',int2str(well),' of ', int2str((1+xy(2)-xy(1)))));
        else
            waitbar(well/(1+xy(2)-xy(1)),bar,append('Processing well ',int2str(well),' of ', int2str((1+xy(2)-xy(1)))));
        end
        %% Define XCoord and YCoord
        XCoord = data{well}.data.XCoord;
        YCoord = data{well}.data.YCoord;
        %% Channel switch
        switch ch
           case 'EKAR'
               channel = data{well}.data.EKAR;
           case 'CFP_Nuc'
               channel = data{well}.data.CFP_Nuc;
           case 'YFP_Nuc'
               channel = data{well}.data.YFP_Nuc;
           case 'nEKAR'
               channel = data{well}.data.nEKAR;
           otherwise
               error('Invalid channel name');
        end
        %% Check inputs
        if isempty(p.upbin)
            warning('Using simple threshold binarization');
            bin = arcos_utils.binarize(channel,p.bin_perc); %Use simple binarization if no user-provided binarized data
        else
            warning('Using user-provided binarization');
            bin = p.upbin{well}; %Use user-provided binarization
        end
        assert(~isempty(XCoord), 'No x coordinate data detected');
        assert(~isempty(YCoord), 'No y coordinate data detected');
        assert(~isempty(bin), 'No binary data detected');
        assert(numel(XCoord) == numel(YCoord), 'XCoords and YCoords not equal')
        if ~isempty(p.eps)
            warning('User-specified epsilon detected');
        end
        %% Preallocate cdata, initialize newmax
        cdata = cell(5,size(XCoord,2)); %Preallocation of output
        newmax = 0;
        %% Loop through time
        for time = 1:size(XCoord,2)
            if isempty(p.eps)
                [eps,minpts] = arcos_utils.prep_dbscan3(XCoord(:,time),YCoord(:,time),bin(:,time)); %Calculate ad hoc epsilon and minpts
            end
            activeXY = [XCoord(bin(:,time)==1,time), YCoord(bin(:,time)==1,time)]; %XY coordinates for active cells
            cdata{1,time} = clustering(activeXY, eps, minpts); %Get untracked data
            %% Track if possible
            if time==1
                [cdata{2,time}, newmax] = tracking(cdata{1,time},cdata{1,time},cdata{1,time},eps, cdata{1,time}(end).id);
            elseif cdata{1,time}(end).id==0
                cdata{2,time} = [];
            else
                [cdata{2,time}, newmax] = tracking(cdata{1,time},cdata{2,time-1},cdata{1,time-1},eps, newmax); %Get 'tracked' data
            end
            cdata{3,time} = eps; %Store epsilon value used
            cdata{4,time} = minpts; %Store minpts value used
            cdata{5,time} = newmax; %Highest cluster ID for this frame
        end
        %% Perform analysis and assign to out cell
        out{1,well} = cdata;
        out{2,well} = arcos_analysis.analyze(XCoord,YCoord,cdata);
    end
    close(bar)
end %wrapper function end
function out = clustering(activeXY, eps, minpts)
    out = struct('points',0,'id',0,'hull',0,'area',0);
    if (isempty(activeXY))
        return
    end
    clusters = dbscan(activeXY, eps, minpts);
    for cl = 1:max(clusters)
        pts = activeXY(clusters==cl,:);
        if size(pts,1)>2
            [hull,area] = boundary(pts);
            out(cl).points = pts; % points in that cluster
            out(cl).id = cl; %cluster identity
            out(cl).hull = hull;
            out(cl).area = area;
        else
            continue
        end 
    end
end %clustering end
function [tracks,newmax] = tracking(sCurr,sPrev,bPrev,eps,newmax)
    %Unpack structs into arrays
    dCurr = unpack(sCurr); %XY and ID for clusters in curr
    dPrev = unpack(sPrev); %XY and ID for clusters in prev
    %Check if the requisite data is present
    if numel(dPrev)<=3 
        dPrev = unpack(bPrev);
    end
    if numel(dCurr)<=3 
        tracks = []; %If not, return empty
        return
    end
    
    [idx,d] = knnsearch(dPrev(:,1:2),dCurr(:,1:2)); %Indices and distances of current's neighbors in previous
    isClose = d <= eps;
    dCurr(isClose,4)= dPrev(idx(isClose),3);
    for i = 1:max(dCurr(:,3))
        cluster = dCurr(:,3)==i;        %Logical map for points in currently evaluated cluster
        id = dCurr(cluster,3:4);        %cluster ids for those points
        n = numel(unique(id(id(:,2)>0,2)));         %Get the number of unique cluster ids for the current cluster >0 (non-ignored)
        if n > 1                            %more than one past cluster
            dCurr(cluster,5) = 0;
        elseif n == 1                       %One past cluster
            dCurr(cluster,4) = unique(id(id(:,2)>0,2));
            dCurr(cluster,5) = 0;
        elseif sum(id(:,2))==0              %No past clusters
            if max(max(dCurr(:,4)))+1 > newmax
                newmax = max(max(dCurr(:,4)))+1;
            end
            dCurr(cluster,4)= newmax;
            dCurr(cluster,5) = 1;
        end
    end
    for i = 1:max(dCurr(:,4))
        %Hull points based on new cluster assignment
        newclust = dCurr(:,4)==i;
        if sum(newclust)>2
            %[hull,area] = convhull(dCurr(newclust,1:2));
            [hull,area] = boundary(dCurr(newclust,1:2));
        else
            hull = [];
            area = [];
        end
        tracks(i).points = dCurr(newclust,1:2);  %#ok<AGROW> %Points that make up the cluster
        tracks(i).id = dCurr(newclust,3:4); %#ok<AGROW> %original cluster ID, new cluster ID
        tracks(i).hull = hull; %#ok<AGROW> %Indices of points that encompass all other points
        tracks(i).area = area; %#ok<AGROW> %Area of the convex hull
        if dCurr(newclust,5) > 0
            tracks(i).new = 1; %#ok<AGROW> %Flag it reassigned cluster is new
        else
            tracks(i).new = 0; %#ok<AGROW>
        end
    end
    map = false(1,size(tracks,2));
    for i = 1:size(tracks,2)
        if isempty(tracks(i).points) || isempty(tracks(i).id)
            map(i) = 1;
        end
    end
    tracks(map) = [];
end %tracking end
function xy = unpack(d)
    xy = zeros(1,2);
    for i = 1:size(d,2) %Loop through struct elements
        for ii = 1:size(d(i).points) %Loop through elements of struct elements
            xy(end+1,1:2) = d(i).points(ii,:);
            if size(d(i).id,2) > 1
                xy(end,3) = d(i).id(ii,2);
                disp('');
            else
                xy(end,3) = d(i).id;
            end
        end
    end
    xy(1,:) = [];
end %unpack function end









%Default for simple binarization - add flag to wrapper for more robust
%binarization

%Input flag to Allow/search for merges/splits
%Take completed dataset and pass over it looking for merges/splits
%Segregate clusters (or don't)

%Set up structure in arcos procedure to do a second pass loop - performs
%actions on output data (cdata).
%QA pass basically. Find out if merge/splitting makes sense. 
%Loop through data cluster-by-cluster rather than time-by-time.





% MP -------------------------------------------------
%Here's the outline and notes that I've got for now.  It's not
%   perfect - let me know if there are things that don't make sense
%   or where you don't see how to implement it.
%   
%   (Maybe just pass in the two "cdata" time slices to compare,
%       since we don't use the rest of cdata in here.)      
%   
%   Get all Past Nearest Neighbors for Current Points
%   Get Past Cluster IDs for Neighbors within Epsilon
%       e.g. posCurr(:,4) = 0;  (or initialize posCurr with zeros)
%            isClose = d <= epsilon;  %Logical index
%            posCurr( isClose ) = posPrev( idx(isClose) ); %Assign
%
        %   FOR each Current Cluster
        %       (Maybe use a separate index for temporary cluster data to write, nextID)
        %       IF maps to more than one Past Cluster
        %           Assign each Current point to its Neighbor's Cluster ID
        %           Store as separate clusters with newly assigned Past IDs
        %               e.g. clust{nextID, 1} = XYs(posCurr(:,4)==PastID); 
        %                    clust{nextID, 2} = PastID; ... 
        %                    nextID = nextID + 1; ...
        %       IF maps to exactly one Past Cluster
        %           Store all points in Current Cluster as Past ID
        %               e.g. clust{nextID, 2} = PastID; ...
        %       IF maps to no Past Clusters
        %           Assign a new cluster ID, and maybe flag as a Start
        %              e.g.  maxID = maxID + 1; %(tdata grows)
        %                    clust{nextID, 2} = maxID; 
        %                    clust{nextID, 5} = 0; ...
        %                    nextID = nextID + 1; ...
        %       (Could sort by ID column after assignment, if helpful)
        %
        %   FOR each Past Cluster (to find Splits)
        %       Find Current Clusters mapped to this Past Cluster
        %       IF more than one, keep the one closest (average distance?),
        %           and make new IDs for the other(s), flagging lineage
        %           with ID from the parent Past Cluster
        %               e.g. tdata{newID, 2} = newID;  %(tdata may not need the self-ID column, if it's in ID order) 
        %                    tdata{newID, 5} = pastID; ...
        %       IF exactly one, copy temporary clust data to tdata
        %           e.g. tdata{pastID, :} = clust{matchingID, :};
        %           (Maybe initialize tdata column ~5 to NaN by default)
        %       IF none, flag past Cluster as an End (this might be useful)
        %           e.g. tdata{pastID, 5} = -1;
        %   
        % ----------------------------------------------------
%Add column to final output that has lineage tracking - flag if it's a new
%cluster started from nothing- split from something else or neither of
%those- was linked
% If new or split - 
%Once tracking loop is done, do a verification / lineage loop to track
%inheritance / etc
%
% Part 1 - 151 - copy in structure of "if more than 1 previous cluster connected,
% assign each point to whatever closest prev clust is. If only 1 prev clust
% connected, assign all points to prev clust.
%if no prev clust connected, keep them all as zeros
%Keep current frame cluster assignment and previous cluster assignment

% Part 2 - lineage tracking

%ARCOS Output
%Row 1: Cluster data - untracked
%Row 2: Cluster data - tracked
%Row 3: Epsilon value used
%Row 4: MinPts value used
%Row 5: Highest cluster value assigned
