function [clust_by_time, clust_by_id, binaries] = arcos(data,xy,ch,varargin)
    %% Optional Parameters
    p.bin = []; %user-provided binarized data %%check if it's the same size as the X and Y coord data
    p.bin_perc = []; %Percentile for threshold binarization
    p.eps = {};
    p.minpts = {};
    %% Prep varargin struct
    nin = length(varargin);
	if rem(nin,2) ~= 0; warning('Additional inputs must be provided as option, value pairs');  end
	for s = 1:2:nin; p.(lower(varargin{s})) = varargin{s+1};   end
	%% Channel prep
    if ischar(ch); ch = {ch}; end % assert ch is cell
    %% Preallocate out
    %out = cell(1,size(data,2));
	clust_by_time = cell(1,length(xy));
	clust_by_id = cell(1,length(xy));
	binaries = cell(1,length(xy));
    %% Check XYs
    goodxys = ~arrayfun(@(x)isempty(data{x}.cellindex),xy);% check to see if the input xys are good
    xy = xy(goodxys);
    %% Loop through XY
	numXYs = numel(xy);
    for iwell = 1:numXYs
        well = xy(iwell);
        %% Define XCoord and YCoord
        XCoord = data{well}.data.XCoord;
        YCoord = data{well}.data.YCoord;
		assert(~isempty(XCoord), 'No x coordinate data detected');
        assert(~isempty(YCoord), 'No y coordinate data detected');
        %% Channel selection
        channel = data{well}.data.(ch{1}); %Create wrapper function to loop through desired channels       
        %% Setup: Binarization
		if isempty(p.bin_perc) && isempty(p.bin)
            warning("Optional parameter 'bin_perc' not set. Binarizing data using 80th percentile threshold")
            p.bin_perc = 80;
		end
        if isempty(p.bin)
            if ~isempty(p.bin_perc)
                bin = arcos_utils.binarize(channel,p.bin_perc); %Use simple binarization if no user-provided binarized data
            else
               error("No binarized data has been provided, so please specify a percentile threshold to binarize data. Ex 'bin_perc', 80"); 
            end
        else
            bin = p.bin{well}; %Use user-provided binarization
        end
        
        %% Format outputs
		if ~isempty(p.eps); eps = p.eps{well}; else eps = []; end 
		if ~isempty(p.minpts); minpts = p.minpts{well}; else minpts = [];end 
        clust_by_time{well} = arcos_core(XCoord,YCoord,bin,'eps',eps,'minpts',minpts);
		clust_by_id{well} = reformat(clust_by_time{well});
		binaries{well} = bin;
        %out{2,well} = reformat(cdata); %Clusters by ID
        %out{3,well}= bin; %Binarization data
    end
end %wrapper function end
function out = reformat(cdata)
    %% Load data into struct
    timerange = size(cdata,2);
    max_id = cdata(end).newmax;
    clusters = repmat(struct('cid',[],'data',struct('time',{},'XYCoord',{},'id',{},'numpts',{},'bounds',{},'inactive',{},'area',{},'compl',{},'rocarea',{},'roccount',{}),'t_start',[],'t_end',[],'dur',[], 'maxarea',[],'maxcount',[]),max_id,1); %set up struct of structs
    for time = 1:timerange
        dtp = cdata(time).tracked; %dtp = data at timepoint
        for cluster = 1:size(dtp,2)
            id = mode(dtp(cluster).id(:,2)); %May not need the mode part... compare timeit results with and without mode.
            clusters(id).cid = id;
            clusters(id).data(time).time = time;
            clusters(id).data(time).XYCoord = dtp(cluster).XYCoord;
            clusters(id).data(time).id = dtp(cluster).id;
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














































%Take completed dataset and pass over it looking for merges/splits
%Segregate clusters (or don't)

%Set up structure in arcos procedure to do a second pass loop - performs
%actions on output data (cdata). 
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

