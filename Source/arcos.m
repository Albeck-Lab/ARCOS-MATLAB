%% ARCOS
% ARCOS, or Automated Recognition of Collective Signaling is a MATLAB-based
% adaptation and expansion of Maciej Dobrzynski's <https://github.com/dmattek/ARCOS R Package> of the same name.
% 
% 
%% Inputs
% * *data* - |cell| - The data to be processed. Each cell is a well/xy.
% * *xy* - |array| - Array of well indices (integers) to process. Can handle
% discontinuous indices. Ex: (1:5, 11:15)
% * *ch* - |char| - The identifier of the channel to process. Ex: 'nEKAR'.
% Must be a 'character', not a "string". 
% * _varargin_ - |option value pairs| - accepts optional inputs as option-value pairs.
%%% Optional Inputs
% * *bin* - |cell| - User-provided binarized data. Each cell is well/xy binary data as logical arrays. Rows are cells (observations) and columns are timepoints. *Default: []*
% * *bin_perc* - |Double|- If auto-binarizing. the percentile by which to threshold and binarize the data. *Default: []*
% * *eps* - |cell| - User-provided epsilon values for dbscan. *Default: []*
% * *minpts* - |cell|- User-provided minpts values for dbscan. *Default: []*
% * *verbose* - |Logical|, |Boolean| - Toggle verbose logging. *Default: true*
%% Outputs
% * *clust_by_time* - |cell| - Clusters organized by timepoint
% * *clust_by_id* - |cell| - Clusters organized by cluster ID
% * *binaries* - |cell| - Binarization data used to cluster
% * *warnings* - |cell| - Warnings that can indicate poor clustering
%% Examples
% See the Demos folder for a variety of examples

function [clust_by_time, clust_by_id, binaries,warnings,labels] = arcos(data,xy,ch,varargin)
	%% Optional Parameters
	p.bin = []; %user-provided binarized data %%check if it's the same size as the X and Y coord data
	p.bin_perc = []; %Percentile for threshold binarization
	p.eps = {[]};
	p.minpts = {[]};
	p.pixsize = [1 1];
	p.verbose = true;
	p.debug = false;
	%% Prep varargin struct
	nin = length(varargin);
	if rem(nin,2) ~= 0; warning('Additional inputs must be provided as option, value pairs');  end
	for s = 1:2:nin; p.(lower(varargin{s})) = varargin{s+1};   end
	%% DBSAN Prep
    if isnumeric(p.eps); p.eps = {p.eps}; end %assert eps is cell
    if isnumeric(p.minpts); p.minpts = {p.minpts}; end %assert eps is cell
	%% XY Prep
	if isempty(xy); xy = (1:size(data,2)); end
	%% Preallocate out
	clust_by_time = cell(1,numel(xy));
	clust_by_id = cell(1,numel(xy));
	binaries = cell(1,numel(xy));
	optionalOuts = cell(1,numel(xy));
	warnings = struct('frame_warnings',{},'excess_nans',{});
	%% Check XYs
	goodxys = ~arrayfun(@(x)isempty(data{x}),xy);% check to see if the input xys are good
	xy = xy(goodxys);
	%% Channel prep
	if ischar(ch); ch = {ch}; end % assert ch is cell
	if isempty(p.bin_perc) && isempty(p.bin)
		warning("Optional parameter 'bin_perc' not set. Binarizing data using 80th percentile threshold")
		p.bin_perc = 80;
	end
	if isempty(p.bin)
		if ~isempty(p.bin_perc)
			bin = arcos_utils.binarize(data,xy,ch,p.bin_perc); %Use simple binarization if no user-provided binarized data
		end
	else
		bin = p.bin;
	end
	%% Loop through XY
	numXYs = numel(xy);
    for iwell = 1:numXYs
		well = xy(iwell);
		if p.debug == true;disp(append('Processing well ',string(well))); end
		%% Define XCoord and YCoord
		XCoord = data{well}.data.XCoord;
		YCoord = data{well}.data.YCoord;
		assert(~isempty(XCoord), 'No x coordinate data detected');
		assert(~isempty(YCoord), 'No y coordinate data detected');
		%% Warn if too many NaNs
		sum_nans = sum(isnan(XCoord),'all');
		sz = numel(XCoord);
		nans_thr = 70; %If this percentage of the data is NaNs it'll get logged in warnings
		if sum_nans/sz*100 > nans_thr
			warnings(well).excess_nans = sum_nans/sz*100;
			if p.verbose == true; warning(append("Excessive NaNs detected in well ", string(well))); end
		end
		%% Format and assign eps and minpts (if given)
		if numel(p.eps) == 1; eps = p.eps{1}; else; eps = p.eps{well}; end 
		if numel(p.minpts) == 1; minpts = p.minpts{1}; else; minpts = p.minpts{well}; end
        %%Do the arcos functions
        %[clust_by_time{well},warnings(well).frame_warnings] = arcos_core(XCoord,YCoord,bin{well},'eps',eps,'minpts',minpts, 'verbose', p.verbose, 'debug', p.debug, 'well', well,'pixsize', p.pixsize);
		
		[labels,warnings{well}.frame_warnings,optionalOut] = arcos_core(XCoord,YCoord,bin{well},'epsilon',eps,'minpts',minpts, 'verbose', p.verbose, 'debug', p.debug, 'well', well,'pixsize', p.pixsize);
		optionalOuts{well} = optionalOut;
		clust_by_time{well} = optionalOut{6};
		clust_by_id{well} = arcos_utils.reformat(clust_by_time{well});
		
		binaries = bin;
    end %well loop
end %wrapper function end
		






















































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

