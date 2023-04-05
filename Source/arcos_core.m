%% ARCOS Core
% The core clustering functions of ARCOS
% 
%% ARCOS_Core
%
% *Inputs*
%
% * *XCoord* - |Double array| - X coordinates of tracked cells
% * *YCoord* - |Double array| - Y coordinates of tracked cells
% * *bin* - |Logical array| - binary data
% * *_varargin_* - |Option-value pairs| - Accepts optional inputs as
% option-value pairs. Ex: 'eps', 60
%
% *Optional Inputs*
%
% * *eps* - |Double| - Optional epsilon value for DBSCAN. Automatically
% calculated if left empty.
% * *minpts* - |Integer| - Optional minpts value for DBSCAN. Automatically
% calculated if left empty.
%
% *Outputs*
%
% * *cdata* - |struct| - Cluster data organized by time
% * *warnings* - |Table| - Warnings indicating where poor clustering may
% have occurred
%
%% Clustering
% Clustering of active cells via DBSCAN
%
% *Inputs*
%
% * *activeXY* - |Double array| - X and Y coordinates for active cells
% * *eps* - |Double| - Epsilon value for DBSCAN
% * *minpts* - |Integer| - Minpts value for DBSCAN
% * _varargin_ - |Option-value pair| - Accepts addition inputs as
% option-value pairs. Ex 'debug', true. 
%
% * Optional Inputs*
%
% * *debug* - |Boolean|, |Logical| - Debug logging toggle. *Default value:
% false*
% * *verbose* - |Boolean|, |Logical| - Verbose logging toggle. *Default
% value: true*
%
% *Ouputs*
%
% * *out* - |Struct| - struct containing XY coordinates of clustered cells
% and their cluster IDs
%
%
%% Tracking
% Tracking method that reassigns cluster IDs based on knnsearch
%
% *Inputs*
%
% * *sCurr* - |Struct| - Struct containing the currently indexed untracked
% cluster data
% * *sPrev* - |Struct| - Struct containing the previously indexed tracked
% cluster data
% * *bPrev* - |Struct| - Struct containing the previously indexed untracked
% cluster data
% * *eps* - |Double| - Epsilon value used to cluster the currently indexed
% data. Used in knnsearch algorithm
% * *newmax* - |Integer| - Max cluster ID assigned in the _previous_ tracking
% iteration
%
% *Outputs*
%
% * *tracks* - |Struct| - Structure containing XY coordinates, cluster ID
% lineage and a flag indicating whether a cluster reassignment is unique
% * *newmax* - |Integer| - Max cluster ID assigned in the _current_ tracking
% iteration
%
%% Unpack
% Helper function to unpack structs into arrays.
%
% *Inputs*
%
% * *d* - |Struct| - Struct to be unpacked into an array
%
% *Outputs*
%
% * *xy* - |Array| - Array of X and Y coordinates unpacked from the
% inputted struct
%
%% Examples
% See "core_only.mlx" in the Demos folder to use ARCOS Core by itself
%
% Otherwise, see default_params.mlx in the Demos folder.
%% See also
% dbscan, knnsearch
function [labelTracked,warnings] = arcos_core(XCoord,YCoord,bin,varargin)
	p.eps = [];
	p.minpts = [];
	p.verbose = true;
	p.debug = false;
	p.well = [];
	p.pixsize = [1 1];
	nin = length(varargin);
	if rem(nin,2) ~= 0; warning('Additional inputs must be provided as option, value pairs'); end  %#ok<WNTAG>
	for s = 1:2:nin; p.(lower(varargin{s})) = varargin{s+1}; end


	runPrep = [isempty(p.eps), isempty(p.minpts)];
    %%Preallocate cdata, initialize newmax, preallocate warnings
    newmax = 0;
	
	warnings = struct('too_sparse',{});


    XCoord = XCoord.*p.pixsize(1);
    YCoord = YCoord.*p.pixsize(2);
	
	allActive = {};
	labelUntracked = zeros(size(XCoord));
	labelTracked = zeros(size(XCoord));
	eps = zeros([size(XCoord,2),1]);
	minpts = zeros([size(XCoord,2),1]);

	%%Time loop
	for time = 1:size(XCoord,2)
		%%Setup: Eps and Minpts
		if any(runPrep); [eps(time),minpts(time)] = arcos_utils.prep_dbscan2(XCoord(:,time),YCoord(:,time),bin(:,time)); end
		if ~runPrep(1);  eps(time) = p.eps;  end  %Override with user-provided eps
		if ~runPrep(2);  minpts(time) = p.minpts; end  %Override with user-provided minpts


		%%Log warning if eps is abnormal
		if eps(time) >150
			warnings(time).too_sparse = "High Epsilon. Check data for low cell density";
		end % Warning if epsilon is abnormally high
		
		
		%%Setup and run Clustering
		
		active = bin(:,time);
        activeXY = [XCoord(active,time), YCoord(active,time)];
		allActive{time} = activeXY;
        labelUntracked(active,time) = clustering(activeXY, eps(time), minpts(time));
		

        %%Track if possible
		%Tracking is retroactive from the current frame. If the current
		%frame is the first frame, untracked becomes tracked.
		if time>1
            [labelTracked(active,time)] = tracking(...
				labelUntracked(:,time),...
				labelTracked(:,time-1),...
				allActive{time},...
				allActive{time-1},...
				eps(time),newmax);
		else
            labelTracked(active,time) = labelUntracked(active,time);
		end
		newmax = max(labelTracked(active,time));
	end


	num_high_eps = numel([warnings(:).too_sparse]);
	if num_high_eps > 0 && p.verbose
	    warning(append("Well: ", string(p.well), " - ", string(num_high_eps), ' timepoints had higher than expected epsilon values and may not cluster well'))
	    disp("See the 'warnings' output for more information")
	end


end


function label = clustering(activeXY, eps, minpts)
    if isempty(activeXY); return; end
    label = dbscan(activeXY, eps, minpts);
    outliers = label == -1;
    maxID = max(label(~outliers));
    if isempty(maxID); maxID=0; end
    outlierclust = maxID+1:maxID+sum(outliers);
    label(outliers) = outlierclust';
end


function labelTracked = tracking(currentUntracked,previousTracked,currentXY,previousXY,eps,newmax)
    labelTracked = currentUntracked;

    [idxPreviousNeighbors,dPreviousNeighbors] = knnsearch(previousXY(:,1:2),currentXY(:,1:2)); %Indices and distances of current's neighbors in previous
    close = dPreviousNeighbors <= eps; %Distances for neighbors within epsilon
	idxPreviousClose = idxPreviousNeighbors(close);
	labelPreviousClose = previousTracked(idxPreviousClose); %Cluster labels of previous points within eps of each current point
	

	currentLabelUnique = unique(currentUntracked)';
    for i = currentLabelUnique
        clusterCurrentMap = currentUntracked==i;
		ilabelCurrent = currentUntracked(clusterCurrentMap);       %cluster ids for those points





        
        n = numel(unique(id));         %Get the number of unique cluster ids for the current cluster
        if sum(id)==0
            if max(labelTracked)+1 > newmax
                newmax = max(max(dCurr(:,4)))+1;
            end
            labelTracked(clusterCurrentMap) = newmax;
        end
    end

end

