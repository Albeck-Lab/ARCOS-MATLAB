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
% option-value pairs. Ex: 'epsilon', 60
%
% *Optional Inputs*
%
% * *epsilon* - |Double| - Optional epsilon value for DBSCAN. Automatically
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
% * *epsilon* - |Double| - epsilon value for DBSCAN
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
% * *epsilon* - |Double| - epsilon value used to cluster the currently indexed
% data. Used in knnsearch algorithm
% * *maxLabel* - |Integer| - Max cluster ID assigned in the _previous_ tracking
% iteration
%
% *Outputs*
%
% * *tracks* - |Struct| - Structure containing XY coordinates, cluster ID
% lineage and a flag indicating whether a cluster reassignment is unique
% * *maxLabel* - |Integer| - Max cluster ID assigned in the _current_ tracking
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
function [labelTracked,warnings,optionalOut] = arcos_core(XCoord,YCoord,bin,varargin)

	%Set default parameters
	p.epsilon = [];
	p.minpts = [];
	p.verbose = true;
	p.debug = false;
	p.well = [];
	p.pixsize = [1 1];

	%Parse inputs
	nin = length(varargin);
	if rem(nin,2) ~= 0; warning('Additional inputs must be provided as option, value pairs'); end  %#ok<WNTAG>
	for s = 1:2:nin; p.(lower(varargin{s})) = varargin{s+1}; end

	[nCells, nTime] = size(XCoord);

	runPrep = [isempty(p.epsilon), isempty(p.minpts)];
    %%Preallocate cdata, initialize maxLabel, preallocate warnings

    maxLabel = 0;
	maxLabels = zeros(nTime,1);
	warnings = struct('too_sparse',{});

	%Scale coordinates to calibrated dimensions
    XCoord = XCoord.*p.pixsize(1);
    YCoord = YCoord.*p.pixsize(2);
	
	allActive = cell(1,nTime);
	labelUntracked = zeros(size(XCoord));
	labelTracked = zeros(size(XCoord));
	epsilon = zeros([nTime,1]);
	minpts = zeros([nTime,1]);

    fHold = tiledlayout(1,2); % NICK FIX
	%%Time loop
	for time = 1:nTime
		%%Setup: epsilon and Minpts
		if any(runPrep); [epsilon(time),minpts(time)] = arcos_utils.prep_dbscan2(XCoord(:,time),YCoord(:,time),bin(:,time)); end
		if ~runPrep(1);  epsilon(time) = p.epsilon;  end  %Override with user-provided epsilon
		if ~runPrep(2);  minpts(time) = p.minpts; end  %Override with user-provided minpts


		%%Log warning if epsilon is abnormal
		if epsilon(time) >150
			warnings(time).too_sparse = "High epsilon. Check data for low cell density";
		end % Warning if epsilon is abnormally high
		
		
		%%Setup and run Clustering
		
		active = bin(:,time);
        activeXY = [XCoord(active,time), YCoord(active,time)];
		allActive{time} = active;
        labelUntracked(active,time) = clustering(activeXY, epsilon(time), minpts(time));
		

        %%Track if possible
		%Tracking is retroactive from the current frame. 
		% If the current frame is the first frame, untracked becomes tracked.
		if time>1
            [labelTracked(:,time),previousTracked] = tracking(...
				labelUntracked(:,time),...
				labelTracked(:,time-1),...
				[XCoord(:,time),YCoord(:,time)],...
				[XCoord(:,time-1),YCoord(:,time-1)],...
				allActive{time-1},...
				epsilon(time),maxLabel);

				%Rectify labelUntracked for previous if start points
				%detected
				delta = labelTracked(:,time-1) ~= previousTracked;
				labelTracked(:,time-1) = previousTracked;

				labelUntracked(delta,time-1) = previousTracked(delta)+max(labelUntracked(:,time-1));
%{ 
for plotting spreads in real time              
                if sum(previousTracked) < 0
                    clrz = hsv(numel(unique(previousTracked))-1);
                    clrz = clrz(randperm(size(clrz,1)),:);
                    clrz = [0.95,0.95,0.95; clrz];
                    nexttile(1);
                    gscatter(labelTracked(:,time-1),previousXY(:,2),previousTracked,clrz,".",15);
                    nexttile(2);
                    gscatter(previousXY(:,1),previousXY(:,2),previousTracked,clrz,".",15);
                    pause(0.1)
                end
%}
		else
            labelTracked(active,time) = labelUntracked(active,time);
		end
		%Update record of maximum cluster ID across this frame
		maxLabel = max(labelTracked(:));
		maxLabels(time) = maxLabel;
	end


	num_high_epsilon = numel([warnings(:).too_sparse]);
	if num_high_epsilon > 0 && p.verbose
	    warning(append("Well: ", string(p.well), " - ", string(num_high_epsilon), ' timepoints had higher than expected epsilon values and may not cluster well'))
	    disp("See the 'warnings' output for more information")
	end
	optionalOut = {};
	optionalOut{1} = allActive; %Per-timepoint vectors 
	optionalOut{2} = labelUntracked;
	optionalOut{3} = epsilon;
	optionalOut{4} = minpts;
	optionalOut{5} = maxLabels;
	
	optionalOut{6} = formatLegacy(labelTracked,labelUntracked,XCoord,YCoord,epsilon,minpts,maxLabels); % This is clust_by_time
end


function label = clustering(activeXY, epsilon, minpts)
	%Short-circuit on empty input
	if isempty(activeXY)
		label = [];
		return; 
	end
	%Perform main dbscan clustering
    label = dbscan(activeXY, epsilon, minpts);
	label(label==-1) = 0; %Change outlier cluster IDs to 0 (unclustered)
end



function [labelTracked,previousTracked] = tracking(currentUntracked,previousTracked,currentXY,previousXY,previousActive,epsilon,maxLabel)
	labelTracked = zeros(size(currentUntracked)); %FIXME - Maybe not initialize, because Untracke labels are not 
	%	constrained with respect to previous tracked labels.  Could randomly equal previous labels?
	
	%Get indices and distances of current points' neighbors in previous frame
    [idxPreviousNeighbors,dPreviousNeighbors] = knnsearch(previousXY,currentXY); 
    isClose = dPreviousNeighbors <= epsilon; 	%Flag if neighbors are within epsilon
	%	Assign current point the cluster label of nearest previous point, IF within epsilon
	isClose(currentUntracked==0) = 0; 
	% FIXME the idx output of knnsearch is the index of the XY coordinates
	% and cannot be used to index into label matrices

	idxPreviousClose = idxPreviousNeighbors(isClose);
	labelTracked(isClose) = previousTracked(idxPreviousClose);

	currentLabelUnique = unique(currentUntracked(currentUntracked > 0))';  %Get list of current cluster labels
    for i = currentLabelUnique  	%Operating cluster-by-cluster
        clusterCurrentMask = currentUntracked==i;  %Mask of points in this Cluster

		%Get (unique) list of new Labels mapped to this Cluster
		newLabelList = unique(labelTracked(clusterCurrentMask));
		newLabelList = newLabelList(newLabelList > 0);

		nNew = length(newLabelList);  %Number of Labels mapped to this Cluster

		%Adjust Label assignments to be consistent across current Clusters
        if nNew == 0  		%IF mapped to no past Clusters
            %{
				Search again for outlier points (individual active cells) that are 
				within some radius of the current cluster, and if found, add this/these 
				point(s) to the previous frame's cluster list, labeled with the 
				current cluster's newly generated label.

				do inpolygon search
				if point in polygon
            %}

            maxLabel = maxLabel + 1;
            if sum(clusterCurrentMask)>2
                currentX = currentXY(clusterCurrentMask,1);
                currentY = currentXY(clusterCurrentMask,2);
                previousX = previousXY(:,1);
                previousY = previousXY(:,2);

                p = boundary(currentX,currentY);
                in = inpolygon(previousX,previousY,currentX(p),currentY(p));
                if numel(previousX(and(in,previousActive))) > 0
                    previousTracked(and(in,previousActive)) = maxLabel;
                end
            end
            %Initialize a new Cluster Label and assign to all points in this Cluster
            labelTracked(clusterCurrentMask) = maxLabel;

        elseif nNew == 1  	%IF only one past Cluster is mapped
            %Assing ALL points in this Cluster to that Label (including points mapped to nothing)
            labelTracked(clusterCurrentMask) = newLabelList;

        elseif nNew > 1		%IF this Cluster maps to multiple past Clusters
            %Make linear index mask of previous Labels mapped to this Cluster
            prevMask = find(ismember(previousTracked,newLabelList));
            %Get indices of (previous) nearest neighbors to points in this Cluster
            idxPrev = knnsearch(previousXY(prevMask,:), currentXY(clusterCurrentMask,:));
            %Assign each point to the cluster of its nearest neighbor (only from newLabelList...)
            labelTracked(clusterCurrentMask) = previousTracked(prevMask(idxPrev));
        end
    end
end

function clust_by_time = formatLegacy(labelTracked,labelUntracked,XCoord,YCoord,epsilon,minpts,maxLabels)
	% formatLegacy is a reformatting function to convert arcos_core data to
	
	[nCells,nTime] = size(labelUntracked);

	% Each row of clust_by_time is a timepoint.
	% At each timepoint we have untracked and tracked structs, eps,
	% minpts and a newmax value.

	clust_by_time = struct("untracked",{},"tracked",{},"eps",{},"minpts",{},"newmax",{});


	for time = 1:nTime
		
		%tTracked = labelTracked(:,time);						%Tracked labels for time
		%tUntracked = labelUntracked(:,time);					%Untracked labels for time
		%active = allActive{time};								%Active cells for time;
		%activeXY = [XCoord(active,time),YCoord(active,time)];	%XYCoords for time

		
																% Build untracked
																% Within untracked each row is a cluster
																% Each cluster has XYCoords and an id


		
		untracked = struct("XYCoord",{},"id",{});

		%mapActiveUntracked = tUntracked ~= 0;					%Mapping of active cells in untracked
		%activeUntracked = tUntracked(mapActiveUntracked);		%Labels of active cells in untracked
		uniqueUntracked = unique(labelUntracked(labelUntracked(:,time)~=0,time));				%Unique  untracked labels ~=0
		for iu = 1:length(uniqueUntracked)						%For every unique untracked label
			u = uniqueUntracked(iu);							%Current indexed unique label
			%clusterXY = activeXY(activeUntracked == u,:);		%Coordinates for that label
			clusterXY = [XCoord(labelUntracked(:,time) == u,time),YCoord(labelUntracked(:,time)==u,time)];
			untracked(iu).XYCoord = clusterXY;
			untracked(iu).id = u;
		end

																% Tracked is similar to untracked with some notable differences
																% In tracked, the id has 2 elements: the untracked id and the
																% tracked id
																% Tracked also has a flag for whether the cluster id 
																% is "new" (less than newmax from the previous timepoint)
																% The first element of id may not be used by any downstream
																% function and is being considered for deprecation. 
																% Newmax is also unused and may be deprecated


		tracked = struct("XYCoord",{},"id",{},"new",{});


		%mapActiveTracked = tTracked ~= 0;
		%activeTracked = tTracked(mapActiveTracked);
		uniqueTracked = unique(labelTracked(labelTracked(:,time)~=0,time));
		
		for it = 1:length(uniqueTracked)
			t = uniqueTracked(it);
			%clusterXY = activeXY(activeTracked == t,:);
			clusterXY = [XCoord(labelTracked(:,time) == t,time),YCoord(labelTracked(:,time)==t,time)];
			tracked(it).XYCoord = clusterXY;
			trackedID = repmat(t,height(clusterXY),1); % Repeat t for the number of cells we have
			% untracked labels for the current tracked labels

			untrackedID = labelUntracked(labelTracked(:,time)==t,time);
			%tracked(it).id = [activeUntracked(activeTracked==t),trackedID];
			tracked(it).id = [untrackedID,trackedID];
			if (time > 1) && (t < maxLabels(time-1))
				tracked(it).new = 1;
			else
				tracked(it).new = 0;
			end
			
		end
		
		


		clust_by_time(time).untracked = untracked;
		clust_by_time(time).tracked = tracked;
		clust_by_time(time).eps = epsilon(time);
		clust_by_time(time).minpts = minpts(time);
		clust_by_time(time).newmax = maxLabels(time);
	end
end





