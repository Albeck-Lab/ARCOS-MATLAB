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
function [cdata,warnings] = arcos_core(XCoord,YCoord,bin,varargin)
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
	cdata = struct('untracked',{},'tracked',{},'eps',{},'minpts',{},'newmax',{});
	warnings = struct('too_sparse',{});
    XCoord = XCoord.*p.pixsize(1);
    YCoord = YCoord.*p.pixsize(2);
	%%Time loop
	for time = 1:size(XCoord,2)
		%%Setup: Eps and Minpts
		if any(runPrep); [epst,minptst] = arcos_utils.prep_dbscan2(XCoord(:,time),YCoord(:,time),bin(:,time)); end
		if ~runPrep(1);  epst = p.eps;  end  %Override with user-provided eps
		if ~runPrep(2);  minptst = p.minpts; end  %Override with user-provided minpts
		%%Log warning if eps is abnormal
		if epst >150
			warnings(time).too_sparse = "High Epsilon. Check data for low cell density";
		end % Warning if epsilon is abnormally high
		%%Setup and run Clustering
<<<<<<< Updated upstream
        activeXY = [XCoord(bin(:,time),time), YCoord(bin(:,time),time)]; %XY coordinates for active cells EDIT FIX ND - this can give a nan?      
        cdata(time).untracked = clustering(activeXY, epst, minptst);
        %%Track if possible
		if time==1
            [cdata(time).tracked,newmax] = tracking(cdata(time).untracked,cdata(time).untracked,cdata(time).untracked,epst,cdata(time).untracked(end).id);
        elseif cdata(time).untracked(end).id==0
            cdata(time).tracked = [];
		else
            [cdata(time).tracked, newmax] = tracking(cdata(time).untracked,cdata(time-1).tracked,cdata(time-1).untracked,epst,newmax);
=======
        activeXY = [XCoord(bin(:,time),time), YCoord(bin(:,time),time)]; %XY coordinates for active cells
        nanmap = isnan(activeXY);
        activeXY = [activeXY(~nanmap(:,1),1),activeXY(~nanmap(:,2),2)];
        [cdata(time).untracked,cdata(time).untracked2] = clustering(activeXY, epst, minptst);
        %%Track if possible
        if cdata(time).untracked(end).id==0 %No active points, no clusters
            cdata(time).tracked = [];
        end
		if time>1
             %{
            %% Reassign cluster IDs
            if isempty(cdata(time-1).tracked)
                for q = 1:size(cdata(time).untracked,2)
                    cdata(time).untracked(q).id = cdata(time).untracked(q).id + newmax;
                end
                newmax = max(unique([cdata(time).untracked.id]));
            end
            %}
            [cdata(time).tracked, newmax] = tracking(cdata(time).untracked,cdata(time-1).tracked,cdata(time-1).untracked,epst,newmax);
            %tracking2(cdata(time).untracked2,cdata(time-1).tracked,cdata(time-1).untracked2,epst,newmax)
        else
            newmax = max([cdata(time).untracked.id]);
            cdata(time).tracked = cdata(time).untracked;
            for ind = 1:size(cdata(time).tracked,2)
                cdata(time).tracked(ind).id = [cdata(time).tracked(ind).id,cdata(time).tracked(ind).id];
            end
>>>>>>> Stashed changes
		end
        cdata(time).eps = epst;
        cdata(time).minpts = minptst;
        cdata(time).newmax = newmax;
	end
	num_high_eps = numel([warnings(:).too_sparse]);
	if num_high_eps > 0 && p.verbose
		warning(append("Well: ", string(p.well), " - ", string(num_high_eps), ' timepoints had higher than expected epsilon values and may not cluster well'))
		disp("See the 'warnings' output for more information")
	end
end
<<<<<<< Updated upstream
function out = clustering(activeXY, eps, minpts)
    out = struct('XYCoord',0,'id',0);
    if isempty(activeXY); return; end
    clusters = dbscan(activeXY, eps, minpts);
    outliers = clusters < 0;
=======
function [out,out2] = clustering(activeXY, eps, minpts)
    out = struct('XYCoord',0,'id',0);
    if isempty(activeXY); return; end
    clusters = dbscan(activeXY, eps, minpts);
    outliers = clusters == -1;
>>>>>>> Stashed changes
    max_clust_id = max(clusters(~outliers));
    if isempty(max_clust_id); max_clust_id=0; end
    sum_outliers = sum(outliers);
    fake_clust = max_clust_id+1:max_clust_id+sum_outliers;
    clusters(outliers) = fake_clust';
<<<<<<< Updated upstream

=======
    out2.id = clusters;
    out2.XYCoord = activeXY;
>>>>>>> Stashed changes
    for cl = 1:max(clusters)
        pts = activeXY(clusters==cl,:);
        out(cl).XYCoord = pts; % XYCoord in that cluster
        out(cl).id = cl; %cluster identity
    end
end %clustering end
function [tracks,newmax] = tracking(sCurr,sPrev,bPrev,eps,newmax)	
    %%Unpack structs into arrays
    dCurr = unpack(sCurr); %XY and ID for clusters in curr
    dPrev = unpack(sPrev); %XY and ID for clusters in prev
    %%Check if the requisite data is present
<<<<<<< Updated upstream
    if numel(dPrev)<=3 
        dPrev = unpack(bPrev);
    end
=======
    
    if numel(dPrev)<3 
        dPrev = unpack(bPrev);
        dPrev(:,3) = dPrev(:,3) + newmax;
    end
    %{
>>>>>>> Stashed changes
    if numel(dCurr)<=3 
        tracks = []; %If no data, return empty
        return
    end
<<<<<<< Updated upstream
=======
    %}
>>>>>>> Stashed changes
    %%Search neighbors and reassign
    [idx,d] = knnsearch(dPrev(:,1:2),dCurr(:,1:2)); %Indices and distances of current's neighbors in previous
    isClose = d <= eps;
    dCurr(isClose,4)= dPrev(idx(isClose),3);
<<<<<<< Updated upstream
    for i = 1:max(dCurr(:,3))
        cluster = dCurr(:,3)==i;        %Logical map for points in currently evaluated cluster
        id = dCurr(cluster,3:4);        %cluster ids for those points
        n = numel(unique(id(id(:,2)>0,2)));         %Get the number of unique cluster ids for the current cluster >0 (non-ignored)
        if n > 1                            %more than one past cluster
            dCurr(cluster,5) = 0; 
=======
    uniqueClusterIDs = unique(dCurr(:,3))';
    for i = uniqueClusterIDs
        cluster = dCurr(:,3)==i;        %Logical map for points in currently evaluated cluster
        id = dCurr(cluster,3:4);        %cluster ids for those points
        n = numel(unique(id(id(:,2)>0,2)));         %Get the number of unique cluster ids for the current cluster >0 (non-ignored)

        % The fifth column is a flag for whether reassignment has occurred.
        % 0 if no reassignment, 1 if reassignment
        if n > 1                            %more than one past cluster
            %Do nothing?
            %Won't the cluster IDs conflict later on?
            dCurr(cluster,5) = 0; 
            
>>>>>>> Stashed changes
        elseif n == 1                       %One past cluster
            dCurr(cluster,4) = unique(id(id(:,2)>0,2));
            dCurr(cluster,5) = 0;
        elseif sum(id(:,2))==0              %No past clusters
            if max(max(dCurr(:,4)))+1 > newmax
                newmax = max(max(dCurr(:,4)))+1;
            end
<<<<<<< Updated upstream
            dCurr(cluster,4)= newmax;
=======
            dCurr(cluster,4) = newmax;
>>>>>>> Stashed changes
            dCurr(cluster,5) = 1;
        end
    end
    %%Format and output new cluster assignments
    for i = 1:max(dCurr(:,4))
        newclust = dCurr(:,4)==i;
        tracks(i).XYCoord = dCurr(newclust,1:2);  %#ok<AGROW> %Points that make up the cluster
        tracks(i).id = dCurr(newclust,3:4); %#ok<AGROW> %original cluster ID, new cluster ID
        if dCurr(newclust,5) > 0
            tracks(i).new = 1; %#ok<AGROW> %Flag it reassigned cluster is new
        else
            tracks(i).new = 0; %#ok<AGROW>
        end
    end
    %%Remove empty entries
    map = false(1,size(tracks,2));
    for i = 1:size(tracks,2)
        if isempty(tracks(i).XYCoord) || isempty(tracks(i).id)
            map(i) = 1;
        end
    end
    tracks(map) = [];
end %tracking end
function xy = unpack(d)
	%Helper function for tracking method - unpacks structs into arrays and
	%returns them
    xy = zeros(1,2);
    for i = 1:size(d,2) %Loop through struct elements
        for ii = 1:size(d(i).XYCoord) %Loop through elements of struct elements
            xy(end+1,1:2) = d(i).XYCoord(ii,:); %#ok<AGROW>
            if size(d(i).id,2) > 1
                xy(end,3) = d(i).id(ii,2);
            else
                xy(end,3) = d(i).id;
            end
        end
    end
    xy(1,:) = [];
end %unpack function end
<<<<<<< Updated upstream



%%Check tracking flag for new cluster assignments
=======
function tracks = tracking2(current,previousTracked,previousUntracked,eps,newmax)
    tracks = current.id;
    if numel(previousTracked) < 2 % FIXME This is a placeholder check that may not work
        previous = previousUntracked;
    else
        previous = previousTracked;
    end

    %%Search neighbors and reassign
    [idx,d] = knnsearch(previous.XYCoord(:,1:2),current.XYCoord(:,1:2)); %Indices and distances of current's neighbors in previous
    close = d <= eps;
    tracks(close) = previous.id(idx(close));
    currentIDsUnique = unique(current.id)';
    for i = currentIDsUnique
        clusterCurrentMap = current.id==i;        %Logical map for points in currently evaluated cluster
        id = tracks(clusterCurrentMap);       %cluster ids for those points
        n = numel(unique(id));         %Get the number of unique cluster ids for the current cluster >0 (non-ignored)
        if sum(id)==0
            if max(tracks)+1 > newmax
                newmax = max(max(dCurr(:,4)))+1;
            end
            tracks(clusterCurrentMap) = newmax;
        end
    end

end
>>>>>>> Stashed changes
