function [cdata,warnings] = arcos_core(XCoord,YCoord,bin,varargin)
	p.eps = [];
	p.minpts = [];
	nin = length(varargin);
	if rem(nin,2) ~= 0; warning('Additional inputs must be provided as option, value pairs'); end  %#ok<WNTAG>
	for s = 1:2:nin; p.(lower(varargin{s})) = varargin{s+1}; end
	runPrep = [isempty(p.eps), isempty(p.minpts)];
    %% Preallocate cdata, initialize newmax, preallocate warnings
    newmax = 0;
	cdata = struct('untracked',{},'tracked',{},'eps',{},'minpts',{},'newmax',{});
	warnings = struct('too_sparse',{});
	%% Time loop
    for time = 1:size(XCoord,2)
		%% Setup: Eps and Minpts
		if any(runPrep); [epst,minptst] = arcos_utils.prep_dbscan3(XCoord(:,time),YCoord(:,time),bin(:,time)); end
		if ~runPrep(1);  epst = p.eps;  end  %Override with user-provided eps
		if ~runPrep(2);  minptst = p.minpts; end  %Override with user-provided minpts
		%% Log warning if eps is abnormal
		if epst >150; warnings(time).too_sparse = "High Epsilon. Check data for low cell density";end % Warning if epsilon is abnormally high
		%% Setup and run Clustering
        activeXY = [XCoord(bin(:,time),time), YCoord(bin(:,time),time)]; %XY coordinates for active cells        
        cdata(time).untracked = clustering(activeXY, epst, minptst);
        %% Track if possible
		if time==1
            [cdata(time).tracked,newmax] = tracking(cdata(time).untracked,cdata(time).untracked,cdata(time).untracked,epst,cdata(time).untracked(end).id);
        elseif cdata(time).untracked(end).id==0
            cdata(time).tracked = [];
		else
            [cdata(time).tracked, newmax] = tracking(cdata(time).untracked,cdata(time-1).tracked,cdata(time-1).untracked,epst,newmax);
		end
        cdata(time).eps = epst;
        cdata(time).minpts = minptst;
        cdata(time).newmax = newmax;
    end
end
function out = clustering(activeXY, eps, minpts)
    out = struct('XYCoord',0,'id',0);
    if (isempty(activeXY))
        return %Return control to main function if no active points
    end
    clusters = dbscan(activeXY, eps, minpts);
    for cl = 1:max(clusters)
        pts = activeXY(clusters==cl,:);
        out(cl).XYCoord = pts; % XYCoord in that cluster
        out(cl).id = cl; %cluster identity
    end
end %clustering end
function [tracks,newmax] = tracking(sCurr,sPrev,bPrev,eps,newmax)
	%% Input explanations
	%sCurr = struct with current frame's untracked data
	%sPrev = struct with previous frame's tracked data
	%bPrev = struct with previous frame's untracked data for backup
	%eps = epsilon value for knnsearch
	%newmax = new maximum cluster id assignment 
	
    %% Unpack structs into arrays
    dCurr = unpack(sCurr); %XY and ID for clusters in curr
    dPrev = unpack(sPrev); %XY and ID for clusters in prev
    %% Check if the requisite data is present
    if numel(dPrev)<=3 
        dPrev = unpack(bPrev);
    end
    if numel(dCurr)<=3 
        tracks = []; %If no data, return empty
        return
    end
    %% Search neighbors and reassign
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
    %% Format and output new cluster assignments
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
    %% Remove empty entries
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



%% Check tracking flag for new cluster assignments