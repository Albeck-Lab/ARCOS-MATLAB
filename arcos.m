%% ARCOS
% Wrapper for arocs_core
% Input "bin" should be a logical array with indices corresponding to cells
% in XCoord and YCoord. 1 if cell state is "active", 0 if "inactive"
% Data for XCoord, YCoord and bin should be formatted such that each row is a
% cell and each column is a timepoint
function cdata = arcos(XCoord, YCoord, bin, varargin)
p.eps = [];
p.minpts = [];
%% Prepare additional inputs
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0
    warning('Additional inputs must be provided as option, value pairs');  
end%Splits pairs to a structure
for s = 1:2:nin
    p.(lower(varargin{s})) = varargin{s+1};   
end

%% Perform checks
assert(~isempty(XCoord), 'No x coordinate data detected');
assert(~isempty(YCoord), 'No y coordinate data detected');
assert(numel(XCoord) == numel(YCoord), 'XCoords and YCoords not equal')
if isempty(bin)
   warning('No binarized data supplied');
end


%% Loop through time
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
end
end %wrapper function end




%% Get minpts and epsilon
function [minpts, eps] = arcos_prep_dbscan(XCoord, YCoord)
%% Calculate minpts
% Default is 4 for 2D data according to Ester et al., 1996)
minpts = ndims(XCoord)*2;
%% Calculate eps
[~,d] = knnsearch([XCoord,YCoord], [XCoord, YCoord],'K', minpts+1); %k-nearest neighbors search
d = sort(d);
d = d(:,2:end);
max_d = max(d,[],2); %Biased toward greater distances as opposed to average of k-nearest
scaled = max_d * length(max_d)/max(max_d); %Scale the data
smoothed = smoothdata(scaled,'gaussian'); %Smooth it
slopes = gradient(smoothed); %Take first derivative
[~,ix]=min(abs(slopes-1)); %Get the index of avg_d for ideal eps (where slope of line tangent to that point is 1);
%eps = avg_d(ix);
eps = max_d(ix);
end %dbscan function end



%% ARCOS Core
% 
function events = arcos_core(activeXY, eps, minpts)
assert(eps>0, 'eps must be greater than 0') %These checks might need to happen earlier
assert(minpts>1, 'minpts must be greater than 1')
%assert(~isempty(activeXY), 'activeXY must contain at least one point')
if (isempty(activeXY))
    events = [];
    return
end
if minpts <=2 
    warning('minpts less than 3 may yield inaccurate results'); 
end
clusters = dbscan(activeXY, eps, minpts);
events = cell(max(clusters),1);
for cl = 1:max(clusters) %Treat clusters differently
    pts = activeXY(clusters==cl,:);
    if size(pts,1)>2
        [hull, area] = convhull(pts);
        cluster.pts = pts;
        cluster.hull = hull;
        cluster.area = area;
        events{cl, 1} = cluster;
    else
        continue
    end
end
end