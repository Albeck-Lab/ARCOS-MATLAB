classdef arcos_analysis
    methods(Static)
        function adata = analyze(XCoord,YCoord,cdata)
            %%Load data into struct
            timerange = size(cdata,2);
            max_id = cdata{5,end};
            clusters = repmat(struct('cid',[],'data',struct('time',{},'points',{},'id',{},'hull',{},'area',{}),'t_start',[],'t_end',[],'dur',[], 'maxarea',[],'rocarea',[],'maxsize',[],'rocsize',[],'compl',[]),max_id,1); %set up struct of structs
            for time = 1:timerange
                dtp = cdata{2,time}; %dtp = data at timepoint
                for cluster = 1:size(dtp,2)
                    id = mode(dtp(cluster).id(:,2)); %May not need the mode part... compare timeit results with and without mode.
                    clusters(id).cid = id;
                    clusters(id).data(time).time = time;
                    clusters(id).data(time).points = dtp(cluster).points;
                    clusters(id).data(time).id = dtp(cluster).id;
                    clusters(id).data(time).hull = dtp(cluster).hull;
                    clusters(id).data(time).area = dtp(cluster).area;
                end
            end
            %Loop through substructs and remove empty entries
            for i = 1:size(clusters,1)  
                map = false(1,size(clusters(i).data,2));
                for ii = 1:size(clusters(i).data,2)
                    if isempty(clusters(i).data(ii).points) || isempty(clusters(i).data(ii).id)
                        map(ii) = 1;
                    end
                end
                clusters(i).data(map) = [];
            end
            %Get start and end timepoints for cluster
            for i = 1:size(clusters,1)
                if size(clusters(i).data,2)>0
                    clusters(i).t_start = clusters(i).data(1).time;
                    clusters(i).t_end = clusters(i).data(end).time;
                    clusters(i).dur = 1+clusters(i).t_end - clusters(i).t_start; 
                    areas = [clusters(i).data(:).area];
                    clusters(i).maxarea = max(areas);
                    clusters(i).rocarea = diff(areas);
                    sizes = [];
                    ratio = [];
                    for ii = 1:size(clusters(i).data,2)
                        sizes(ii) = size(clusters(i).data(ii).points,1); %#ok<AGROW>
                        xq = XCoord(:,clusters(i).data(ii).time); %Query points x values
                        yq = YCoord(:,clusters(i).data(ii).time); %Query points y values
                        xv = clusters(i).data(ii).points(clusters(i).data(ii).hull,1); %Polygon x values
                        yv = clusters(i).data(ii).points(clusters(i).data(ii).hull,2); %Polygon y values
                        [in,on] = inpolygon(xq,yq,xv,yv); %Logical map of points in and on the polygon (separately)
                        all = logical(in+on); %Logical map of points in and on the polygon (combined)
                        nactive = size(clusters(i).data(ii).points,1);
                        nall = sum(all);
                        ratio(ii,1) = nactive/nall; %#ok<AGROW>
                        %Ratio of active:inactive 
                        %Inactive = all-active
                    end
                    clusters(i).maxsize = max(sizes);
                    clusters(i).rocsize = diff(sizes);
                    clusters(i).compl = ratio;
                else
                    continue
                end
            end
            adata = clusters;
        end
        function fdata = filter(adata,varargin)
            p.filt_dur = []; %Filter data by duration. Specify upper and lower duration thresholds
            p.filt_roc = []; %Filter data by rate of change. Specify upper and lower roc thresholds
            nin = length(varargin);     %Check for even number of add'l inputs
            if rem(nin,2) ~= 0
                warning('Additional inputs must be provided as option, value pairs');  
            end%Splits pairs to a structure
            for s = 1:2:nin
                p.(lower(varargin{s})) = varargin{s+1};   
            end
            for i = 1:size(adata,1)
                if adata(i).dur < p.filt_dur(1) || adata(i).dur > p.filt_dur(2)
                    adata(i) = []; %Removing an element of a struct doesn't leave gaps, so the index gets thrown off
                end
            end
            fdata = adata;
        end
    end
end
%adata format
% cid: cluster id
% data: time series data for cluster
% t_start: time point when cluster first appeared
% t_end: time point when the cluster disappeared
% dur: duration during which the cluster was tracked
% maxarea: maximum area of the convex hull encompassing cluster
% maxsize: the maximum number of cells
% roc: rate of chance of maxarea
% compl : completion - active cells / all cells within the cluster

%Output inactive cells from inpolygon into data field

% Donut versus Disk of activity
% 
% Measurement of Merge/split scenarios? - Dealing with them?
%Use negative roc and duration of activity to filter data


%Prop - Velocity with units - um/sec
%How to get scale?
%Have arcos_analysis take extra optional params - um/pixel
%Wrapper can pull layer of dataset that includes objective mag, other
%relevant params for unit calculation

%How to take mean across single event? Point-to-point? 


%When filtering, save filter params somewhere

%Define properties for class - maintains filter dataset and filter props
%
