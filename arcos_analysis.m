classdef arcos_analysis
    methods(Static)
        function out = analyze(arcos_data,raw_data)
            %% Time Start
            run_start = now;
            for well = 1:size(arcos_data,2)
                %% Progress bar setup
                if well==1
                    bar = waitbar(well/size(arcos_data,2),append('Analyzing well ',int2str(well),' of ', int2str(size(arcos_data,2))));
                else
                    waitbar(well/size(arcos_data,2),bar,append('Analyzing well ',int2str(well),' of ', int2str(size(arcos_data,2))));
                end
                clusters = arcos_data{2,well};
                XCoord = raw_data{well}.data.XCoord;
                YCoord = raw_data{well}.data.YCoord;
                %% Loop through clusters
                for i = 1:size(clusters,1)
                    if size(clusters(i).data,2)>0
                        for ii = 1:size(clusters(i).data,2)
                            %% Boundary
                            [bounds,area] = boundary(clusters(i).data(ii).points);
                            clusters(i).data(ii).bounds = bounds;
                            clusters(i).data(ii).area = area;
                            %% Completion
                            xq = XCoord(:,clusters(i).data(ii).time); %Query points x values
                            yq = YCoord(:,clusters(i).data(ii).time); %Query points y values
                            xv = clusters(i).data(ii).points(clusters(i).data(ii).bounds,1); %Polygon x values
                            yv = clusters(i).data(ii).points(clusters(i).data(ii).bounds,2); %Polygon y values
                            [in,on] = inpolygon(xq,yq,xv,yv); %Logical map of points in and on the polygon (separately)
                            all = logical(in+on); %Logical map of points in and on the polygon (combined)
                            nactive = size(clusters(i).data(ii).points,1);
                            nall = sum(all);
                            clusters(i).data(ii).compl = nactive/nall; %Inf where no boundary can be drawn
                            %% Size
                            clusters(i).data(ii).numpts = size(clusters(i).data(ii).points,1);
                        end
                        clusters(i).t_start = clusters(i).data(1).time; %First appearance of cluster
                        clusters(i).t_end = clusters(i).data(end).time; %Last appearance of cluster
                        clusters(i).dur = 1+clusters(i).t_end - clusters(i).t_start; %Total duration the cluster existed
                        clusters(i).maxarea = max([clusters(i).data(:).area]); %Maximum area of cluster (should this include index of max?)
                        rocarea = diff([clusters(i).data(:).area]); %Rate of change in area
                        clusters(i).maxsize = max([clusters(i).data(:).numpts]);
                        rocsize = diff([clusters(i).data(:).numpts]); %Rate of change in size
                        for ii = 1:size(clusters(i).data,2)
                            if ii == 1
                                clusters(i).data(ii).rocarea = 0; %Rate of change since previous timepoint
                                clusters(i).data(ii).rocsize = 0; %Rate of change since previous timepoint
                            else
                                clusters(i).data(ii).rocarea = rocarea(ii-1); %Rate of change since previous timepoint
                                clusters(i).data(ii).rocsize = rocsize(ii-1); %Rate of change since previous timepoint
                            end
                        end
                    else
                        continue
                    end
                end
                arcos_data{2,well} = clusters;
            end
            out = arcos_data;
            close(bar) %Close progress bar
            run_end = now;
            elapsed = datestr(run_end - run_start,'HH:MM:SS FFF');
            disp(append('Elapsed time: ', elapsed));
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
% maxarea: maximum area of the boundary encompassing cluster
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
