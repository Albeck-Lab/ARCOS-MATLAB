%% ARCOS
% A limited feature set port of the R package "ARCOS"
%
% "Automated Recognition of Collective Signalling (ARCOS) is an R package
% to identify collective spatial events in time series data."
%
% <https://github.com/dmattek/ARCOS ARCOS GitHub>
%
% <<https://github.com/dmattek/ARCOS/raw/main/man/figures/README-ex1plotTS-1.png>>
%
% <<https://github.com/dmattek/ARCOS/raw/main/man/figures/README-ex1plotColl-1.png>> 
% 
%% Inputs
% * *filename* - |string| - full path to data file containing time series data.
% * _varargin_ - |string| - accepts optional inputs as option-value pairs.
%%% Optional Inputs
% * *epsilon* - |single| , |double| - Search radius used to find points in the
% dbscan clustering algorithm. Must be greater than zero. *Default: 0.3*
% * *minpts* - |single|, |double| - Minimum number of neighbors required for a
% core point, specified as a positive integer. Must be a positive integer greater than zero. 
% *Default: 1*
% * *time* - |1x2 array of ints| - Desired start and end time points. First
% element is start time, second element is end time. *Default: Earliest
% time to latest time index* 
% * *genplot* - |boolean| - Specify whether to plot the results. *Default:
% false*
%% Outputs
%
%% Examples
% *Using default parameters*
%
%   output = arcos('filepath\file.csv')
%
% *Using optional parameters, specifying epsilon of 0.25, minpts 2 and
% genplot to true.*
%
%   output = arcos('filepath\file.csv', 'epsilon', 0.25, 'minpts', 2, 'genplot', 'true')
%% See Also
% * <https://www.mathworks.com/help/stats/dbscan.html#mw_a35e3831-bd7a-437f-b128-889a3a444aa2
% Density-based spatial clustering of applications with noise (DBSCAN)>
% * <https://www.mathworks.com/help/matlab/ref/convhull.html Convex Hull>
%% To Do
% * Test for false-positive rate
% * Test for divergence events
% * Test for convergence events
% * Static x and y lims for plot


function dataloc = arcos(filename, varargin)
    p.epsilon = [];
    p.minpts = [];
    p.time = [];
    p.genplot = false;
    %%
    
    % Prepare additional inputs
    nin = length(varargin);     %Check for even number of add'l inputs
    if rem(nin,2) ~= 0
        warning('Additional inputs must be provided as option, value pairs');  
    end%Splits pairs to a structure
    for s = 1:2:nin
        p.(lower(varargin{s})) = varargin{s+1};   
    end
    %%
    
    % Set parameters/defaults
    if isempty(p.epsilon)
        epsilon = 0.3;
    else
        epsilon = p.epsilon;
    end
    if (isempty(p.minpts))
        minpts = 1;
    else
        minpts = p.minpts;
    end
    %%
    
    % Read Table and Perform Checks
    data = readmatrix(filename); %Read the file as a table
    assert(~isempty(data),'Table is empty');
    assert(epsilon>0,'epsilon must be greater than 0');
    assert(minpts>=1, 'minpts must be greater than 1');
    %%
    
    % Loop through time points
    if (isempty(p.time))
        p.time = [min(data(:,2)),max(data(:,2))];
    end
    %dataloc.s{iXY}.
    for time = p.time(1):p.time(2)
            row_time = data(:,2) == time; %1 if t column = time, else 0
            row_on = data(:,5) == 1 & row_time(:,1)==1; %1 if t column = time & flagged active, else 0
            xy_on = data(row_on==1, [3 4]);
            %%

            % Cluster points and generate convex hulls

            % Make spreads a cell array
            if size(xy_on,1)>0
                %clusters = xy_on(dbscan(xy_on,epsilon, minpts), :); %Get cluster indices for filtered data
                clusters = dbscan(xy_on,epsilon, minpts);
            else
                disp('Too few pts for cluster')
                continue
            end
            %%

            % Group clustered points together as separate arrays in output
            if size(clusters,1)>2 & ~ismember(clusters, -1)
                for cluster = 1:max(clusters)
                    pts = xy_on(clusters==cluster,:);
                    [hull, area] = convhull(pts);
                end
            else
                disp('Too few pts for convhull')
                continue %Too few points in cluster to form convex hull
            end
        end
    %%
    
    % Plot data
%     if p.genplot == true
%         %xlim([min(data(:,3))-1,max(data(:,3))+1]);
%         %ylim([min(data(:,4))-1,max(data(:,4))+1]);
%         
%         for ind = 1:size(out,1)
%             o = out{ind};
%             clf
%             hold on
%             plot(o.pts(:,1), o.pts(:,2), '-o', 'Color',[0.5,0.5,0.5], 'LineStyle', 'none')
%             
%             if ~isempty(o.active)
%                 plot(o.active(:,1),o.active(:,2),'.r', 'LineStyle', 'none');
%                 if ~isempty(o.spread)
%                     plot(o.active(o.spread, 1),o.active(o.spread, 2), '.r', 'LineStyle', 'none');
%                     if ~isempty(o.hull) %Is there a hull?
%                         plot(o.active(o.hull,1), o.active(o.hull,2));
%                     end  
%                 end
%             end
%             hold off
%             xlim([min(data(:,3))-1,max(data(:,3))+1]);
%             ylim([min(data(:,4))-1,max(data(:,4))+1]);
%             saveas(gcf,append(int2str(ind), '.png'))
%         end
%     end
end
% function dataloc = ARCOS(dataloc, inp)
% if ~isempty(inp.arcos)
%     NumXYs = numel(inp.arcos);
%     XYnums = inp.arcos;
% else
%     %get the xy names from the data
%     NumXYs = numel(dataloc.xys);
%     XYnums = extractBetween(dataloc.xys,'xy','.mat'); XYnums = str2double(XYnums);
% end
% 
% if ~isequal(NumXYs,numel(XYnums)); fprintf('Somehow the number of Xys you is not equal to the number of data XYs. \n'); end %FIX BREAK
% if isempty(inp.arcosepsilon)
%     epsilon = 0.15; % Should be calculated as a function of a point's distance to k-nearest neighbors
% else
%     epsilon = inp.arcosepsilon;
% end
% if (isempty(inp.arcosminpts))
%     minpts = 4;
% else
%     if inp.arcosminpts < 3
%         disp('Warning: Setting ARCOS minpts to a value less than three may not yield accurate results')
%     end
%     minpts = inp.arcosminpts;
% end
% assert(epsilon>0,'epsilon must be greater than 0');
% assert(minpts>=1, 'minpts must be greater than 1');
% PercentBar = floor(NumXYs/4); PercentBar = [1:PercentBar:PercentBar*4];
% Quarters = 0;
% 
% fprintf('Spreading event quantification will be performed on %s. \n', inp.arcos);
% 
% for iiXY = 1:NumXYs
%     iXY = XYnums(iiXY);
%     if any(iXY == PercentBar)
%        fprintf('Spatial cross correlations data processed: %d/4th done \n', Quarters); 
%        Quarters = Quarters +1;
%     end
%     if ~isempty(dataloc.d{iXY}) %check the XY isn't empty
%     if isfield(dataloc.d{iXY}, 'data') %check the XY has data
%     if isfield(dataloc.d{iXY}.data, 'XCoord') %check for XCoord
%     if ~isempty(dataloc.d{iXY}.data.XCoord) %check that there is actually data
%     if isfield(dataloc.d{iXY}.data, 'YCoord') %check for YCoord 
%         for time = 1:size(dataloc.d{iXY}.data.XCoord,2)
%             threshold = mean(prctile(dataloc.d{iXY}.data.EKAR,25)); %Cells with signal intensity > threshold are "active". Else "inactive"
%             activeCells = dataloc.d{iXY}.data.EKAR(:,time)>threshold;
%             activeCoord = [dataloc.d{iXY}.data.XCoord(activeCells==1,time), dataloc.d{iXY}.data.YCoord(activeCells==1,time)];
%             %%
% 
%             % Cluster active cells
%             if size(activeCoord,1)>0
%                 clusters = dbscan(activeCoord,epsilon, minpts);
%             else
%                 continue % Too few points to cluster
%             end
%             %%
% 
%             % Form convex hulls around clusters
%             if size(clusters,1)>2 & ~ismember(clusters, -1)
%                 for cluster = 1:max(clusters)
%                     pts = activeCoord(clusters==cluster,:);
%                     [hull, area] = convhull(pts);
%                 end
%             else
%                 continue %Too few points to form convex hull
%             end
%         end
%         dataloc.s{iXY}.data = out; % 
%     end %check for y coord
%     end %check for data existing
%     end %check for x coord
%     end %if data exists
%     end %if dataloc.d is not empty
%        
% end % iXY loop
% 
% fprintf('DatalocHandler has finished AROCS, saved as: %s. \n', dataloc.file.proc)
%dataloc.s
% Arrange data as fields in a matrix
% Want output as cell array per xy
% Take input 
% Open dataloc handler - 
% Implement into dataloc handler
% Output as struct - s - cell array with each cell being image xys with
% channel, cell array with every instance of spread, within that all the
% data associated
% Copy DatalocHandler line 665-680
% 700-720