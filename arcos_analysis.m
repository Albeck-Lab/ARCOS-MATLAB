classdef arcos_analysis
    methods(Static)
        function out = analyze(arcos_data,raw_data)
            %% Time Start
            run_start = now;
            %%Check for XYs with data and get their indicies
            xy = 1:numel(arcos_data);
            goodxys = ~cellfun(@isempty,arcos_data);% check to see if the input xys are good
            xy = xy(goodxys);
            %% Loop through XYs
            numXYs = numel(xy);
            for iwell = 1:numXYs
                well = xy(iwell);
                %% Progress bar setup
                if iwell==1
                    bar = waitbar(iwell/numXYs,append('Analyzing well ',int2str(iwell),' of ', int2str(numXYs)));
                else
                    waitbar(iwell/numXYs,bar,append('Analyzing well ',int2str(iwell),' of ', int2str(numXYs)));
                end
                %% Store clusters and coords
                clusters = arcos_data{well}; %fixed - change this for new input structure (single row cell array)
                XCoord = raw_data{well}.data.XCoord;
                YCoord = raw_data{well}.data.YCoord;
                %% Loop through clusters
                for i = 1:size(clusters,1)
                    if size(clusters(i).data,2)>0
                        for ii = 1:size(clusters(i).data,2)
                            %% Boundary
                            [bounds,area] = boundary(clusters(i).data(ii).XYCoord);
                            clusters(i).data(ii).bounds = bounds;
                            clusters(i).data(ii).area = area;
                            %% Completion
                            xq = XCoord(:,clusters(i).data(ii).time); %Query points x values (all x coords for time)
                            yq = YCoord(:,clusters(i).data(ii).time); %Query points y values (all y coords for time)
                            xv = clusters(i).data(ii).XYCoord(clusters(i).data(ii).bounds,1); %Polygon x values
                            yv = clusters(i).data(ii).XYCoord(clusters(i).data(ii).bounds,2); %Polygon y values
                            [in,on] = inpolygon(xq,yq,xv,yv); %Separate logical maps of active and inactive XYCoord in and on the polygon
                            all = logical(in+on); %Combined logical map of active and inactive XYCoord in and on the polygon
                            nactive = size(clusters(i).data(ii).XYCoord,1);
                            nall = sum(all);
                            clusters(i).data(ii).compl = nactive/nall; %Inf where no boundary can be drawn
                            %% Inactive points in bounds
                            if nactive/nall ~= 1 || nactive/nall ~= Inf %Skip if all cells within bounds are active or ratio is inf
                                xy_all = [XCoord(all,clusters(i).data(ii).time),YCoord(all,clusters(i).data(ii).time)]; %X and Y Coords for all points within and on the bounds
                                active = clusters(i).data(ii).XYCoord; %Active cells within the bounds
                                lmap = ~ismember(xy_all,active); %Logical map of inactive cells (2-dimensional logical array)
                                lmap = logical(lmap(:,1) + lmap(:,2)); %Collapsed 1D logical array
                                clusters(i).data(ii).inactive = [xy_all(lmap,1),xy_all(lmap,2)]; %X and Y coordinates of inactive cells
                            end
                            %% Size
                            clusters(i).data(ii).numpts = size(clusters(i).data(ii).XYCoord,1); %How many active points in the cluster
                        end
                        clusters(i).t_start = clusters(i).data(1).time; %First appearance of cluster
                        clusters(i).t_end = clusters(i).data(end).time; %Last appearance of cluster
                        clusters(i).dur = 1+clusters(i).t_end - clusters(i).t_start; %Total duration the cluster existed
                        clusters(i).maxarea = max([clusters(i).data(:).area]); %Maximum area of cluster (should this include index of max?)
                        rocarea = diff([clusters(i).data(:).area]); %Rate of change in area
                        clusters(i).maxcount = max([clusters(i).data(:).numpts]);
                        roccount = diff([clusters(i).data(:).numpts]); %Rate of change in count
                        for ii = 1:size(clusters(i).data,2)
                            if ii == 1
                                clusters(i).data(ii).rocarea = 0; %Rate of change since previous timepoint
                                clusters(i).data(ii).roccount = 0; %Rate of change since previous timepoint
                            else
                                clusters(i).data(ii).rocarea = rocarea(ii-1); %Rate of change since previous timepoint
                                clusters(i).data(ii).roccount = roccount(ii-1); %Rate of change since previous timepoint
                            end
                        end
                    else
                        continue
                    end
                end
                arcos_data{well} = clusters; %fixed - change this for new input structure (single row cell array)
            end
            out = arcos_data;
            %% Clean up and time reporting
            close(bar) %Close progress bar
            run_end = now;
            elapsed = datestr(run_end - run_start,'HH:MM:SS FFF');
            disp(append('Elapsed time: ', elapsed));
        end
        function out = filter(arcos_data,fts,varargin)
            % example input 
            % fts(1).t = 'Min'; fts(1).c = 'dur'; fts(1).p = 5;
            % this would filter and keep all the data where the spread
            % lasts at least 5 tps. 
            
            %% Default params
            p.xy = []; %which xys to process (if left empty it will do all)
           
            %% Process additional inputs
            nin = length(varargin);    
            if rem(nin,2) ~= 0; warning('Additional inputs must be provided as option, value pairs'); end  %Check for even number of add'l inputs
            for s = 1:2:nin; p.(lower(varargin{s})) = varargin{s+1}; end %Splits pairs to a structure

            %process either certain xys or all of them
            if ~isempty(p.xy)
                goodXYs = ~cellfun(@isempty,arcos_data(p.xy));
                XYnums = p.xy(goodXYs);
            else % figure out how many xys there are and which have data
                XYnums = 1:size(arcos_data,2);
                goodXYs = ~cellfun(@isempty,arcos_data);
                XYnums = XYnums(goodXYs);
            end      
            
            % List possible filter types
            fPars = {'max', 'min', 'custom'}; %add more? but idk what to add lol
            
            %set up possible fields for filtering, assert proper filter names / grab their function
            fNames = { ...
                'dur',                  []; ...
                't_start',              []; ...
                'maxarea',              []; ...
                'mean_roc_count',        '@(x)nanmean([x.data.roccount])';...
                'mean_roc_area',        '@(x)nanmean([x.data.rocarea])';...
                };    
            
            % Check that you gave filters
            nFts = numel(fts); % how many filters are there? %FIX to filter for and warn people of bad filters
                
            if nFts < 1; warning('Give at least one filter set... to filter the data...'); else            
            % Check there is a .c and .p for every .t
            if ~isstruct(fts) || any(cellfun(@isempty,struct2cell(fts)),'all'); warning('One of your filter parameters are empty, check them'); else

            % Find the filter names and their functions
            trueFts = arrayfun(@(x)find(~cellfun(@isempty,regexpi(x.c,fNames(:,1)))),fts,'UniformOutput',false)';
            badFts = cellfun(@isempty,trueFts); % toss bad filter names and warn them
            if any(badFts); sprintf('Failed to find proper filter info for %s. Disposing of it. \n',fts(cell2mat(trueFts(badFts))).c); fts = fts(~badFts); trueFts = trueFts(~badFts); end
            
            %Tell us once if we need to pull and run a function from the filter list above (fNames)
            functLog = ~arrayfun(@(x)isempty(fNames{x,2}),cell2mat(trueFts));
            %% Loop over all the data and see what you want to keep
            for iXY = 1:numel(XYnums)
                tXY = XYnums(iXY);
                keepThese = false(size(arcos_data{tXY},1),nFts); % set up the logicals for the filers
                for iFt = 1:nFts %hehe NFTs
                    % See if special filtering is needed and do it if so
                    if functLog(iFt)
                        aData = arrayfun(str2func(fNames{trueFts{iFt},2}),arcos_data{tXY});                      
                    else
                        aData = [arcos_data{tXY}.(fts(iFt).c)]';
                    end
                    
                    %fill in the logicals for keeping the data
                    switch lower(fts(iFt).t)
                        case 'max'
                            keepThese(:,iFt) = aData < fts(iFt).p;
                        case 'min'
                            keepThese(:,iFt) = aData > fts(iFt).p;
                        case 'custom'
                            keepThese(:,iFt) = any( fts(iFt).p(aData(fts(iFt).c)) , 2);
                    end
                end
                
                % now find the indicies that match all of your parameters
                keepThese = all(keepThese,2);
                
                % keep the data you actually want
                arcos_data{tXY} = arcos_data{tXY}(keepThese);
                
            end %xy loop
            out = arcos_data;
            
            end %if check for matching parameters 
            end %if check for at least one filter
        end %filtering function 
    end
end













































%adata format
% cid: cluster id
% data: time series data for cluster
% t_start: time point when cluster first appeared
% t_end: time point when the cluster disappeared
% dur: duration during which the cluster was tracked
% maxarea: maximum area of the boundary encompassing cluster
% maxcount: the maximum number of cells
% roc: rate of chance of maxarea
% compl : completion - active cells / all cells within the cluster























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




% %process either certain xys or all of them
% if ~isempty(inp.certainxys)
%     NumXYs = numel(inp.certainxys);
%     XYnums = inp.certainxys;
% else
%     %get the xy names from the data
%     NumXYs = numel(dataloc.xys);
%     XYnums = extractBetween(dataloc.xys,'xy','.mat'); change this 
%     XYnums = str2double(XYnums); 
% end