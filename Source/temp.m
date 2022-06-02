% MakeMoviesND(nd2path, videoXY, varargin)
% Make greyscale movies with data (or not) overtop of single or multiple
% XYs (in the same movie)
% 
% Required input parameters:
%
% nd2path - full path and file of the ND2 file of the video you want to use.
% videoXY - the xy(s) in your movie that you want to plot
%
% Important optional parameters:
% moviechan - the channel (number) to use to make the greyscale movie
%   Default = 1; (the first channel from the movie)
%
% overlaytype - blank (no overlay - default), spreads
% under construction - livecellbars, livecellplots
%
% data - if you want to overlay data on the cells, include that data here.
%   data must be a cell array of the relevent data for the overlay type
%   desired (see below for overlay types)
%
% chans - a string or cell array of the channels (as they are called in
%       your data) that you want on the cells
%
% chanminprctl - Min percentile of data to use as a "bottom" for the bar
%           default = 2
% chanmaxprctl - Max percentile of data to use as a "bottom" for the bar
%           default = 97
%
% Optional parameters:
% fps - frames per second you want for your video (default = 3 fps)
%
% spike - a time point where you added drugs, text will be added to the top
%   left corner at that time
%
% xystogether - put all the xys you gave into one movie? or make each one
%   seperate (default = true)
%
% tsamp - how many minutes does it take for a movie loop? (default = 6 (min))
%
% quality - output video quality desired (default = 95 (out of 100))
%
% twin - time frame you want the video from 
%   (default = all frames, input as [tStart, tFinish])
%
% movieoutfolder - where do you want the movie saved? (default - will save
%   the movie in a new folder within the same folder as your ND2 file)
%
% Features under construction: 
% filetype - what type of file you want output (default = .avi)
%
% N. DeCuzzi, Albeck Lab, UC Davis. 2020-2022


function arcos_id_apoptosis(nd2path, videoXY, varargin)
%% Normal video parameters
p.fps = 3; p.filetype = 'AVI';p.tsamp = 6; p.quality = 95;
p.moviechan = 1; p.overlaytype = 'blankmovies'; p.twin = [];
p.spike = []; p.data = []; p.spreadsbytime = []; p.spreadsbyid = []; p.bindata = [];
p.markersize = 10; p.pmd = [];
p.xystogether = true; p.movieoutfolder = []; 
p.chans = {}; %chans for live cell data
p.debug = false;

%% Bar parameters 
p.chanminprctl = 1; p.chanmaxprctl = 99;
p.UpDown = false; %for later to add updown bars instead FIX
p.BarColors = [0.0588 1 1; 1 0 0] ;
p.BarHeight = 10; % bar pixel height
p.BarGap = 5; % gap between bars
p.BarWidth = 50;  %bar pixel width

%% Check input varargin parameters
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end%Splits pairs to a structure
for s = 1:2:nin;   p.(lower(varargin{s})) = varargin{s+1};   end

%% Assert Channels as cell array and save your old path
if isstring(p.chans)
    p.chans{1} = p.chans;
end

%% Video and Data Info
[VideoPath, VideoFile] = fileparts(nd2path); %get the video folder
% make a new folder for the output movies (or find it)
if isempty(p.movieoutfolder)
    p.movieoutfolder = [VideoPath, '\Movies'];
end
if ~exist(p.movieoutfolder, 'dir')
    mkdir(p.movieoutfolder)
end

%% Get Video Object
dao = ImageAccess(nd2path);

%% add platemap data
if ~isempty(p.pmd)
    p.pmd = ct_maptx(p.pmd);
end

p.overlaytype = lower(strrep(p.overlaytype,' ',''));

% pass everything along
if p.xystogether; videoXY = {videoXY}; end % trick the xy loop into thinking there is only 1 xy

for iiXY = 1:numel(videoXY)
    if p.xystogether; ttXY = videoXY{iiXY}; else; ttXY = videoXY(iiXY); end % now do a backflip
    %% Set up the video to make
    newVideoName = [VideoFile ' XY ', strrep(strrep(num2str(ttXY),'   ',' '),'  ', ' '), ' chan ', num2str(p.moviechan),'.avi'];
    newVideoFile = fullfile(p.movieoutfolder, newVideoName);
    OverlayVideo = VideoWriter(newVideoFile); %#ok<TNMLP> 
    OverlayVideo.Quality = p.quality;
    OverlayVideo.FrameRate = p.fps;
    open(OverlayVideo);
    
    % Assert a time window (if desired)
    if isempty(p.twin)
        p.twin = [1, dao.imax.t];
    end
    
    % Get the frames
    VideoFrames = p.twin(1):p.twin(2);
    
    % Set up the figure
    figgy = tiledlayout('flow','Padding','tight','TileSpacing','tight');
    %% Loop through every time frame
    for Frame = 1:numel(VideoFrames)
        iFrame = VideoFrames(Frame);
        if p.debug; fprintf('Processing frame: %d\n', iFrame); end 
        % Loop throught the XYs
        for iXY = 1:numel(ttXY)
            tHands = nexttile(figgy,iXY); % get the tile
            tXY = ttXY(iXY); % get the xy
            if p.debug; printf('Processing XY: %d\n', tXY); end 
            %get the image and adjust it
            CurrFrame = dao.GetFrame([VideoFrames(iFrame), p.moviechan, tXY, 1]);
            CurrFrame = uint16(CurrFrame);
            adjust = stretchlim(CurrFrame);
            CurrFrame = imadjust(CurrFrame,adjust,[]);
            CurrFrame = double(CurrFrame);
            CurrFrame = uint8(round(CurrFrame/256));
            CurrFrame = repmat(CurrFrame, [1, 1, 3]);
            imshow(CurrFrame)

            % Call the right imageappending function
            switch p.overlaytype 
                case {'blank','blankmovies'}
                    % Do nothing
                case {'spreads','spread','arcos'}
                    MovieSpreads(tHands, iFrame, tXY, p)
                case {'detailedspreads','detailedspread','detailedarcos','spreadsbytime'}
                    InDepthSpreadMovies(tHands, iFrame, tXY, p)
                case {'livecellbars', 'bars'}
                    % MovieBars(dao, tXY, VideoFile, p)
                case 'livecellplots'
                    % make one for single cell plots next to the movies (dao, videoXY, VideoFile, p)
            end

            if ~isempty(p.pmd)
                %add a label
                whichLabel = structfun(@(x)any(ismember(x.xy,tXY)),p.pmd);
                fldz = fieldnames(p.pmd);
                fldz = fldz{whichLabel};
                emptyTpData = cellfun(@isempty,{p.pmd.(fldz).tx.time});
                [p.pmd.(fldz).tx.time] = deal(0);
                whichLabel2 = [p.pmd.(fldz).tx.time]' >= iFrame;
                
            else
                xylab = ['XY ', num2str(tXY)];
            end

            title([xylab,' t = ', num2str(iFrame*p.tsamp), ' min']) %Fix to add titles for treatments
        end
        set(figgy,'Position',[0,0,1920,1080])
        thisFrame = getframe(figgy);
        writeVideo(OverlayVideo,thisFrame)
    end
    
    close(OverlayVideo)
    close(gcf)
end %first xy loop
end %MakeMovies ND function

%% Spread Movies using clust_by_ID
function MovieSpreads(tHands, iFrame, tXY, p)
    hold on 
    spreadData = [p.spreadsbyid{tXY}.data]; %collect the spread data
    if ~isempty(spreadData)
        tData = [spreadData.time] == iFrame; %get the right time
        spreadData = spreadData(tData);
        for iSpread = 1:numel(spreadData)
            bounds = spreadData(iSpread).bounds;
            points = spreadData(iSpread).XYCoord;
            plot(tHands,points(bounds,1),points(bounds,2),'g','LineWidth',1.5)
        end
    end % spread check
    hold off
end % Spread appender

%% Spread Movies using clust_by_time
function InDepthSpreadMovies(tHands, iFrame, tXY, p)
    hold on
    binData = p.bindata{tXY}(:,iFrame); % get the logical for the xy coords of positive cells    spreadData = spreadData(tData);
    XYdata = [p.data{tXY}.data.XCoord(:,iFrame), p.data{tXY}.data.YCoord(:,iFrame)]; %get the xy coords
    tXYdata = XYdata(binData,:);
    if ~isempty(tXYdata)
        plot(tHands,tXYdata(:,1),tXYdata(:,2),'o','MarkerEdgeColor','r','MarkerSize',p.markersize,'LineWidth',1.5,'LineStyle','none') %clusters = red, thanks Daniel!
    end %xy data check

    spreadData = [p.spreadsbyid{tXY}.data]; %collect the spread data
    if ~isempty(spreadData)
        tData = [spreadData.time] == iFrame; %get the right time
        spreadData = spreadData(tData);
        for iSpread = 1:numel(spreadData)
            bounds = spreadData(iSpread).bounds;
            points = spreadData(iSpread).XYCoord;
            plot(tHands,points(bounds,1),points(bounds,2),'g','LineWidth',1.5)
        end
    end % spread check

    hold off
end % Detailed spread movies

%% Movie Bars -- FIX THIS
function MovieBars(dao, videoXY, VideoFile, p)
    %% Ask what cells you want to visualize
    pickTheCells = iman_getframe(dao, [p.twin(1), p.moviechan, tXY, 1]);
    pickTheCells = uint16(pickTheCells);
    adjust = stretchlim(pickTheCells);
    pickTheCells = imadjust(pickTheCells,adjust,[]);
    
    cellRange = size((p.d{tXY}.data.XCoord),1);
    cellNumbersText = cell(cellRange,1);
    
    for iCell = 1:cellRange
        cellNumbersText{iCell} = ['\color{cyan} ' num2str(iCell)];
    end        
    
    fig1 = figure;
    imshow(pickTheCells);
    text(p.d{tXY}.data.XCoord(:,p.twin(1)), p.d{tXY}.data.YCoord(:,p.twin(1)), cellNumbersText);
    
    cellsToPlot = input('What cells would you like to see? (ex: [1 2 3 4 5])\n');
    
    if isempty(cellsToPlot)
       cellsToPlot = 1:5;
    end
    
    numCellsToPlot = numel(cellsToPlot);
    clear PickTheCells;
    close(fig1)
    
    %% Set up the bars
    
    % Set up the background Bar (do not change)
    WhiteBar = uint8(ones(p.BarHeight, p.BarWidth, 3)*255);
    
    for iChan = 1:NumChans
        BarName = p.barchans{iChan};
        LegendBar.(sprintf('%s',BarName)) = zeros(p.BarHeight, p.BarWidth, 3);
        LegendBar.(sprintf('%s',BarName))(:,:,1) = uint8(round(p.BarColors(iChan, 1)*255));
        LegendBar.(sprintf('%s',BarName))(:,:,2) = uint8(round(p.BarColors(iChan, 2)*255));
        LegendBar.(sprintf('%s',BarName))(:,:,3) = uint8(round(p.BarColors(iChan, 3)*255));  
    end
    
    % Figure out how many bars you need (number of channels to plot) and set
    % them up to be used 
    
    % Assign values to use as manipulators for the cell's position
    BarYMod = NaN(NumChans, 1);
    UpDownSplit = ceil(NumChans/2);
    
    if NumChans == 1
        BarYMod = (p.BarHeight/2);
    
    elseif rem(NumChans, 2) == 0 %even number of bars !! DETERMINE TOP OF BOX Y LOCATION
        
        for iChan = 1:NumChans
            if iChan <= UpDownSplit
               BarYMod(iChan) = ((p.BarHeight + (p.BarGap/2)) + ((p.BarHeight + p.BarGap) * (UpDownSplit - iChan)));         
               
            else
               BarYMod(iChan) = (((p.BarHeight + p.BarGap) * (UpDownSplit - (iChan - 1))) - (p.BarGap/2)); % y pos modifier for top of bar
            end
        end
                
        
    elseif rem(NumChans, 2) == 1 %odd number of bars 
        
        for iChan = 1:NumChans %for all bars
            BarYMod(iChan) =  (((p.BarHeight + p.BarGap) * (UpDownSplit - iChan)) + (p.BarHeight/2)); %y pos modifier for top of bar 
        end
    end
    
    
    % X Bar modifier
    BarXMod = -(p.BarWidth/2);
    
    
    %% Set up the video to make
    
    OverlayVideo = VideoWriter(newVideoFile);
    OverlayVideo.Quality = p.quality;
    OverlayVideo.FrameRate = p.fps;
    open(OverlayVideo);
    
    % Now get an image to overlay all this daaaattttaaaa
    if ~isempty(p.twin)
        VideoFrames = p.twin;
    else
        VideoFrames = 1:TimePts;
    end
    
    for iFrame = 1:numel(VideoFrames)
        
        CurrFrame = iman_getframe(dao, [VideoFrames(iFrame), VideoChan, VideoXY, 1]);
        CurrFrame = uint16(CurrFrame);
        fprintf('Processing frame: %d\n', iFrame);
        adjust = stretchlim(CurrFrame);
        CurrFrame = imadjust(CurrFrame,adjust,[]);
        CurrFrame = double(CurrFrame);
        CurrFrame = uint8(round(CurrFrame/256));
        CurrFrame = repmat(CurrFrame, [1, 1, 3]);
        
        %imshow(CurrFrame)%test feature
    
        for iCell = 1:NumCellsToPlot %Take each cell you want to plot
            
            CurrXCoord = round(CellPositionData.XCoords(iCell, iFrame)); %get the cell's coordinates at that frame
            CurrYCoord = round(CellPositionData.YCoords(iCell, iFrame));
            
            if isnan(CurrXCoord) == true
            else
            %make the bars
            for iChan = 1:NumChans
                
                CurrBarYPosTop = (CurrYCoord + BarYMod(iChan)); 
                CurrBarYPosTop = floor(CurrBarYPosTop);
                CurrBarYPosBottom = ((CurrBarYPosTop - (p.BarHeight-1))); 
                CurrBarYPosBottom = floor(CurrBarYPosBottom);
                
                CurrBarXPosLeft = (CurrXCoord + (BarXMod+1)); 
                CurrBarXPosRight = (CurrXCoord - (BarXMod)); 
                
                CurrFrame(CurrBarYPosBottom:CurrBarYPosTop, CurrBarXPosLeft:CurrBarXPosRight,:) = WhiteBar(:,:,:); %make the white background bar
                
                DataBarWidth = floor((RelDataForVideo(iChan).Data(iCell, iFrame)*p.BarWidth));
                
                DataBar = zeros((p.BarHeight-2), DataBarWidth, 3);
                DataBar(:,:,1) = (p.BarColors(iChan, 1)*255);
                DataBar(:,:,2) = (p.BarColors(iChan, 2)*255);
                DataBar(:,:,3) = (p.BarColors(iChan, 3)*255);
                DataBar = uint8(DataBar); 
                
                CurrFrame((CurrBarYPosBottom+1):(CurrBarYPosTop-1), CurrBarXPosLeft:(CurrBarXPosLeft+DataBarWidth-1),:) = DataBar(:,:,:); %add in the data!     
            end
            
            end
        %picc = imshow(CurrFrame);    
        end
        
        
        
        %FIX: Add color bar legends
    %     for iChan = 1:NumChans
    %     BarName = p.BarChans{iChan};
    %     CurrLedgYPosTop = (100 + BarYMod(iChan)); 
    %     CurrLedgYPosTop = floor(CurrLedgYPosTop);
    %     CurrLedgYPosBottom = ((CurrLedgYPosTop - (p.BarHeight-1)));
    %     LegendBar.(sprintf('%s',BarName)) 
    %     CurrFrame((100):(CurrLedgYPosTop), 100:(100 + DataBarWidth-1),:)
    %     end
    
        
        writeVideo(OverlayVideo, CurrFrame);
     
    end
    
    close(OverlayVideo)
       
    fprintf('Your video is now finished processing, the video can be found\n in the folder: %s \n and is named: %s \n', DataPath, newVideoName); 
    
    
    
    path(oldpath)
end

