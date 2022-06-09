%% ARCOS Plot
% Class with methods for plotting cluster activity.
% 
%% Plot
% Plots collective activity in red and displays boundaries for visualization
% of clustering. Plots are saved as .png files in the current working
% directory or in the user-specified output path using 'outpath', 'path'
%
% *Inputs*
%
% * *clust_by_time - |Cell| - The clust_by_time cell array output of ARCOS
% * *clust_by_id* - |Cell| - The clust_by_id cell array output of ARCOS
% * *binarized_data* - |Cell| - The binarized data cell array output of
% ARCOS
% * *raw_data* - |Cell| - The raw_data cell array input of ARCOS
% * *xy* - |Array| - An array of XY/well indices to plot. Can be
% discontinuous Ex (1:5, 7:15)
% * *t* - |Array| - An array of timepoints to plot. Can be discontinuous Ex
% (1:5, 7:15)
% * _varargin_ - |Option-value pair| - Accepts additional inputs as
% option-value pairs. Ex 'usebin',true
%
% *Optional inputs*
%
% * *usebin* - |Boolean|, |Logical| - Flag to use binarized data for
% visualizing active unclustered cells. *Default value: true*
% * *usebounds* - |Boolean|, |Logical| - Flag to plot boundaries
% encompassing clusters. *Default value: true*
% * *tracked* - |Boolean|, |Logical| - Flag to plot tracked (true) and
% untracked (false) clusters. *Default value: true*
% * *outpath* - |String|, |Char| - Path to desired output directory.
% *Default value: pwd (current working directory)*
% * *gif* - |Boolean|, |Logical| - Flag to create animated gif from the
% saved plot images. *Default value: false*
% * *bin_real* - |Boolean|, |Logical| - If using synthetic data from
% arcos_utils layered onto real binary data (see arcos_utils.gensynth), supplying
% the real data here and it will be colored differently. *Default value:
% []*
%
%% GIF
% Method for loading images into an animated gif.
%
% *Inputs*
%
% * *path - |String|, |Char| - Path to the directory containing the images
% you wish to animate. Images will be loaded in alpha-numeric order. This
% will also serve as the output folder for the finished gif.
% * *name* - |String|, |Char| - Desired name for the gif file. Must contain
% legal characters only. Do not include the file extension in the name. Ex:
% 'my_gif'. 
%
%% GIF SBS (Side-by-side)
% Method for creating gifs for side-by-side comparison. Useful for
% validating clustering, tracking, filtering, etc. 
%
% Images in both directories _must_ be the same resolution, ex 640x480 == 640x480.
%  If dimensions are not equal it will break, ex 640x480 ~= 640x640.
%
% There _must_ be the same number of images in both directories. 
% 
% *Inputs*
%
% * *path1* - |String|, |Char| - Path to first directory (will appear on
% the left).
% * *path2* - |String|, |Char| - Path to the second directory (will appear
% on the right). 
% * *name* - |String|, |Char| - Desired name for the gif file. Must contain
% legal characters only. Do not include the file extension in the name. Ex:
% 'my_gif'.
%
%% Clusters
% This script plots all clusters by their cluster ID and overlays cluster
% data onto microscope images. 
% It includes a few frames prior to a cluster/spread's start so the user
% can view the environment just prior to a spread's occurrance.
% These frames are loaded into an animated gif file and saved to the user's
% current working directory. 
%
%The name of this script it misleading as it does *nothing* to identify
%apoptosis by itself. Rather, it provides the user with the data they need
%to identify apoptosis events.
%
% *Inputs*
%
% * *clust_by_id* - |Struct| - ARCOS output containing clusters organized
% by cluster ID. *Must be run through ARCOS Analysis "Analyze" first
% * *xy* - |Array| - Array of well indices (integers) to process. Can handle
% discontinuous indices. Ex: (1:5, 11:15)
% * *ch* - |Integer| - nd2 channel to view
% * *nd2path* - |String|, |Char| - Path to the nd2 file for the data.
classdef arcos_plot
    methods(Static)
        function plot(clust_by_time,clust_by_id,binarized_data,raw_data,xy,t,varargin)
            %%Optional parameters
            p.usebin = true; %Flag to use binarized data for visualizing active cells
            p.usebounds = true;
            p.tracked = true; %Flag to plot tracked clusters or untracked clusters
            p.outpath = pwd;
            p.gif = false; %Flag to output images to gif
            p.bin_real = []; %If using synthetic data + real data, supply real data here and it'll be colored differently
            %%Prep varargin structure
            nin = length(varargin);     %Check for even number of add'l inputs
            if rem(nin,2) ~= 0
                warning('Additional inputs must be provided as option, value pairs');  
            end%Splits pairs to a structure
            for s = 1:2:nin
                p.(lower(varargin{s})) = varargin{s+1};   
            end
            if p.usebounds == true
                if isempty(clust_by_id{xy(1)}(1).data(1).bounds) %This will bug out if there's no data in clust_by_id{1}
                    error('No bounds detected in data. Have you run analysis yet?')
                end
            end
            %%Set up environment
            mkdir(p.outpath); %Make outpath if it doesn't exist
            fh = figure(); %Figure handle
            set(fh,'WindowStyle','Normal') %Set figure window behavior
            set(fh,'Resize','off') %Lock figure dimensions
            warning("Do not close the figure until the process has finished");
			warning('off','MATLAB:legend:IgnoringExtraEntries');
            %%Loop through XY
            for iwell = 1:numel(xy)
				well = xy(iwell);
                mkdir(append(p.outpath,'\XY ',int2str(well))); %Make output directory
                pwell_data = clust_by_time{well}; %Processed data for the indexed well
                cwell_data = clust_by_id{well}; %Cluster-wise data for the indexed well
                rwell_data = raw_data{well}.data; %Raw data for the indexed well
                rXCoord = rwell_data.XCoord; %XCoords for all points, indexed well
                rYCoord = rwell_data.YCoord; %YCoords for all points, indexed well
                binw = binarized_data{well};
                %%Loop through time
                for itime = 1:numel(t)
					time = t(itime);
					if p.tracked == true
						cdata = pwell_data(time).tracked;
					else
						cdata = pwell_data(time).untracked;
					end
                    XCoord = rXCoord(:,time);
                    YCoord = rYCoord(:,time);
                    %image = plotter(XCoord(:,time), YCoord(:,time), cdata{rw,time},time, p.usebounds, p.bin(:,time),p.bin_real(:,time));
                    clf %Clear current figure
                    plot(XCoord,YCoord,'o', 'MarkerEdgeColor','k', 'MarkerSize', 3, 'LineStyle', 'none' ); %plot all cells in frame
                    hold on;
                    axis square;
                    if p.usebin %Plot active cells
                        bin = binw(:,time);
                        hold on;
                        if ~isempty(p.bin_real)
                            bin = logical(bin-p.bin_real);
                            plot(XCoord(p.bin_real),YCoord(p.bin_real),'o','MarkerEdgeColor','g', 'MarkerSize', 8,'LineStyle','none'); %Real data = green
                        end
                        plot(XCoord(bin),YCoord(bin),'o','MarkerEdgeColor','r', 'MarkerSize', 3,'LineStyle','none'); %Synthetic data = red
                        hold on;
                    end
                    xlim([0 inf])
                    ylim([0 inf])
                    for event = 1:size(cdata,2) %Plot collective events
                        if cdata(event).XYCoord > 0
                            xycoord = cdata(event).XYCoord;
                            plot(xycoord(:,1),xycoord(:,2), 'o','MarkerEdgeColor','b','MarkerSize', 3, 'LineStyle', 'none') %clusters = blue
                            hold on
                        end
                    end
                    if p.usebounds              %Plot bounds if specified
                        for cluster = 1:size(cwell_data,1)
                            for itime2 = 1:size(cwell_data(cluster).data,2)
                                if cwell_data(cluster).data(itime2).time == time
                                    hold on
                                    bounds = cwell_data(cluster).data(itime2).bounds;
                                    points = cwell_data(cluster).data(itime2).XYCoord;
                                    plot(points(bounds,1),points(bounds,2),'b');
                                end
                            end
                        end
                    end
                    title(append("Well: ",int2str(well)," Time: ",int2str(time)));
					if isempty(p.bin_real)
						legend('inactive','active','cluster','bounds','Location','northeastoutside');
					else
						legend('inactive','active_real','ative_synth','cluster','bounds','Location','northeastoutside');
					end
                    set(gca,'ydir','reverse') %Reverse Y axis (image origin is top left)
                    image = gcf;
                    saveas(image,append(p.outpath,'/','XY ',int2str(well),'/',sprintf('%04d',time), '.png'))
                end
                close gcf
                %%Prepare GIF if specified
                if p.gif == true
                    arcos_plot.gif(p.outpath,'GIF');
                end
            end
        end 
        function gif(path,name)
			disp('An empty figure will appear. Do not close it until processing completes')
            files = dir(append(path,'/*.png')); %Get list of filenames
            filename = append(path,'/',name,'.gif'); % Specify the output file name
			for f = 1:length(files) %Loop through files list
                im = imread(append(files(f).folder,'/',files(f).name)); %Read image into memory
                [A,map] = rgb2ind(im,256); %Convert from rgb
                if f == 1
                    imwrite(A,map,filename,'gif','LoopCount',Inf); %Initialize gif
                else
                    imwrite(A,map,filename,'gif','WriteMode','append'); %Append gif
                end
                hold off
			end
			disp('Processing complete')
			close(gcf)
        end 
        function gif_sbs(path1,path2,name)
            warning('Assuming equal numbers of files in both directories');
            files1 = dir(append(path1,'/*.png'));
            files2 = dir(append(path2,'/*.png'));
            filename = append(path1,'/',name,'.gif'); % Specify the output file name
            for f = 1:length(files1)
                im1 = imread(append(files1(f).folder,'/',files1(f).name));
                im2 = imread(append(files2(f).folder,'/',files2(f).name));
                im = cat(2,im1,im2);
                [A,map] = rgb2ind(im,256);
                if f == 1
                    imwrite(A,map,filename,'gif','LoopCount',Inf);
                    bar = waitbar(f/length(files1),append('Assembling GIF ',int2str(f),' of ', int2str(length(files1))));
                else
                    imwrite(A,map,filename,'gif','WriteMode','append');
                    waitbar(f/length(files1),bar,append('Assembling GIF ',int2str(f),' of ', int2str(length(files1))));
                end
                hold off
            end
            close(bar);
            close gcf
            disp(append('GIF Assembly complete. File saved in ',path1))
		end
		function clusters(clust_by_id,xy,ch,nd2path)
			handle=3; % The number of frames before and after a spread's start
			dao = ImageAccess(nd2path); %Data access object for the nd2 file
			if ~isfield(clust_by_id{xy(1)}, 't_start'); error("Insufficient data. Analysis must be run first."); end
			for ixy = 1:numel(xy)
				well = xy(ixy);
				spreads = clust_by_id{well};
				for spr = 1:size(spreads,1)
					spread = spreads(spr);
					if spread.t_start-handle>0
						first_frame = spread.t_start-handle;
					else
						first_frame = spread.t_start;
					end
					if spread.t_start+handle>dao.imax.t
						last_frame = dao.imax.t;
					else
						last_frame = spread.t_start+handle;
					end
					frames = first_frame:last_frame;
					for iframe = frames(1):frames(end)
						clf
						CurrFrame = dao.GetFrame([iframe,ch,well,1]);
						CurrFrame = uint16(CurrFrame);
            			adjust = stretchlim(CurrFrame);
            			CurrFrame = imadjust(CurrFrame,adjust,[]);
            			CurrFrame = double(CurrFrame);
            			CurrFrame = uint8(round(CurrFrame/256));
            			CurrFrame = repmat(CurrFrame, [1, 1, 3]);
            			imshow(CurrFrame);
						hold on
						plot(1,1)
						title(append('XY ',string(well),' Spread ',string(spr),' Time: ',string(iframe)))
						spreadData = [];
						for n = 1:size(spread.data,2)
							if iframe == spread.data(n).time
								spreadData = spread.data(n);
								break
							end
						end
						if ~isempty(spreadData)
							bounds = spreadData.bounds;
							points = spreadData.XYCoord;
							hold on
							plot(points(bounds,1),points(bounds,2),'g','LineWidth',1.5)
						end
						im = getframe(gcf);
						
						[A,map] = rgb2ind(im.cdata,256); %Convert from rgb
						filename = append('XY ',string(well),' Spread ',string(spr),'.gif');
						if iframe == frames(1)
                    		imwrite(A,map,filename,'gif','LoopCount',Inf); %Initialize gif
                		else
                    		imwrite(A,map,filename,'gif','WriteMode','append'); %Append gif
						end
						hold off
					end
				end %End spread loop
			end %End XY loop
		end
    end 
end


%tiledlayout - make figure, tell it dimensions of output or put in "flow"
%to dynamically change. Have to say "nexttile" even for first tile