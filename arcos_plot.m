classdef arcos_plot
    methods(Static)
        function plot(clust_by_time,clust_by_id,binarized_data,raw_data,xy,t,varargin)
            %% Optional parameters
            p.usebin = true; %Flag to use binarized data for visualizing active cells
            p.usebounds = true;
            p.tracked = true; %Flag to plot tracked clusters or untracked clusters
            p.outpath = pwd;
            p.gif = false; %Flag to output images to gif
            p.bin_real = []; %If using synthetic data + real data, supply real data here and it'll be colored differently
            %% Prep varargin structure
            nin = length(varargin);     %Check for even number of add'l inputs
            if rem(nin,2) ~= 0
                warning('Additional inputs must be provided as option, value pairs');  
            end%Splits pairs to a structure
            for s = 1:2:nin
                p.(lower(varargin{s})) = varargin{s+1};   
            end
            if p.usebounds == true
                if isempty(clust_by_id{1}(1).data(1).bounds) %This will bug out if there's no data in clust_by_id{1}
                    error('No bounds detected in data. Have you run analysis yet?')
                end
            end
            %% Set up environment
            mkdir(p.outpath); %Make outpath if it doesn't exist
            fh = figure(); %Figure handle
            set(fh,'WindowStyle','Normal') %Set figure window behavior
            set(fh,'Resize','off') %Lock figure dimensions
            warning("Warning: Do not close the figure until the process has finished");
            %% Loop through XY
            for iwell = 1:numel(xy)
				well = xy(iwell);
                mkdir(append(p.outpath,'\XY ',int2str(well))); %Make output directory
                pwell_data = clust_by_time{well}; %Processed data for the indexed well
                cwell_data = clust_by_id{well}; %Cluster-wise data for the indexed well
                rwell_data = raw_data{well}.data; %Raw data for the indexed well
                rXCoord = rwell_data.XCoord; %XCoords for all points, indexed well
                rYCoord = rwell_data.YCoord; %YCoords for all points, indexed well
                binw = binarized_data{well};
                %% Loop through time
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
                    if p.usebin == true %Plot active cells
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
                    if p.usebounds == true              %Plot bounds if specified
                        for cluster = 1:size(cwell_data,1)
                            for itime = 1:size(cwell_data(cluster).data,2)
                                if cwell_data(cluster).data(itime).time == time
                                    hold on
                                    bounds = cwell_data(cluster).data(itime).bounds;
                                    points = cwell_data(cluster).data(itime).XYCoord;
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
                %% Prepare GIF if specified
                if p.gif == true
                    arcos_plot.gif(p.outpath,'GIF');
                end
            end
        end 
        function gif(path,name)
            files = dir(append(path,'/*.png')); %Get list of filenames
            filename = append(path,'/',name,'.gif'); % Specify the output file name
            for f = 1:length(files) %Loop through files list
                im = imread(append(files(f).folder,'/',files(f).name)); %Read image into memory
                [A,map] = rgb2ind(im,256); %Convert from rgb
                if f == 1
                    imwrite(A,map,filename,'gif','LoopCount',Inf); %Initialize gif
                    bar = waitbar(f/length(files),append('Assembling GIF ',int2str(f),' of ', int2str(length(files)))); %Initialize progress bar
                else
                    imwrite(A,map,filename,'gif','WriteMode','append'); %Append gif
                    waitbar(f/length(files),bar,append('Assembling GIF ',int2str(f),' of ', int2str(length(files)))); %Update progress bar
                end
                hold off
            end
            close(bar); %Close progress bar
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
    end 
end


%tiledlayout - make figure, tell it dimensions of output or put in "flow"
%to dynamically change. Have to say "nexttile" even for first tile