classdef arcos_plot
    methods(Static)
        function plot(arcos_data,raw_data,xy,t,varargin)
            %% Optional parameters
            p.usebin = true; %Flag to use binarized data for visualizing active cells
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
            %% Set up environment
            mkdir(p.outpath); %Make outpath if it doesn't exist
            fh = figure(); %Figure handle
            set(fh,'WindowStyle','Normal') %Set figure window behavior
            set(fh,'Resize','off') %Lock figure dimensions
            if p.tracked == true %Select tracked or untracked dataset
                rw = 2; %Row 2 = tracked cluster data
            else
                rw = 1; %Row 1 = untracked cluster data
            end
            warning("Warning: Do not close the figure until the process has finished");
            %% Loop through XY
            for well = xy(1):xy(2)
                mkdir(append(p.outpath,'/','XY ',int2str(well))); %Make output directory
                pwell_data = arcos_data{1,well}; %Processed data for the indexed well
                rwell_data = raw_data{well}.data; %Raw data for the indexed well
                rXCoord = rwell_data.XCoord; %XCoords for all points, indexed well
                rYCoord = rwell_data.YCoord; %YCoords for all points, indexed well
                binw = arcos_data{3,well};
                
                %% Loop through time
                for time = t(1):t(end)
                    cdata = pwell_data{rw,time};
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
                        plot(XCoord(bin),YCoord(bin),'o','MarkerEdgeColor','r', 'MarkerSize', 8,'LineStyle','none'); %Synthetic data = red
                        hold on;
                    end
                    xlim([0 inf])
                    ylim([0 inf])
                    for event = 1:size(cdata,2) %Plot collective events
                        if cdata(event).points > 0
                            xy = cdata(event).points;
                            plot(xy(:,1),xy(:,2), 'o','MarkerEdgeColor','b','MarkerSize', 8, 'LineStyle', 'none') %clusters = blue
%                             if p.usebounds == true              %Plot bounds if specified
%                                 bounds = cdata(event).bounds;
%                                 plot(xy(bounds,1),xy(bounds,2), 'b');
%                             end
                            hold on
                        end
                    end
                    title(append("Well: ",int2str(well)," Time: ",int2str(time)));
                    legend('inactive','active real','active synth','cluster','bounds','Location','northeastoutside');
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