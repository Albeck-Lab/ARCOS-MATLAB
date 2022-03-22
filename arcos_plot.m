classdef arcos_plot
    methods(Static)
        function plot(XCoord, YCoord, cdata,t,varargin)
        p.tracked = true; %Flag to plot tracked clusters or untracked clusters
        p.usehull = true; %Flag to use convex hulls or ignore them
        p.outpath = pwd;
        
        p.gif = false; %Flag to output images to gif
        p.bin = []; %Binarized data
        p.bbin = []; %If using synthetic data + real data, supply real data here and it'll be colored differently
        nin = length(varargin);     %Check for even number of add'l inputs
        if rem(nin,2) ~= 0
            warning('Additional inputs must be provided as option, value pairs');  
        end%Splits pairs to a structure
        for s = 1:2:nin
            p.(lower(varargin{s})) = varargin{s+1};   
        end
        mkdir(p.outpath); %Make outpath if it doesn't exist
        fh = figure(); %Figure handle
        set(fh,'WindowStyle','Normal') %Set figure window behavior
        set(fh,'Resize','off') %Lock figure dimensions
        if p.tracked == true %Select tracked or untracked dataset
            rw = 2;
        else
            rw = 1;
        end
        warning("Warning: Do not close the figure until the process has finished");
        for time = t(1):t(end)
            image = plotter(XCoord(:,time), YCoord(:,time), cdata{rw,time},time, p.usehull, p.bin(:,time),p.bbin(:,time));
            saveas(image,append(p.outpath,'/',sprintf('%04d',time), '.png'))
        end
        close gcf
        if p.gif == true
            arcos_plot.gif(p.outpath,'GIF');
        end
        end 
        function gif(path,name)
            files = dir(append(path,'/*.png'));
            filename = append(path,'/',name,'.gif'); % Specify the output file name
            for f = 1:length(files)
                im = imread(append(files(f).folder,'/',files(f).name));
                [A,map] = rgb2ind(im,256);
                if f == 1
                    imwrite(A,map,filename,'gif','LoopCount',Inf);
                    bar = waitbar(f/length(files),append('Assembling GIF ',int2str(f),' of ', int2str(length(files))));
                else
                    imwrite(A,map,filename,'gif','WriteMode','append');
                    waitbar(f/length(files),bar,append('Assembling GIF ',int2str(f),' of ', int2str(length(files))));
                end
                hold off
            end
            close(bar);
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
function image = plotter(XCoord, YCoord,cdata,time,usehull,bin,bbin)
            clf %Clear current figure
            plot(XCoord,YCoord,'o', 'MarkerEdgeColor','k', 'MarkerSize', 3, 'LineStyle', 'none' ); %plot all cells in frame
            hold on;
            axis square;
            if ~isempty(bin) %Plot active cells
                hold on;
                if ~isempty(bbin)
                    bin = logical(bin-bbin);
                    plot(XCoord(bbin),YCoord(bbin),'o','MarkerEdgeColor','g', 'MarkerSize', 8,'LineStyle','none'); %Real data = green
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
                    if usehull == true              %Plot hulls if specified
                        hull = cdata(event).hull;
                        plot(xy(hull,1),xy(hull,2), 'b');
                    end
                    hold on
                end
            end
            title(time);
            legend('inactive','active real','active synth','cluster','hull','Location','northeastoutside');
            set(gca,'ydir','reverse') %Reverse Y axis (image origin is top left)
            image = gcf;
end 
%Merge gif and gif_sbs into one function that can accept an arbitrary
%number of paths