%% ARCOS Plot
%  Methods for plotting collective events
%% Plot
% Main plotting method for ARCOS.
%
% *Inputs*
%
% * XCoord
% * YCoord
% * cdata
% * t
%
% *Optional Inputs*
%
% * *tracked*
% * *usehull*
% * *outpath*
%
% *Examples*
%
%% GIF
% Method for creating animated gifs
%
% *Inputs*
%
% *Examples*
%

classdef arcos_plot
    methods(Static)
        function plot(XCoord, YCoord, cdata,t,varargin) %Wrapper for plotter
        p.tracked = true;
        p.usehull = true;
        p.outpath = pwd;
        p.bin = [];
        p.gif = false;
        nin = length(varargin);     %Check for even number of add'l inputs
        if rem(nin,2) ~= 0
            warning('Additional inputs must be provided as option, value pairs');  
        end%Splits pairs to a structure
        for s = 1:2:nin
            p.(lower(varargin{s})) = varargin{s+1};   
        end
        mkdir(p.outpath);
        fh = figure();
        fh.WindowState = 'maximized';
        for time = t(1):t(end)
            image = plotter(XCoord, YCoord, cdata, time,p.tracked, p.usehull, p.bin);
            path = append(p.outpath,append('/',int2str(time)));
            saveas(image,append(path, '.png'))
        end
        close gcf
        if p.gif == true
            arcos_plotter.gif(p.outpath,p.outpath);
        end
        end %EOF
        function gif(inpath, outpath)
            files = dir(append(inpath,'/*.png'));
            path = append(outpath,'/');
            filename = append(path,'gif.gif'); % Specify the output file name
            for f = 1:length(files)                
                %im = frame2im(getframe(image));
                im = imread(append(files(f).folder,'/',files(f).name));
                [A,map] = rgb2ind(im,256);
                if f == 1
                    imwrite(A,map,filename,'gif','LoopCount',Inf);
                else
                    imwrite(A,map,filename,'gif','WriteMode','append');
                end
                hold off
            end
        end %EOF
    end %End of methods
end %End of class
function image = plotter(XCoord, YCoord,cdata, time,tracked,usehull,bin)
            if tracked == true
                rw = 3;
            else
                rw = 1;
            end
            clf
            plot(XCoord(:,time),YCoord(:,time),'o', 'MarkerSize', 4, 'LineStyle', 'none' );
            hold on;
            axis square;
            if ~isempty(bin) %Plot active cells
                hold on;
                plot(XCoord(bin(:,time)==1,time),YCoord(bin(:,time)==1,time),'o','MarkerEdgeColor','k', 'MarkerSize', 7,'LineStyle','none');
                hold on;
            end
            xlim([0 inf])
            ylim([0 inf])
            for event = 1:size(cdata{rw,time},1) %Plot collective events
                if ~isempty(cdata{rw,time}{event,1})
                    xy = cell2mat(cdata{rw,time}{event,1});
                    plot(xy(:,1),xy(:,2), '.r','MarkerSize', 12, 'LineStyle', 'none')
                    if usehull == true              
                        hull = cell2mat(cdata{rw,time}{event,3});
                        plot(xy(hull,1),xy(hull,2), 'r');
                    end
                    hold on
                end
            end
            set(gca,'ydir','reverse') %Reverse Y axis (image origin is top left)
            image = gcf;
        end %EOF