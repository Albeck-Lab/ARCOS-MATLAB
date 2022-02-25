%% arcos_plot
% Description
%
% <https://google.com Web Link>
%
% <<Link to image>>
%
% 
%% Inputs
% * *XCoord* - |2D Matrix| - X coordinates. Rows are cells, columns are
% timepoints
% * *YCoord* - |2D Matrix| - Y coordinates. Rows are cells, columns are
% timepoints
% * *cdata* - |Cell array| - The output of the arcos script
% * t* - |2D array| - Beginning and ending timepoints. Ex [1,10]
% * _varargin_ - |option value pairs| - accepts optional inputs as option-value pairs.
%%% Optional Inputs
% * *save* - |boolean| - Set to true to save individual figures to file. *Default: false*
% * *gif* - |boolean| - Set to true to save the timeseries as an animated gif. *Default: false*
%% Examples
% *Using default parameters*
%
%   arcos_plot(XCoord,YCoord,cdata,[1,5]);
%
% *Using optional parameters*
%
%   arcos_plot(XCoord,YCoord,cdata,[1 10],'save',true);
%
function arcos_plot(XCoord, YCoord, cdata,t,varargin)
p.save = false;
p.gif = false;
p.tracked = true;
p.usehull = true;
p.outpath = pwd;
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
    image = plotter(XCoord, YCoord, cdata, time,p.tracked, p.usehull);
    if p.save==true
        path = append(p.outpath,append('/',int2str(time)));
        saveas(image,append(path, '.png'))
    end
    if p.gif == true
        gif(image,time, t(1))
    end
end
end %EOF
%% Plot
function image = plotter(XCoord, YCoord,cdata, time,tracked,usehull)
    
    if tracked == true
        rw = 3;
    else
        rw = 1;
    end
    clf
    %assert(~isempty(cdata{rw, time}), sprintf('No data for given time %d', time));
    plot(XCoord(:,time),YCoord(:,time),'o', 'MarkerSize', 4, 'LineStyle', 'none' );
    hold on;
    axis square;
    %xlim([min(XCoord,[],'all'),max(XCoord,[],'all')]);
    %ylim([min(XCoord,[],'all'),max(XCoord,[],'all')]);
    xlim([0 inf])
    ylim([0 inf])
    for event = 1:size(cdata{rw,time},1)
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
    
    
end
%% GIF Maker
function gif(image,idx, start)
    filename = 'GIF.gif'; % Specify the output file name
    im = frame2im(getframe(image));
    [A,map] = rgb2ind(im,256);
    if idx == start
        imwrite(A,map,filename,'gif','LoopCount',Inf);
    else
        imwrite(A,map,filename,'gif','WriteMode','append');
    end
    hold off
end