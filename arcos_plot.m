%% arcos_plot
% Description
%
% <https://google.com Web Link>
%
% <<Link to image>>
%
% 
%% Inputs
% * *XCoord* - |Data type| - description of input
% * *YCoord* - |Data type| - description of input
% * *cdata* - |Data type| - description of input
% * t* - |Data type| - description of input
% * _varargin_ - |option value pairs| - accepts optional inputs as option-value pairs.
%%% Optional Inputs
% * *save* - |boolean| - Description. *Default: default value*
% * *gif* - |boolean| - Description. *Default: default value*
%% Outputs
% *output* - |Data type| - description of output
%
%% Examples
% *Using default parameters*
%
%   arcos_plot(XCoord,YCoord,cdata,t);
%
% *Using optional parameters*
%
%   arcos_plot(XCoord,YCoord,cdata,t,'save',true);
%
%% See Also
% * Item 1
% * Item 2
%% To Do
% * Item 1
% * Item 2
% * Item 3
function arcos_plot(XCoord, YCoord, cdata,t,varargin)
p.save = false;
p.gif = false;
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0
    warning('Additional inputs must be provided as option, value pairs');  
end%Splits pairs to a structure
for s = 1:2:nin
    p.(lower(varargin{s})) = varargin{s+1};   
end

fh = figure();
fh.WindowState = 'maximized';
for time = t(1):t(end)
    image = plotter(XCoord, YCoord, cdata, time);
    if p.save==true
        saveas(image,append(int2str(time), '.png'))
    end
    if p.gif == true
        gif(image,time, t(1))
    end
end

end %EOF
%% Plot
function image = plotter(XCoord, YCoord,cdata, time)
    assert(~isempty(cdata{time}), sprintf('No data for given time %d', time));
    clf
    
    plot(XCoord(:,time),YCoord(:,time),'o', 'MarkerSize', 4, 'LineStyle', 'none' );
    hold on;
    axis square;
    xlim([min(XCoord,[],'all'),max(XCoord,[],'all')]);
    ylim([min(XCoord,[],'all'),max(XCoord,[],'all')]);
    for event = 1:size(cdata{time},1)
        if ~isempty(cdata{time}{event})
            xy = cdata{time}{event}.pts;
            hull = cdata{time}{event}.hull;
            plot(xy(:,1),xy(:,2), '.r','MarkerSize', 4, 'LineStyle', 'none')
            plot(xy(hull,1),xy(hull,2), 'r');
            hold on
        end
    end
    set(gca,'ydir','reverse') %Reverse Y axis (image origin is top left)
    
    image = gcf;
    
    
end
%% GIF Maker
function gif(image,idx, start)
    filename = 'testAnimated.gif'; % Specify the output file name
    im = frame2im(getframe(image));
    [A,map] = rgb2ind(im,256);
    if idx == start
        imwrite(A,map,filename,'gif','LoopCount',Inf);
    else
        imwrite(A,map,filename,'gif','WriteMode','append');
    end
    hold off
end