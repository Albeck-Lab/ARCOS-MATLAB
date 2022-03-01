%% ARCOS genSynth
% Synthetic spread generator for ARCOS
%
% 
%% Inputs
% * *XCoord* - |Data type| - description of input
% * *YCoord* - |Data type| - description of input

% * _varargin_ - |option value pairs| - accepts optional inputs as option-value pairs.
%%% Optional Inputs
% * *numSpreads* - |Integer| - the desired number of spreads. *Default: 5*
% * *bin* - |cell x time logical array| - Baseline binarized data. *Default: None*
% * *dist* - |double| , |int| - Distance the spread grows each frame. *Default: 1*
% * *lifetime* - |Integer| - How many frames a point will be active. *Default: 2*
% * *seed* - |Positive integer| - A seed value for the random number generator that selects start points. *Default: None*
% * *maxsize* - |Integer| , |Double| - The max size multiplier. Specifies how large a spread can be before it dissipates. This value is multiplied by *dist*. *Default: 2.8*
%% Outputs
% *bin* - |cells x time Logical array| - A logical array where rows are
% cells and columns are timepoints. Array size is determined by the
% dimensions of XCoord input.
%
%% Examples
% *Using default parameters*
%
%   output = function(input);
%
% *Using optional parameters*
%
%   output = function(input, 'optional 1', value);
%
%% See Also
% * Item 1
% * Item 2
%% To Do
% * Item 1
% * Item 2
% * Item 3
function bin = arcos_genSynth(XCoord,YCoord,interval,varargin)
    p.numspreads = 7;
    p.bin = [];
    p.dist = 1;
    p.lifetime = 2;
    p.seed = [];
    p.maxsize = 2.8;
    
    %%Setup
    nin = length(varargin);     %Check for even number of add'l inputs
    if rem(nin,2) ~= 0
        warning('Additional inputs must be provided as option-value pairs');  
    end%Splits pairs to a structure
    for s = 1:2:nin
        p.(lower(varargin{s})) = varargin{s+1};   
    end
    assert(sum(size(XCoord)==size(YCoord),'all')==2,'XCoord and YCoord must be equal size.');
    binIn = boolean(zeros(size(XCoord,1), size(XCoord,2)));
    if ~isempty(p.seed) %Use user-provided seed
    rng(p.seed);
    end
    if isempty(interval)
        interval(1)=2;
        interval(2)=size(XCoord,2);
    end
    for q = 1:p.numspreads %Generate random indices for starting pts
        istartPts(q,1) = randi(size(XCoord,1));
    end
    binIn(istartPts,1) = 1; %Set start pts to active
    startPts = [XCoord(istartPts,1),YCoord(istartPts,1)]; %Get coords for start pts
    for time = interval(1):interval(2)
        pts = [XCoord(:,time),YCoord(:,time)]; %Get all xy coords at current time
        [idx,d] = rangesearch(pts,startPts,p.dist*time); %Get points within dist*time of start pts
        for i = 1:size(idx)
            binIn(idx{i}(d{i}<= p.dist*p.maxsize),time) = 1; %Set points within dist*time of start pts to active
            for row = 1:size(binIn,1)
                if sum(binIn(row,:))>=p.lifetime
                binIn(row,time)=0; %Set cells to inactive if they've been active > lifetime
                end
            end
        end

    end
    if ~isempty(p.bin) %Use user-provided binary data
        bin = binIn+p.bin;
    else
        bin = binIn;
    end
end %EOF