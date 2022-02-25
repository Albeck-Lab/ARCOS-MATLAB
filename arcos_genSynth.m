function bin = arcos_genSynth(XCoord,YCoord,numSpreads,varargin)
    p.bin = [];
    p.eps = 1;
    p.lifetime = 2;
    p.seed = [];
    p.maxsize = 2.8;
    
    nin = length(varargin);     %Check for even number of add'l inputs
    if rem(nin,2) ~= 0
        warning('Additional inputs must be provided as option-value pairs');  
    end%Splits pairs to a structure
    for s = 1:2:nin
        p.(lower(varargin{s})) = varargin{s+1};   
    end
    
    if ~isempty(p.bin)
        binIn = p.bin;
    else
        binIn = boolean(zeros(size(XCoord,1), size(XCoord,2)));
    end
    if ~isempty(p.seed) %Use user-provided seed
    rng(p.seed);
    end
    for q = 1:numSpreads %Generate random indices for starting pts
        istartPts(q,1) = randi(size(XCoord,1));
    end
    binIn(istartPts,1) = 1;
    startPts = [XCoord(istartPts,1),YCoord(istartPts,1)];
    for time = 1:size(XCoord,2)
        pts = [XCoord(:,time),YCoord(:,time)]; %Get all xy coords at current time
        [idx,d] = rangesearch(pts,startPts,p.eps*time); %Get points within eps*time of start pts
        for i = 1:size(idx)
            binIn(idx{i}(d{i}<= p.eps*p.maxsize),time) = 1;
        end
        for row = 1:size(binIn,1)
            if sum(binIn(row,:))==p.lifetime
                binIn(row,time)=0;
            end
        end
    end
    bin = binIn;
end %EOF