%% ARCOS Binarize
% Binarizes data using channel intensity values
function [bin,threshold] = arcos_binarize(ch,varargin)
p.percentile = 25;
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0
    warning('Additional inputs must be provided as option, value pairs');  
end%Splits pairs to a structure
for s = 1:2:nin
    p.(lower(varargin{s})) = varargin{s+1};   
end

threshold = mean(prctile(ch,p.percentile));
bin = ch>threshold;
end