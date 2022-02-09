%% ARCOS Binarize
% Binarizes data using channel intensity values
%% Inputs
% * *ch* - |2D Matrix| - channel intensity data where rows are cells and
% columns are timepoints
% * _varargin_ - |option value pairs| - accepts optional inputs as option-value pairs.
%%% Optional Inputs
% * *percentile* - |int| - 1:99. The percentile of signal intensities above which will be 
% 'active' and below which will be ' inactive'. *Default: 75*
% 
%% Outputs
% * *bin* - |Matrix| - logical array where 1 is active and 0 is inactive.
% * *threshold* - |double| - The calculated threshold value used to determine 'active' or 'inactive' cells.
% 
%% Examples
% *Using default parameters*
%   bin = arcos_binarize(ch);
%
%   [bin,thr] = arcos_binarize(ch);
%
%   [~,thr] = arcos_binarize(ch);
%
% *Using optional parameters*
%
%   bin = arcos_binarize(ch,'percentile',80);
%
%   [bin,thr] = arcos_binarize(ch,'percentile',25);
%
%   [~,thr] = arcos_binarize(ch,'percentile',80);
%
function [bin,threshold] = arcos_binarize(ch,varargin)
p.percentile = 75;
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