%% ARCOS Binarize
% Binarizes data using channel intensity values
%
% <https://google.com Web Link>
%
% <<Link to image>>
%
% 
%% Inputs
% * *Input 1* - |Data type| - description of input
% * *Input 2* - |Data type| - description of input
% * _varargin_ - |option value pairs| - accepts optional inputs as option-value pairs.
%%% Optional Inputs
% * *Optional 1* - |Data type| , |Data type| - Description. *Default: default value*
% * *Optional 2* - |Data type| , |Data type| - Description. *Default: default value*
%% Outputs
% *output* - |Data type| - description of output
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