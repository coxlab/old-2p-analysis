function result = calc_gamma(value,coef,varargin)
% CALC_GAMMA returns the gamma corrected values for VALUE
% using the specified coefficient COEF. Note that the values
% will be normalized to 1.
% OPTION: specify MAXVAL to change point
% where the curve reaches a height of 1
%
% EXAMPLE:
% calc_gamma(127,.45)
% CALC_GAMMA also deals with vectors/matrices
% plot(calc_gamma(0:255,.45,255))
% calc_gamma(repmat((0:255)',1,3),.45,255)

if nargin==1
    error('At least 2 input arguments should be specified...')
elseif nargin==2
    maxVal=255;
elseif nargin==3
    maxVal=varargin{1};
elseif nargin>3
    error('At most 3 input arguments may be specified...')
end

if coef>1
    %disp('Warning: gamma-value exceeds 1. Typically this should be lower than 1...')
end

result = value.^(coef) ./ (maxVal.^(coef));
