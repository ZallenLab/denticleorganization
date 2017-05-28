

function savstr = SpiffyName(type, basename, varargin)
% makeafilename generates string to save genotype, etc into the data file
%   spiffyname(identifier, basename,varargin(1-4)

% fprintf('Number of arguments: %d\n',nargin)
% size(varargin,2)

narginchk(2,6)

if nargin == 3,
    savstr = sprintf('%3$s_%2$s.%1$s', type, basename, varargin{1});
elseif nargin == 4,
    savstr = sprintf('%3$s_%4$s_%2$s.%1$s', type, basename, varargin{1}, varargin{2});
elseif nargin == 5,
    savstr = sprintf('%3$s_%4$s_%5$s_%2$s.%1$s', type, basename, varargin{1}, varargin{2}, varargin{3});
elseif nargin == 6,
    savstr = sprintf('%3$s_%4$s_%5$s_%6$s_%2$s.%1$s', type, basename, varargin{1}, varargin{2}, varargin{3}, varargin{4});
    
elseif nargin > 5
    error('Too many arguments!!!')
elseif nargin < 2
    error('Not enough arguments!!!')
end

end

