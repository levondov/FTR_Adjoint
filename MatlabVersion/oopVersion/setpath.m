function flag = setpath()
% Setup controls system directories
%
%
%

% add everything in UMER_CONTROL_2016 to matlab path
ROOT = fileparts(mfilename('fullpath'));
addpath(genpath(ROOT),'-begin');

end