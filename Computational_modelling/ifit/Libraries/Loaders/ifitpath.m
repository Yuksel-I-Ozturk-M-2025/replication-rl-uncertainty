function d = ifitpath
% ifitpath iFit library location
%
% Version: Aug. 22, 2017
% (c) E.Farhi, ILL. License: EUPL.

d = [ fileparts(which('iData/version')) filesep '..' filesep '..' filesep ];


