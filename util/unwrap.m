%% unwrap.m
% *Summary:* Extract the numerical values from $s$ into the column vector $v$.
% The variable $sS can be of any type, including struct and cell array.
% Non-numerical elements are ignored. See also the reverse rewrap.m.
%
%    v = unwrap(s)
%
% *Input arguments:*
%
%   s     structure, cell, or numeric values
%
%
% *Output arguments:*
%
%   v     structure, cell, or numeric values
%
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-03-25

function v = unwrap(s)
%% Code

v = [];
if isnumeric(s)
  v = s(:);                        % numeric values are recast to column vector
elseif isstruct(s)
  v = unwrap(struct2cell(orderfields(s))); % alphabetize, conv to cell, recurse
elseif iscell(s)
  for i = 1:numel(s)             % cell array elements are handled sequentially
    v = [v; unwrap(s{i})];
  end
end