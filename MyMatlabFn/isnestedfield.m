function yn = isnestedfield(s, fields)
%ISNESTEDFIELD True if nested fields are in a structure.
%   ISNESTEDFIELD(S,FIELDs) returns true if the string FIELDs is containes
%   nested fields of S in the format of a single string of fieldnames 
%   separated by '.'
%
%Example: 
%   Check if s.f1.f11.f111 is a field:
%   isfield(s, 'f1.f11.f111')
%
%   Thorsten.Hansen@psychol.uni-giessen.de  2016-11-10

% extract fields
field = textscan(fields, '%s', 'delimiter', '.');
field = field{1};

yn = true;
structstr = 's'; % name of first input parameter
i = 1;
while yn && i <= numel(field)
  % disp(['yn = isstruct(' structstr ')'])
  eval(['yn = isstruct(' structstr ');']);
  if yn
    % disp(['yn = isfield(' structstr ', ' field{i} ')'])
    eval(['yn = isfield(' structstr ', field{i});']);
  end
  if yn
    structstr = strcat(structstr, '.', field{i});
    i = i + 1;
  end
end