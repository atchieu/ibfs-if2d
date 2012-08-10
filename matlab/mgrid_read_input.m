%% ------------------------------------------------------------------------
% function read_input(file)
%
% Reads the input file from the two dimensional MGRID code.
%
%   file        Input file name.
%   names       Parameter names.
%   values      Parameter values.
%--------------------------------------------------------------------------
function s = mgrid_read_input(file)

% Check the correct input arguments
if nargin < 1
  error('Function requires one input argument');
elseif ~ischar(file)
  error('Input must be a string representing a filename');
end

% Number of fields in the input file (should make this automatically read)
% Release 2.10 has 16 fields, my hodge-podge version has 20 fields
n_fields = 16;

% Open file
fid = fopen(file);
if fid == -1
    error('Could not find the file %s', file)
end

% Discard &READ_PARAMETERS
discard = fscanf(fid, '%s', 1); %#ok<NASGU> 

for i = 1:n_fields
    % Get the name of the parameter
    eval_string = sprintf('name%d = fscanf(fid, ''%s'', 1);', i, '%s');
    eval(eval_string);
    
    % Discard =
    discard = fscanf(fid, '%s', 1); %#ok<NASGU> 
    
    % Get the number associated with it
    eval_string = sprintf('number%d = fscanf(fid, ''%s'', 1);', i, '%s');
    eval(eval_string);
    
    % Delete the comma at the end of the number
    eval_string = sprintf('number%d(length(number%d)) = '''';', i, i);
    eval_string2 = sprintf('number%d', i);
    eval(eval_string);
    
    % Check whether it is true, false, or a number
    eval_string1 = sprintf('number%d(1) == ''F''', i);
    eval_string2 = sprintf('number%d(1) == ''T''', i);
    
    if eval(eval_string1) % Check if false
        eval_string = sprintf('number%d = false;', i);
        eval(eval_string);
    elseif eval(eval_string2) % Check if true
        eval_string = sprintf('number%d = true;', i);
        eval(eval_string);
    else % Otherwise it is a number
        eval_string = sprintf('number%d = str2num(number%d);', i, i);
        eval(eval_string);
    end
end

% Start writing it into cell/array so it is more easily accessible
eval_string1 = '{name1';
eval_string2 = '{number1';

for i = 2:n_fields
    eval_string1 = sprintf('%s; name%d', eval_string1, i);
    eval_string2 = sprintf('%s; number%d', eval_string2, i);
end

eval_string1 = sprintf('%s};', eval_string1);
eval_string2 = sprintf('%s};', eval_string2);

% Create the cell (array of strings) and the value array
names  = eval(eval_string1);
values = eval(eval_string2);

% Create a structure
s = cell2struct(values', names', 2);

% Close files out
fclose(fid);

display('Done reading input file...')


return;