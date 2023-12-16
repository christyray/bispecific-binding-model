function [values, fields] = field_select(struc, pattern)

% FIELD_SELECT	Create structure with selected fields from input structure.
%
%   Uses a pattern to select fields from an input structure and returns
%   arrays with the values of the selected fields and the selected field
%   names.
%
%   USAGE:
%       [VALUES, FIELDS] = FIELD_SELECT(STRUC, PATTERN)
%
%   INPUT:
%       STRUC = structure to select fields from
%
%       PATTERN = regular expression pattern(s) to match to select fields,
%       given as a string array
%
%   OUTPUT:
%       VALUES = the values originally assigned to the selected fields
%
%       FIELDS = the names of the selected fields
%
%   NOTES:
%       PATTERN can be any regular expression pattern recognized by MATLAB.
%       MATLAB's regular expressions are similar in functionality to
%       Perl-Compatible Regular Expressions and allow meta-characters,
%       shorthand character classes, and look-arounds.
%
%       FIELD_SELECT will return the VALUES as a matrix if it is possible
%       to vertically concatenate the values; otherwise, it will return
%       VALUES as a cell array.
%
%   See also FIELD_ASSIGN, REGEXP.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    struc struct
    pattern string
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select Matching Field Values and Names

% If multiple patterns were provided, join them with "|"
pattern = join(pattern, "|");

% Get field names
all_fields = names(struc)';

% Find and select all field names that match the pattern
matches = regexpl(all_fields, pattern);
fields = all_fields(matches);

% Select values of each field that matches the pattern
values = arrayfun(@(field) struc.(field), fields, 'UniformOutput', false);

% Concatenate the cell array into a matrix if possible
try
    values = vertcat(values{:});
catch ME
    switch ME.identifier
        % If the values cannot be concatenated
        case "MATLAB:catenate:dimensionMismatch"
            % Reshape into a column vector
            values = reshape(values, [], 1);

        % If there was a different error, return the error
        otherwise
            rethrow(ME)
    end
end
