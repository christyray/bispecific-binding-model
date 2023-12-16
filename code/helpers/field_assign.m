function struc = field_assign(struc, pattern, value)

% FIELD_ASSIGN	Assign values to structure fields by field name.
%
%   Assigns the provided value to each field inside the structure whose
%   field name matches the regular expressions pattern given. Will
%   overwrite all fields that match the given pattern.
%
%   USAGE:
%       STRUC = FIELD_ASSIGN(STRUC, PATTERN, VALUE)
%
%   INPUT:
%       STRUC = structure with fields to be changed
%
%       PATTERN = string array with regular expression pattern(s) to search
%       for in the field names
%
%       VALUE = value to assign to matched fields
%
%   OUTPUT:
%       STRUC = Input structure with matching fields reassigned
%
%   NOTES:
%       PATTERN can be any regular expression pattern recognized by MATLAB.
%       MATLAB's regular expressions are similar in functionality to
%       Perl-Compatible Regular Expressions and allow meta-characters,
%       shorthand character classes, and look-arounds.
%
%   See also FIELD_SELECT, STRUCT_INIT, REGEXP, NAMES.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    struc {mustBeA(struc, "struct")}
    pattern {mustBeA(pattern, "string")}
    value
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign Values to Matching Fields

% If multiple patterns were provided, join them with "|"
pattern = join(pattern, "|");

% Get field names
fields = names(struc);

% Find all field names that match the pattern
matches = regexp(fields, pattern);

% Reassign values to each field that matches the pattern
for j = 1:length(fields)
    if matches{j}
        struc.(fields(j)) = value;
    end
end
