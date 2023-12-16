function data_names = names(object)

% NAMES     Names of fields, data columns, or variables.
%
%   Returns the names of structure fields or table columns from the input
%   object.
%
%   USAGE:
%       DATA_NAMES = NAMES(OBJECT)
%
%   INPUT:
%       OBJECT = a structure, table, or other object
%
%   OUTPUT:
%       DATA_NAMES = a string vector containing the field names or table
%       columns, or an empty string vector if a different object type was
%       provided
%
%   NOTES:
%       For an OBJECT that is not class "struct" or class "table", NAMES
%       will return an empty string vector for compatibility with other
%       functions.
%
%   See also RENAME.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    object
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Names of Input Object

% Different method to get the names depending on the object class
switch class(object)
    case "struct"
        % Transpose the names to make a row vector, like the table names
        data_names = string(fieldnames(object)');

    case "table"
        data_names = string(object.Properties.VariableNames);

    otherwise
        % For all other object types, return an empty string vector
        data_names = strings(0);
end
