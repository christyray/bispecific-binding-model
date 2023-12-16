function struct_out = struct_init(struct_in, value)

% STRUCT_INIT	Initialize a structure from an existing structure.
%
%	Copies all field names to a new structure and sets all values in that
%	structure to the provided input value (default 0). Used for
%   initializing a structure without needing to recreate the field names or
%   for copying field names to another structure.
%
%	USAGE:
%		STRUCT_OUT = STRUCT_INIT(STRUCT_IN, VALUE)
%
%	INPUT:
%		STRUCT_IN = existing structure with correct fields and unwanted
%       values
%
%	OUTPUT:
%		STRUCT_OUT = new structure with unchanged field names and all
%       values set to VALUE
%
%       VALUE = the value(s) to set all fields to; must be type double,
%       string, character, or cell; default is 0
%
%	NOTES:
%		Particularly useful for reinitializing the parameter structure
%		between runs because it does not require re-reading the parameter
%		names with READ_PARAMS.
%
%       Can also be used to create a new structure that has the same fields
%       as an existing structure, e.g., copying all of the species names
%       from the molecule structure to the mass balance structure.
%
%	See also READ_PARAMS, FIELD_SELECT, FIELD_ASSIGN.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    struct_in struct
    value {mustBeA(value, ["double", "string", "char", "cell"])} = 0
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Replace All Values with Initialization Value

% Find current field names
pnames = fieldnames(struct_in);

% Create a cell array of the given value with the same length as the names
pvalues = repmat({value}, length(pnames), 1);

% Match the field names to the values to create an initialized structure
struct_out = cell2struct(pvalues, pnames);

end
