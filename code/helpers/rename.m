function object = rename(object, objnames)

% RENAME     Change names of fields or data columns.
%
%   Returns the input object with names of structure fields or table
%   columns changed to provided names.
%
%   USAGE:
%       OBJECT = RENAME(OBJECT, OLD=old_names, NEW=new_names)
%
%   INPUT:
%       OBJECT = a structure or table with named fields or columns
%
%       OLD = (optional); a string vector of the original names of the
%       object in the order that the new names should be applied; default
%       is the internal object name order
%
%       NEW = (optional); a string vector of the new names for the object;
%       default is the value of the OLD argument
%
%   OUTPUT:
%       OBJECT = the original input object with new names for the fields or
%       data columns
%
%   NOTES:
%       If OLD is specified and NEW is not specified, the OBJECT names will
%       remain the same and will be ordered by the order in old. This is
%       convenient shortcut to reorder the columns in a table by name.
%
%       If fewer names are provided in OLD than the length of the input
%       OBJECT, the OBJECT will be subsetted by names in OLD. This is a
%       convenient shortcut to select only specific fields or columns from
%       OBJECT.
%
%       OLD provides the order to use for renaming, and NEW provides the
%       new names to use.
%
%   See also NAMES, MOVEVARS, RENAMEVARS.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    object {mustBeA(object, ["struct", "table"])}
    objnames.old string = names(object)
    objnames.new string
end

% If new names were not provided, use the old names
if ~isfield(objnames, 'new')
    objnames.new = objnames.old;
end

% If the new names and old names are the same as the object names, generate
% a warning
if isequal(names(object), objnames.old, objnames.new)
    wid = 'Input:identical';
    msg1 = 'Both the new and old names match the original names.';
    msg2 = 'No changes will be made to the object.';
    warning(wid, '%s\n%s', msg1, msg2)
end

% Throw error if any of the old names are not present in the input object
if ~all(matches(objnames.old, names(object)))
    missing = objnames.old(~matches(objnames.old, names(object)));
    missing = append("'", missing, "'");
    missing = "Missing: " + strjoin(missing, ", ");
    eid = 'Input:missing';
    msg = 'Input names are not present in input object.';
    msg = sprintf('%s\n%s', msg, missing);
    throwAsCaller(MException(eid, msg))
end

% Test for equal size
if ~isequal(length(objnames.old),length(objnames.new))
    eid = 'Size:notEqual';
    msg = 'Length of the names to replace must be the same as new names.';
    throwAsCaller(MException(eid,msg))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rename Fields or Columns

% Initialize new object to hold transferred data
newobj = struct();

% Iterate through the original object names
for i = 1:length(objnames.old)
    new = objnames.new(i);
    old = objnames.old(i);

    % Move data to new names
    newobj.(new) = object.(old);
end

% If input was a table, need to convert output back to table
switch class(object)
    case "struct"
        % Output updated structure
        object = newobj;

    case "table"
        % Convert structure back to table and output
        object = struct2table(newobj);
end
