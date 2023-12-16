function obj_split = split_object(n, obj)

% SPLIT_OBJECT	Split input objects(s) into smaller parts.
%
%	Splits the input objects(s) into the provided number of smaller tables
%   by row. Useful for splitting the input into a PARFOR loop.
%
%	USAGE:
%		OBJ_SPLIT = SPLIT_OBJECT(N, OBJ1, OBJ2, ...)
%
%	INPUT:
%		N = number of divisions to split the objects(s) into
%
%       OBJ = object(s) to split into smaller objects
%
%	OUTPUT:
%		OBJ_SPLIT = cell array where the rows correspond to the different
%		input OBJs and the columns are the split input objects
%
%	NOTES:
%		If the input OBJ has an ID column, it will be used for indexing the
%		division. Otherwise, the object will be split based on the row
%		number.
%
%	See also PARFOR.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    n double
end

arguments (Repeating)
    obj {mustBeA(obj, ["double", "cell", "table", "string", "char"])}
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Split Input Object(s)

% Initialize output
obj_split = cell(length(obj), n);

% For each input object
for i = 1:length(obj)
    obj_i = obj{i};

    % Determine if there is an ID column
    cols = names(obj_i);
    id = regexpl(cols, "ID", 'ignorecase');

    % If there is an ID column, split based on it
    if any(id)

        % If there are multiple ID columns, split based on the first one
        id = find(id, 1);

        % Determine how many IDs should be in each division
        nrow = max(obj_i{:,id});
        subrow = (nrow)/n;

        % For each division
        for j = 1:n

            % Split into divisions based on ID values
            idx = obj_i{:,id} > round(subrow * (j - 1)) & ...
                obj_i{:,id} <= round(subrow * j);
            obj_split{i,j} = obj_i(idx, :);
        end

    % If there is not an ID column, split based on the number of rows
    else
        % Determine how many rows should be in each division
        nrow = size(obj_i, 1);
        subrow = (nrow)/n;

        % For each division
        for j = 1:n

            % Split into divisions based on row number
            idx = round(subrow * (j - 1) + 1):round(subrow * j);
            obj_split{i,j} = obj_i(idx, :);
        end
    end
end
