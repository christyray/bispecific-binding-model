function tblout = combine(varargin)

% COMBINE	All combinations of the elements in the arrays A1, A2, ..., AN.
%
%   Returns all combinations of the elements in the input objects. Outputs
%   a P-by-N table where P is the product of the number of elements of the
%   N inputs. The column names in the output table correspond to the
%   variable names passed into the function.
%
%	USAGE:
%		TBLOUT = COMBINE(A1, A2, A3, ... AN)
%
%	INPUT:
%		A1, A2, ..., AN = input objects to be combined together; input must
%		be type double, string, character, cell, or table
%
%	OUTPUT:
%		TBLOUT = table of all combinations of the elements of the input
%		objects; column names of the output table correspond to the
%		variable names of the input arguments
%
%	NOTES:
%		If no variable name is provided, COMBINE will assign an arbitrary
%		name ("arg" + a number) to the output column. If there are
%		duplicate names, COMBINE will append a number to the duplicates to
%		make unique names.
%
%       Cell array input must be coercible to a matrix of type double or a
%       string array.
%
%       Tablular input will be merged into a single column for the
%       combination instead of being treating as individual variables;
%       e.g., a table with 3 columns and 5 rows will be treated as 5
%       elements for the combination instead of 15 elements.
%
%       Matrices are reshaped into vectors by column before combination.
%
%	See also NDGRID, RESHAPE, CELL2MAT, CONVERTCHARSTOSTRINGS, TABLE2CELL.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments (Repeating)
    varargin {mustBeA(varargin, ...
        ["double", "string", "char", "cell", "table"])}
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convert Inputs

% Get the number of input arguments
numvec = nargin;

% Get the input arguments
args = varargin;

% Iterate through the input arguments to get the input variable names and
% convert the input arguments

% Initialize empty string array for names
in_names = strings(1,numvec);

% Initialize variable to hold original tables
argstbl = args;

% Determine indices of table inputs
idxtbl = find(cellfun('isclass',args,'table'));

for j = 1:numvec

    % If the input was not named, assign it an arbitrary name
    if isempty(inputname(j))
        in_names(j) = "arg" + string(j);

    % If the input is a repeat name, append a number to the end
    elseif any(matches(in_names, inputname(j)))
        in_names(j) = inputname(j) + string(j);

    % Otherwise, just collect the input name
    else
        in_names(j) = inputname(j);
    end

    % Convert cell arrays to regular arrays and store tables for later use
    type = class(args{j});
    switch type
        case "cell"
            % Convert cell array to matrix, throw error if it cannot be
            % converted
            try
                % Convert character vectors to strings and numbers to
                % matrices
                if iscellstr(args{j})
                    args{j} = convertCharsToStrings(args{j});
                else
                    args{j} = cell2mat(args{j});
                end
            catch
                eid = 'Input:cellArrayCannotConvert';
                msg = 'The input cell array cannot be converted to a matrix.';
                throwAsCaller(MException(eid, msg))
            end

        case "table"
            % For each table, merge all of the columns into one column
            % Prevents combining over each column of the table, instead the
            % table as a whole is combined
            name = in_names(j);
            tbl = args{j};
            argstbl{j} = mergevars(tbl, 1:width(tbl), ...
                "MergeAsTable", 1, "NewVariableName", name);

            % For the actual combine step, insert row numbers in place of
            % the tabular data
            args{j} = 1:size(tbl,1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Combine Inputs

% Enter arguments backwards, so last one (AN) is changing fastest
ii = numvec:-1:1;           % Reverses indices
args = args(1:numvec);      % Uses the reversed indices to reverse the args

% Create a cell array of each input argument repeated to be the same length
[tblout{ii}] = ndgrid(args{ii});

% Reshape all of the repeated input arguments into single columns where
% each column corresponds to an input argument
tblout = cellfun(@(x) reshape(x, [], 1), tblout, 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Table

% Create table of the combined variables, named by the input variables
tblout = table(tblout{:}, 'VariableNames', in_names);

% Replace the numeric stand-ins with the original tabular data
% For each input argument that was a table
for i = 1:length(idxtbl)
    name = in_names(idxtbl(i));
    tbl = argstbl{idxtbl(i)};

    % Use the row numbers from the combine step to repeat the table
    % accordingly
    tbl = tbl(tblout{:,idxtbl(i)},1);

    % Place the repeated table back in the output
    tblout.(name) = tbl.(name);
end

% Split all of the combined tables back into individual columns
tblout = splitvars(tblout, in_names);
