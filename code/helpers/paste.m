function pasted = paste(input, sep)

% PASTE	Concatenate columns of an object into strings with a separator.
%
%   Takes values from each row of an input object and pastes the values
%   into a single string. Useful in combination with combine() to create a
%   string array of all possible combinations of two strings.
%
%   USAGE:
%       PASTED = PASTE(INPUT, SEP)
%
%   INPUT:
%       INPUT = object with columns to paste together
%
%       SEP = separator to use between strings, default is "_"
%
%   OUTPUT:
%       PASTED = string array of pasted table values
%
%   NOTES:
%       TBL can be any data object that is coercible to a string array,
%       e.g., a structure, matrix, vector, or a table. If the coerced array
%       only has a single column, all values of that column will be pasted
%       together.
%
%       The result is that PASTE is equivalent to STRJOIN for string
%       vectors but has the added functionality of joining table columns or
%       structures.
%
%       Any numeric values in TBL will be coerced to strings before
%       pasting.
%
%   See also COMBINE, STRJOIN.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    % Input only needs to be coercible to a string array
    input {mustBeA(input, ...
        ["struct", "table", "double", "char", "string", "cell"])}
    sep string = "_"
end

% Coerce input object to a string array; throw an error if it cannot be
% coerced
type = class(input);
try
    switch type
        % Convert tables to string arrays
        case "table"
            % Convert all table columns to strings
            input = convertvars(input, 1:width(input), "string");

            % Convert the table to a string array
            input = table2array(input);

        % For structures, first convert to cell arrays
        case "struct"
            input = struct2cell(input);

            % Find the most common dimension across the structure fields to
            % use as the joining dimension
            sz = cellfun(@size, input, 'UniformOutput', false);
            sz = horzcat(sz{:});
            n = histcounts(sz);         % Occurrences of each size
            dim = find(n == max(n), 1, 'last');
            % Determine the most common dimension length; if there are
            % multiple of the same length, take the largest

            % For each cell, rotate the values so the rows are the joining
            % dimension
            for i = 1:length(input)
                sz = size(input{i}, [1 2 3]);

                % If the rows match the common dimension
                if sz(1) == dim
                    input{i} = string(input{i});

                % If the columns match the common dimension
                elseif sz(2) == dim
                    % Flip the rows and columns
                    input{i} = string(permute(input{i}, [2 1 3]));

                % If the depth matches the common dimension
                elseif sz(3) == dim
                    % Flip the depth and rows
                    input{i} = string(permute(input{i}, [3 2 1]));
                end
            end

            % Convert the cell array into a string array
            input = horzcat(input{:});
    end
catch
    eid = 'String:couldNotCoerce';
    msg = 'Input object could not be coerced into a string array.';
    throwAsCaller(MException(eid,msg))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Paste Row Values

% If the input is a vector, can simply use strjoin
if any(size(input) == 1)
    row = string(input);
    pasted = strjoin(row, sep);

% Otherwise, for each row in the input, convert all values to strings and
% join the values with the separator
else
    % Initialize empty string vector
    pasted = strings(size(input, 1), 1);

    % Paste each row together
    for i = 1:size(input, 1)
        row = string(input(i,:));
        pasted(i) = strjoin(row, sep);
    end
end
