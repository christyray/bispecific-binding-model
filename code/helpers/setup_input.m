function [tbl, chk_tbl] = setup_input(id, varargin)

% SETUP_INPUT	Generate table of input species and values for simulation.
%
%	Takes input of molecules or parameters used in the system and the
%	values that will be used for those inputs. Accepts either the actual
%	molecule or parameter name or a regular expression to select multiple
%	values. Generates a table of all combinations of the species and
%	values.
%
%	USAGE:
%		[TBL, CHK_TBL] = SETUP_INPUT(ID, AB, RECEP, PERIODS, TIMES, ...
%           SPECIES1, VALUES1, ...)
%       [TBL, CHK_TBL] = SETUP_INPUT(ID, SPECIES1, VALUES1, ...)
%
%	INPUT:
%       ID = the number to start the ID column with
%
%		AB = a string vector with the antibodies to use in the complexes;
%		required to generate initial concentations
%
%       RECEP = a cell array of string vectors with the receptor(s) to use
%       in the complexes; the rows of the array should correspond to the
%       antibodies in AB; required to generate initial concentations
%
%       PERIODS = string or cell array of time period types for the
%       simulations; types should correspond to the TIMES argument values;
%       required to generate initial concentations
%
%       TIMES = matrix or cell array of time values for the simulations;
%       should be given in a vector where each value is the start of a new
%       time period and the final value is the end of the simulation;
%       required to generate initial concentations
%
%		SPECIES = string array for the species to use in the table; species
%       that are paired and should not be iterated togther should be input
%       as a single argument; species that should be iterated together
%       should be input in separate arguments
%
%       VALUES = matrix or cell array of values to use in the table; VALUES
%       should be paired with their corresponding SPECIES argument; cells
%       or matrix columns/depth should correspond to the order of species
%       in the SPECIES argument
%
%	OUTPUT:
%		TBL = table with columns for the simulation ID, time period type,
%		starting and ending times for the time period, species name, and
%		input value for that species
%
%       CHK_TBL = cell array with tables for each set of input SPECIES and
%       VALUES arguments; each table contains a row for each set of
%       concentrations/values for all of the input species in that
%       argument; shows how the inputs are paired together and how the
%       values are assigned to time periods; used for debugging TBL
%
%	NOTES:
%		This function generates the yin and p tables to be combined
%		together and used as inputs to the BINDING_SIM function.
%
%		SETUP_INPUT generates all possible combinations of SPECIES and
%       VALUES that are provided in separate arguments. Inputs that should
%       have all possible combinations generated should be given together
%       in a single argument.
%
%       The AB, RECEP, PERIODS, and TIMES arguments are required to
%       generate the initial concentrations table, and they are not
%       required to generate the parameters table.
%
%	See also COMBINE, SETUP_ID, BINDING_SIM, CREATE_PARAMS.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    id (1,1) double
end

arguments (Repeating)
    varargin
end

% For both the `recep` case and the `values` case, the number of input
% arguments must be odd because SPECIES and VALUES are paired
if mod(nargin, 2) ~= 1
    eid = 'Arguments:incorrectNumber';
    msg = 'Incorrect number of input arguments provided.';
    throwAsCaller(MException(eid, msg))
end

% Determine which type of arguments were passed in
% The second argument is either `recep`, which can be a string or a cell
% array of strings, or `values`, which can be double or a cell array of
% double
if class(varargin{2}) == "cell"
    arg2 = class(varargin{2}{1});
else
    arg2 = class(varargin{2});
end

% If the second argument is recep, then ab, periods, and times were also
% provided
if arg2 == "string"
    [ab, recep, periods, times] = deal(varargin{1:4});
    species = varargin(5:2:nargin-1);
    values = varargin(6:2:nargin-1);
% If the second argument is values, then only species and values pairs were
% given
else
    [ab, recep] = deal("");
    periods = "params";
    times = [0; 0];
    species = varargin(1:2:nargin-1);
    values = varargin(2:2:nargin-1);
end

% Validate that the antibodies and receptors are the same length and
% convert the receptors to a cell array for correct indexing
[ab, recep] = validate_ab(ab, recep);

% Ab and recep are concatenated together so the column dimensions will
% be consistent when they are stacked on top of each other in a table
% Different antibodies are separated by underscores
% Receptors for different antibodies are separated by an underscore,
% receptors for the same antibody are separated by a dash
ab = join(ab, "_");
recep = cellfun(@(x) join(x, "-"), recep, 'UniformOutput', false);
recep = join(horzcat(recep{:}), "_");

% ==================== PERIODS ====================
% Standarized format: each set of time periods in a separate cell,
% individual period names and time points in column vectors inside the
% cells, no duplicate period names

% -------------------- ERRORS --------------------
% Determine the number of input species and dimensions of input values
nper = size(periods);
ntime = size(times);

% Error if period and time set dimensions do not match
if strcmp(class(periods), class(times)) && ...
        ~(all(nper == ntime) || all(nper == flip(ntime)))
    eid = 'Size:notEqual';
    msg = 'Provided time and period sets are not same size.';
    throwAsCaller(MException(eid, msg))
elseif ~(any(nper == ntime) || any(nper == flip(ntime)))
    eid = 'Size:notEqual';
    msg = 'Provided time and period sets are not same size.';
    throwAsCaller(MException(eid, msg))
end

% -------------------- CONVERSION --------------------
% Convert periods and times to cell arrays for correct indexing
if class(periods) == "string" && class(times) == "double"

    % If there is one set of periods, just wrap inputs in cell array
    if any(nper == 1) && any(ntime == 1)
        periods = {periods};
        times = {times};
    else
        % Determine which dimensions of each array correspond to different
        % period sets and reshape so rows = sets and cols = periods
        % Both arrays will have the same number of sets, but the time array
        % will have one more entry for the periods (since it has start and
        % end times)
        if nper(1) == ntime(1) && (nper(2) + 1) == ntime(2)
            % No transposing needed, everything in correct orientation
            % Included in order to bypass the other ifs
        elseif nper(1) == ntime(2) && (nper(2) + 1) == ntime(1)
            times = times';
        elseif nper(2) == ntime(1) && (nper(1) + 1) == ntime(2)
            periods = periods';
        elseif nper(2) == ntime(2) && (nper(1) + 1) == ntime(1)
            periods = periods';
            times = times';
        end

        % Convert times and periods to cell arrays with each period in a
        % different cell
        rowdim = ones(1, size(periods,1));
        periods = mat2cell(periods, rowdim);
        times = mat2cell(times, rowdim);
    end

elseif class(times) == "double"

    % Reshape so the rows correspond to different sets of periods
    if ntime(1) ~= nper(1)
        times = times';
    end

    % Convert to a cell array with each period set in a different cell
    rowdim = ones(1, size(times, 1));
    times = mat2cell(times, rowdim);

elseif class(periods) == "string"

    % Reshape so the rows correspond to different sets of periods
    if nper(1) ~= ntime(1)
        periods = periods';
    end

    % Convert to a cell array with each period set in a different cell
    rowdim = ones(1, size(periods, 1));
    periods = mat2cell(periods, rowdim);
end

% Check individual time period names and values
for i = 1:length(periods)

    % Reshape the times and periods into column vectors
    times{i} = reshape(times{i}, [], 1);
    periods{i} = reshape(periods{i}, [], 1);

    % Make all of the period names unique values so they can be used as
    % table column names
    periods{i} = make_unique(periods{i});

    % Error if individual period and time dimensions do not match
    if length(times{i}) ~= (length(periods{i}) + 1)
        eid = 'Size:notEqual';
        msg = 'Provided times and periods are not same length.';
        throwAsCaller(MException(eid, msg))
    end
end

% ==================== VALUES ====================
% Standardized format: each separate species argument in a separate cell,
% each species within a single argument in a cell within the argument cell,
% values within the cells shaped so the rows = different simulations,
% columns = time periods
% Top level = separate arguments, second level = separate species in same
% argument, third level = values

% Fix the species input if receptors in a cell array were given
for i = 1:length(species)
    if class(species{i}) == "cell"
        % Make all cells row vectors
        species{i} = cellfun(@(x) reshape(x, 1, []), ...
            species{i}, 'UniformOutput', false);
        % Concatenate into a non-cell array
        species{i} = horzcat(species{i}{:});
    end

    % Keep only unique values
    species{i} = unique(species{i}, "stable");

    % Add "IL" to receptor names if it is missing
    species{i} = regexprep(species{i}, "^([68]R)$", "IL$1");
end

% Iterate through each pair of inputs to confirm sizes match
for i = 1:length(values)

    % -------------------- ERRORS --------------------
    % Determine the number of input species and dimensions of input values
    nspec = length(species{i});
    nval = size(values{i}, [1 2 3]);

    % Error if number of input species does not match value dimensions
    if nspec == 1 && class(values{i}) == "cell" && any(nval > 1)
        eid = 'Size:notEqual';
        msg = 'Too many input values provided; only one species given.';
        throwAsCaller(MException(eid, msg))
    elseif nspec == 1 && ~any(nval == 1)
        eid = 'Size:notEqual';
        msg = 'Too many input values provided; only one species given.';
        throwAsCaller(MException(eid, msg))
    elseif nspec > 1 && ~any(nval == nspec)
        eid = 'Size:notEqual';
        msg = 'Provided species and values are not same size.';
        throwAsCaller(MException(eid, msg))
    end

    % -------------------- CONVERSION --------------------
    % Remove any extra dimensions from the values and recalculate the size
    values{i} = squeeze(values{i});
    nval = size(values{i}, [1 2 3]);

    % All given numbers of time periods
    nper = cellfun(@length, periods);

    % Convert input matrices to cell arrays in correct shape for indexing
    % If there is one species, just wrap inputs in cell array
    if class(values{i}) == "double" && nspec == 1

        % If values are a vector, reshape to a column vector
        if isvector(values{i})
            values{i} = reshape(values{i}, [], 1);
        % If values are a matrix, make columns correspond to time periods
        elseif ~any(nval(2) == nper)
            values{i} = values{i}';
        end

        % Convert to cell array
        values{i} = values(i);

    % If the input is two-dimensional
    elseif class(values{i}) == "double" && nval(3) == 1

        % Reshape so the rows correspond to species
        if nval(1) ~= nspec
            values{i} = values{i}';
        end

        % Convert to cell array with each species in different cell
        nrow = size(values{i}, 1);
        celldim = ones(nrow, 1);
        values{i} = mat2cell(values{i}, celldim);

        % Reshape the values in the cell array so the columns correspond to
        % time periods
        for j = 1:length(values{i})
            if ~any(size(values{i}{j}, 2) == nper)
                values{i}{j} = values{i}{j}';
            end
        end

        % Reshape the resulting cell array to a row vector
        values{i} = reshape(values{i}, 1, []);

    % If the input is three-dimensional
    elseif class(values{i}) == "double"

        % Reshape so the depth corresponds to species
        if nval(3) ~= nspec && nval(2) == nspec
            values{i} = permute(values{i}, [1 3 2]);
        elseif nval(3) ~= nspec && nval(1) == nspec
            values{i} = permute(values{i}, [3 2 1]);
        end

        % Reshape so the columns correspond to time periods
        if ~any(size(values{i}, 2) == nper)
            values{i} = permute(values{i}, [2 1 3]);
        end

        % Convert to cell array with each species in different cell
        [nrow, ncol, ndep] = size(values{i}, [1 2 3]);
        celldim = ones(ndep, 1);
        values{i} = mat2cell(values{i}, nrow, ncol, celldim);

        % Reshape the resulting cell array to a row vector
        values{i} = reshape(values{i}, 1, []);

    % If the input is a cell array
    elseif class(values{i}) == "cell"

        % For each species in the cell array
        for j = 1:length(values{i})
            % Remove any extra dimensions from the values
            values{i}{j} = squeeze(values{i}{j});
            nval = size(values{i}{j});

            % If the values are a vector, reshape to a column vector
            if any(nval == 1)
                values{i}{j} = reshape(values{i}{j}, [], 1);

            % If the values are a matrix, reshape so the columns correspond
            % to time periods
            elseif ~any(nval(2) == nper)
                values{i}{j} = values{i}{j}';
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Tables of Converted Values for Verification

% Find which set of time periods has the most periods to use as the names
% for columns of the displayed table
[~, max_idx] = max(cellfun(@length, periods));
per_names = periods{max_idx};

% Initialize output cell array
chk_tbl = cell(length(species), 1);

% For each separate species argument, create a table where each row is one
% set of concentrations for all of the input species in that argument -
% shows how the inputs are paired together and assigned to time periods
for i = 1:length(species)

    % Initialize cell array for tables for each species in one argument
    chk_tbl{i} = cell(length(species{i}), 1);

    % For each species in a single argument, make a table where the first
    % column is the species name and the other columns contain the values
    % for each time period
    for j = 1:length(species{i})
        spe = string(species{i}{j});
        val = values{i}{j};
        per = per_names(1:size(val, 2)) + "_" + j;

        % Convert the concentration values to a table with a column for
        % each time period and add a column for the species name
        chk_tbl{i}{j} = array2table(val, 'VariableNames', per);
        chk_tbl{i}{j} = addvars(chk_tbl{i}{j}, ...
            repmat(spe, size(val, 1), 1), ...
            'NewVariableNames', "Species_" + j, 'Before', 1);
    end

    % Concatenate all of the species in a single argument into one table
    chk_tbl{i} = horzcat(chk_tbl{i}{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Combine Input Values into Table

% Initialize output variable
tbl = cell(length(periods), 1);

% Initialize ID counter
id_start = id;

% For each set of time periods given
for i = 1:length(periods)

    % Current time vector and period type
    times_i = times{i};
    periods_i = periods{i};

    % Reshape the times to two columns with starting and ending time
    times_i = [times_i(1), repmat(times_i(2:end-1)', 1, 2), times_i(end)];
    times_i = reshape(times_i, [], 2);

    % Create a table with the time points and period types
    times_i = array2table(times_i);
    times_i = addvars(times_i, periods_i);
    times_i = rename(times_i, new=["Start", "End", "Period"]);

    % Initialize empty cell array for tables for each species
    tbl_i = cell(size(species));

    % Create tables from the input species and names
    for j = 1:length(tbl_i)

        % Make a table for all of the species for this time period
        tbl_ij = maketbl(species{j}, values{j}, periods_i);

        % Concatenate outputs
        tbl_i{j} = tbl_ij;
    end

    % Generate a table of all combinations of the values in the inputs
    tbl_i = combine(tbl_i{:});

    % Un-merge the species columns into separate columns
    tbl_i = splitvars(tbl_i, 1:width(tbl_i));

    % Create an ID column
    id = (id_start:height(tbl_i) + id_start - 1)';
    id_start = height(tbl_i) + id_start;
    tbl_i = addvars(tbl_i, id, 'Before', 1);

    % Create columns for the antibodies and receptors
    ab_col = repmat(ab, height(tbl_i), 1);
    recep_col = repmat(recep, height(tbl_i), 1);
    tbl_i = addvars(tbl_i, ab_col, recep_col, 'After', 1, ...
        'NewVariableNames', ["ab", "recep"]);

    % Convert the table from wide form to long form
    tbl_i = stack(tbl_i, 4:width(tbl_i), ...
        'IndexVariableName', "Label", ...
        'NewDataVariableName', "Value");

    % Split the labels column into columns for periods and species
    tbl_i = convertvars(tbl_i, "Label", "string");

    % Replace the final underscore with a dash to prevent extra matching
    % when the label is a regular expression
    label = regexprep(tbl_i.Label, "_([^_]*)$", "-$1");
    label = split(label, "-");

    % Add an additional column for the species if species name was not
    % attached to the period type (i.e., if there is only one species)
    if size(label, 2) == 1
        label(:,2) = label;     % Period type should be second column
        label(:,1) = repmat(species{1}, size(label, 1), 1);
    end

    % Rebuild the table
    tbl_i = table(tbl_i.id, tbl_i.ab, tbl_i.recep, ...
        label(:,2), label(:,1), tbl_i.Value, 'VariableNames', ...
        ["ID", "Ab", "Recep", "Period", "Species", "Value"]);

    % Join the time period labels with the start and end times
    tbl_i = join(tbl_i, times_i);

    % Reorganize the table
    tbl_i = movevars(tbl_i, ["Start", "End"], 'After', "Period");
    tbl_i = sortrows(tbl_i, ["ID", "Start"]);

    % Clean up labels
    tbl_i.Period = regexprep(tbl_i.Period, "[0-9]+$", "");
    if periods{i} == "params"
        tbl_i = removevars(tbl_i, ...
            ["Ab", "Recep", "Period", "Start", "End"]);
        tbl_i = renamevars(tbl_i, ["ID", "Species"], ["ParamID", "Param"]);
    else
        tbl_i = renamevars(tbl_i, ["ID", "Value"], ["ConcID", "Conc"]);
    end

    % Concatenate output
    tbl{i} = tbl_i;
end

% Concatenate output
tbl = vertcat(tbl{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Function

% Helper function for combining the species and values into a table
function tbl = maketbl(species, values, type)

% Initialize cell array for output
tbl = cell(length(species),1);

% For each species
for i = 1:length(species)

    % Remove extraneous periods
    nper = length(type);
    nval = size(values{i}, 2);

    % Keep the smaller number of dimensions between period types and
    % concentrations
    if nper > nval
        type = type(1:nval);
    else
        values{i} = values{i}(:,1:nper);
    end

    % Make a table of the input values
    tbl_i = array2table(values{i}, 'VariableNames', type);

    % Combine all of the columns into one table
    tbl_i = mergevars(tbl_i, 1:length(type), ...
        'MergeAsTable', true, 'NewVariableName', species(i));

    % Store table for output
    tbl{i} = tbl_i;
end

% Combine all of the species tables into a single table
tbl = horzcat(tbl{:});

end
