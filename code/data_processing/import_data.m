function data = import_data(file, type, option)

% IMPORT_DATA	Import experimental data from CSV files.
%
%	Import experimental data from CSV files into MATLAB tables to use in
%	the model code.
%
%	USAGE:
%		DATA = IMPORT_DATA(FILE, TYPE)
%
%	INPUT:
%		FILE = name of the file to import data from; if no file extension
%		is provided, the function defaults to ".csv"
%
%       TYPE = string with the type of data in the file, used to determine
%       how to process the data into a table after importing; defaults to
%       "binding"
%
%       OPTION = string array with additional options for output; options
%       are ["noneg"]
%
%	OUTPUT:
%		DATA = a long-form table containing the imported data
%
%	NOTES:
%		Designed to import data into MATLAB that has been copied from the
%       original Prism files into CSV files. The resulting MATLAB tables
%       can be used directly with BINDING_COST and BINDING_DRIVER.
%
%	See also BINDING_DRIVER, BINDING_SIM, BINDING_COST.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    file string
    type string = "binding"
    option string = ""
end

% If the file does not end in an extension type, add the CSV extension
if ~regexpl(file, "\.[a-zA-Z0-9]{2,4}$")
    file = file + ".csv";
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read Data Table from File

% Read in the .csv file as a table
data = readtable(file);

% Process the data table based on which type of data it is
switch type

    % Binding assay data - normalized mean fluorescent intensity
    case "binding"

        % Convert data from wide form to long form
        data = stack(data, 2:width(data), ...
            'IndexVariableName', "Species", ...
            'NewDataVariableName', "Value");

        % Split the species column by underscores and convert into a table
        % The species column has the antibody, cell line, and replicate
        data = convertvars(data, "Species", "string");
        species = split(data.Species, "_");
        varnames = ["Species", "Cell", "Rep"];
        species = array2table(species, 'VariableNames', varnames);
        species = convertvars(species, "Rep", "double");

        % Merge the species columns back with the data
        data = [data(:,1), species, data(:,3)];

        % Determine how the data was normalized
        data = find_norm(data);

        % Convert the concentration to nM if it is not given in nM already
        if any(data.Conc < 0)
            data.Conc = 10.^data.Conc * 1e9;    % Conc given as log10(M)
        end

        % Convert normalized values from percent to decimal
        if any(data.Value > 100)
            data.Value = data.Value / 100;
        end

    % Maximum mean fluorescent intensity data for each antibody
    case "max"

        % Data is already in long form and does not need to be stacked

        % Some data contains replicates and some data does not, need to
        % standardize the species column
        data = convertvars(data, "Species", "string");

        % Count the number of underscores in each string
        n = count(data.Species, "_");

        % If there are not the same number of underscores in each string
        if length(unique(n)) > 1

            % Find strings that need underscores
            idx = n ~= max(n);

            % Append underscore to short strings
            data.Species(idx) = append(data.Species(idx), "_1");
        end

        % Split the species column by underscores and convert into a table
        species = split(data.Species, "_");

        % Name the columns in the species data based on number of fields
        if size(species, 2) == 3
            varnames = ["Species", "Cell", "Rep"];
            species = array2table(species, 'VariableNames', varnames);
            species = convertvars(species, "Rep", "double");
        else
            varnames = ["Species", "Cell"];
            species = array2table(species, 'VariableNames', varnames);
        end

        % Merge the species columns back with the data
        data = [species, data(:,2)];
end

% Remove negative binding if requested
if any(option == "noneg")
    data = remove_neg(data);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions

% -------------------- FIND_NORM --------------------
% Determine how the data was normalized
function data = find_norm(data)

% Sort rows by species and cell so the concatenated indices will match the
% order of the data
data = sortrows(data, ["Species", "Cell", "Conc"]);

% Determine which antibody-cell combinations are in the data
combos = unique(data(:,["Species", "Cell"]), "rows", "stable");

% Initialize cell array to hold indices of normalization
normidx = cell(height(combos), 1);

% Each antibody-cell combination is normalized separately
for i = 1:height(combos)

    % Select current data
    data_i = innerjoin(data, combos(i,:));

    % Initialize normalized flag
    norm = 0;

    % How many values to attempt to average together
    for j = 1:height(data_i)

        % Last start index to use for average
        idxN = height(data_i) - j + 1;

        % For each possible start index
        for k = 1:idxN

            % Find indices to average and calculate the average
            idx = k:k + j - 1;
            avg = mean(data_i.Value(idx));

            % If the values average up to 100, they were used for the
            % normalization and the loop can be stopped
            if ismembertol(avg, 100, 1e-6)
                norm = 1;
                break
            end
        end

        % If normalization was found, exit outer loop
        if norm == 1
            break
        end
    end

    % Indices of normalized values - 0 = not used for normalization; 1 =
    % used for normalization
    normidx{i} = zeros(height(data_i), 1);
    if norm == 1
        normidx{i}(idx) = 1;
    end
end

% Concatenate the normalization indices onto the original data table
normidx = vertcat(normidx{:});
data = addvars(data, normidx, 'NewVariableNames', "Norm");

end

% -------------------- REMOVE_NEG --------------------
% Remove antibody-cell combinations that are not able to bind
function data = remove_neg(data)

% Remove the double-negative cells from the data
data = data(data.Cell ~= "Neg", :);

% Combinations of antibodies and cells that do not result in binding
combos = {["Toci", "8R"], ["H2", "6R"]};

% Remove each unnecessary combination from the data table
for i = 1:length(combos)

    % Find rows that match the combination
    idx_remove = data(data.Species == combos{i}(1) & ...
        data.Cell == combos{i}(2), :);

    % Remove those rows from the data table
    data = setdiff(data, idx_remove, "stable");
end

end
