function [norm, diff, denom] = norm_output(out, id, yin, data, ...
    basis, calc, time, option, denom_in)

% NORM_OUTPUT	Normalize model output based on input data normalization.
%
%	Normalize the model output to convert it to fractional binding for
%	comparison with experimental data. Accepts different normalization
%	types and uses the experimental data to determine which concentrations
%	to use as the basis for the normalization.
%
%	USAGE:
%       [NORM, DIFF, DENOM] = NORM_OUTPUT(OUT, ID, YIN, DATA, BASIS, CALC,
%           TIME, OPTION, DENOM_IN)
%
%   INPUT:
%		OUT = long form table with calculated simulation output; contains
%       columns for input and output ID (corresponding to INPUT and
%       OUTPUT), which calculation was performed, and the time points,
%       species name, and value for each output point; generated by
%       BINDING_SIM
%
%       ID = table with columns for the simulation ID, initial
%       concentration ID, and parameter value ID; used for connecting the
%       initial concentrations used for the simulation with the parameter
%       values used; generated by SETUP_ID
%
%		YIN = table with columns for the simulation ID, time period type,
%		starting and ending times for the time period, species name, and
%		initial concentration for that species for each time period;
%		generated by SETUP_INPUT
%
%       DATA = normalized data points for the experimental data, input as a
%       long form table with columns for the concentration, species, cell
%       type, replicate number, value, and normalization flag; generated by
%       IMPORT_DATA
%
%       BASIS = string with which species were used as the basis for
%       normalizing the data; options are ["Ab", "BS1"]
%
%       CALC = string with which normalization calculation was used for the
%       data; options are ["data", "max"]
%
%       TIME = time point to use for the normalization, defaults to the
%       final time point in the simulation output
%
%       OPTION = string array with additional options for output; options
%       are ["table", "matrix"]
%
%       DENOM_IN = table with the normalization denominators to use for
%       each antibody-cell line combination; generated as an output from
%       this function and used for normalizing additional data and species
%       types that do not directly correspond to the experimental data
%
%	OUTPUT:
%		NORM = long form table with normalized simulation output;
%		contains columns for species name, cell line, and value
%
%       DIFF = long form table with difference between the normalized
%       simulation output and the experimental data, contains columns for
%       species name, cell line, replicate number, and value
%
%       DENOM = table with the denominator used for the normalization of
%       each species and cell line combination
%
%	NOTES:
%		The OUT table needs to share the same time points and
%		initial concentrations for each simulation performed. If there are
%		multiple time points in the table, the final time point will be
%		used for determining the denominator for the normalization. There
%		must a single output value for each time point - most commonly,
%		this will be the concentration for the species matching the
%		experimental data (e.g., the sum of all bound receptors).
%
%       The INPUT table is used to determine which antibodies and cell
%       lines were used for each simulation.
%
%       All of the simulations must be run with the same set of initial
%       concentrations. All of the experimental data points must also share
%       the same initial concentrations, but the simulations do not have to
%       use the same set of initial concentrations as the experimental
%       data. NORM_OUTPUT will compare the simulation output from the
%       closest initial concentration to each data point.
%
%       The "Ab" BASIS normalizes each antibody-cell line
%       combination separately. The "BS1" BASIS normalizes each
%       antibody to the BS1 concentrations in that cell line.
%
%       The "data" CALC normalizes the model output to the same
%       concentrations that were used in the experimental data (determined
%       by the normalization flag assigned when the data is imported by
%       IMPORT_DATA). The "max" CALC normalizes the model output to the
%       maximum binding of the input concentrations present in the
%       experimental data (e.g., if simulations were done at higher
%       concentrations than are present in the experimental data, the
%       normalization will still use the output from the highest
%       concentration present in the experimental data to do the
%       normalization).
%
%	See also BINDING_COST, BINDING_OPT, BINDING_SIM.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    out table
    id table
    yin table
    data table
    basis string {mustBeMember(basis, ["Ab", "BS1"])}
    calc string {mustBeMember(calc, ["data", "max"])}
    time double
    option string {mustBeMember(option, ["table", "matrix"])}
    denom_in table = table()
end

% Remove unnecessary columns from the output table
out = removevars(out, ["OutputID", "Calc", "Species"]);

% Remove duplicate time points from output table (at period junctions)
[~, idxu] = unique(out(:,["ID", "Time"]), "rows", "stable");
out = out(idxu, :);

% Shrink the input table to just the IDs that are present in the simulation
% data and just the relevant columns to streamline the later filtering and
% joining
idxo = unique(out.ID);
idxi = id(ismember(id.ID, idxo), ["ID", "ConcID"]);
input = innerjoin(idxi, yin);
input = input(:, ["ID", "Ab", "Recep", "Period", "Species", "Conc"]);

% Sort data table to ensure it is in a consistent order
data = sortrows(data, ["Species", "Cell", "Conc", "Rep"]);

% If time argument was empty, set the time equal to the last time point
if isempty(time)
    time = max(out.Time);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Organize Model Output and Experimental Data

% -------------------- ANTIBODY-CELL LINE COMBOS --------------------
% Filter input table into a single row for each input ID to make it faster
% to determine the antibodies and cell line for each simulation
[~, idxu] = unique(input.ID);
input_u = input(idxu, :);

% Split all of the receptors in each row into an array and concatenate the
% unique values into a single string to make a column vector of cell types
cat_row = @(x) strjoin(unique(split(x, ["-", "_"])), "");
in_recep = arrayfun(cat_row, input_u.Recep);

% Select the input antibodies as a column vector
in_ab = input_u.Ab;

% Store the unique antibodies separately to use for indexing the inputs
% When multiple antibodies are present, split the string into individual
% antibodies, then keep the unique antibodies for indexing
ab = unique(in_ab);
ab = arrayfun(@(x) unique(split(x, "_")), ab, 'UniformOutput', false);
ab = unique(vertcat(ab{:}));

% Remove M or B for valency from the antibody names to match the
% experimental data
% Using a positive lookahead to match all M or B before an underscore or at
% the end of the string
in_ab = regexprep(in_ab, "(M|B)(?=_|$)", "");

% Isolate the antibody and cell type for each simulation
combos = table(input_u.ID, in_ab, in_recep, ...
    'VariableNames', ["ID", "Species", "Cell"]);

% Add columns to the output table with the antibody-cell combination
% matching the format in the data table, key variable = ID
out = innerjoin(combos, out);

% Determine the unique antibody-cell combinations in the simulations
combos = unique(combos(:, ["Species", "Cell"]), "rows");

% Select only the data that corresponds to one of the input combinations
% The sort doesn't matter because combos is joined onto the sorted data
data = innerjoin(data, combos);

% -------------------- INITIAL CONCENTRATIONS --------------------
% Determine the initial concentration of each antibody for each simulation
initial = unique(input(...
    input.Period == "initial" & ismember(input.Species, ab), ...
    ["ID", "Species", "Conc"]), "rows");

% If there were multiple antibodies present, sum their concentrations
% groupsummary groups the table by the simulation ID, then adds all values
% in the Conc column from the same simulation ID; then rename the columns
% to match the experimental data
initial = groupsummary(initial, "ID", "sum", "Conc");
initial = renamevars(initial(:, ["ID", "sum_Conc"]), "sum_Conc", "Conc");

% Join the total initial antibody concentrations onto the output table
out = innerjoin(initial, out);

% Sort the out table so that it matches the data table
out = sortrows(out, ["Species", "Cell", "Conc", "Time"]);

% Determine the unique concentrations for the model output and the
% experimental data
% Using uniquetol() (default tolerance = 1e-12) so floating-point error
% differences will not cause duplicate concentrations to be returned
% uniquetol() does not have "stable" option, so need to manually sort the
% unique indices to get the data returned in the original order
[~, idxu] = uniquetol(out.Conc);
conc_out = out.Conc(sort(idxu));
[~, idxu] = uniquetol(data.Conc);
conc_exp = data.Conc(sort(idxu));

% Determine which simulations are closest to the initial concentrations
% used in the experimental data
idx = dsearchn(conc_out, conc_exp);
% Output order of dsearchn depends on the order of the second
% argument (it gives the indices in the first array that are
% the closest matches to each value of the second array)

% -------------------- TIME POINT --------------------
% Determine the time point in the simulation output that is closest to the
% requested time for the normalization
time_out = unique(out.Time, "stable");
idx_time = dsearchn(time_out, time);

% -------------------- EXPERIMENTAL DATA --------------------
% Determine if any of the model antibodies are not present in the
% experimental data (e.g., the combination of tocilizumab and 10H2)
in_ab = unique(in_ab);          % Names in in_ab were fixed to match data
ab_exp = unique(data.Species);
missing_ab = in_ab(~ismember(in_ab, ab_exp));

% If any of the model antibodies are not present in the experimental data,
% create "dummy" data for the antibodies to keep dimensions consistent in
% the normalization
if ~isempty(missing_ab)

    % Build the columns of dummy data for the missing antibodies
    cell_exp = unique(data.Cell);
    rep_exp = unique(data.Rep);

    % Create a table of dummy data using the experimental concentrations,
    % cell lines, and replicates, and using 0 for Value and Norm
    dummy = combine(conc_exp, missing_ab, cell_exp, rep_exp, 0, 0);
    dummy = rename(dummy, new=names(data));

    % Add the dummy data onto the original data
    data = [data; dummy];

    % Re-sort the data table to ensure it is still in a consistent order
    data = sortrows(data, ["Species", "Cell", "Conc", "Rep"]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Normalize Values and Calculate Cost

% -------------------- INITIALIZATION --------------------
% Count number of input variables
nconcO = length(conc_out);          % Model output concentrations
nconcE = length(conc_exp);          % Experimental data concentrations
nrep = length(unique(data.Rep));    % Replicates
ntime = length(unique(out.Time));   % Time points
ncombo = height(combos);            % Antibody-cell line combinations

% Initialize output matrices
diff = zeros(nconcE,nrep,ncombo);   % x = conc; y = reps; z = combos
diffID = zeros(nconcE,nrep,ncombo);
norm = zeros(nconcO,ntime,ncombo);  % x = conc; y = time; z = combos
normID = zeros(nconcO,ntime,ncombo);
denom = zeros(ncombo, 1);           % One denominator per combination
colnames = strings(ncombo, 1);

% Initialize counter for filling cost matrix
idxz = 1;

% If a pre-calculated normalization denominator was not provided, create
% an empty table with the same syntax as the output denominator table
if isempty(denom_in)
    denom_in = addvars(combos, zeros(ncombo, 0), ...
        'NewVariableNames', "Value");
end

% -------------------- NORMALIZATION BASIS --------------------
% Perform normalization based on which type of experimental data was given
switch basis
    % Normalize each antibody-cell combination separately
    case "Ab"
        % Determine the unique antibody-cell combinations
        norms = combos;

    % Normalize each antibody to the BS1 values in the same cell line
    case "BS1"
        % Determine the unique cell lines used
        % Allowing norms to sort because then it will match the sorting
        % used by data and out
        norms = unique(combos(:, "Cell"), "rows");
end

% -------------------- CALCULATION --------------------
% Calculate for each normalization basis
for i = 1:height(norms)

    % Filter data (values and normalization indices) and model
    % output for current antibody-cell combination
    % Normalization indices are only used for the "data" calculation
    out_i = innerjoin(out, norms(i,:));
    data_i = innerjoin(data, norms(i,:));
    denom_i = innerjoin(denom_in, norms(i,:));
    Y = out_i.Value;
    ID = out_i.ID;
    exp_i = data_i.Value;
    nrm_i = data_i.Norm;
    denom_i = denom_i.Value;

    % Reshape so the antibodies are in separate pages
    % Time is changing fastest in model output, then concentration
    % Replicates are changing fastest in data, then concentration

    % Reshape fills down columns first, so the fastest changing variables
    % will be in the rows, the second fastest in columns, and the third
    % fastest in the depth
    Y = reshape(Y, ntime, nconcO, []);
    ID = reshape(ID, ntime, nconcO, []);
    exp_i = reshape(exp_i, nrep, nconcE, []);
    nrm_i = reshape(nrm_i, nrep, nconcE, []);

    % Flip the rows and columns to get the desired shape for the outputs
    Y = permute(Y, [2 1 3]);
    ID = permute(ID, [2 1 3]);
    exp_i = permute(exp_i, [2 1 3]);
    nrm_i = permute(nrm_i, [2 1 3]);
    nrm_i = logical(nrm_i(:,1,:));  % All replicates are normalized same
    % rows = conc; cols = time or reps; depth = combos

    % Select the specified time points from the simulations that are
    % closest to the initial concentrations used in the experimental data
    % to calculate the denominator for the normalization
    Yconc = Y(idx,idx_time,:);
    IDconc = ID(idx,idx_time,:);
    % cols = time, so the time index will select the corresponding column

    % If no pre-calculated denominator was provided, calculate the
    % denominator using the provided calculation type
    if isempty(denom_i)

        % Determine the normalization denominator based on the calculation
        switch calc
            case "data"
                % Normalize the model output using the average values for
                % the same concentrations that were normalized in the exp
                % data
                denom_i = mean(Yconc(nrm_i));

            case "max"
                % Normalize the model output using the maximum binding for
                % the antibodies selected as the basis
                % Using maximum binding from same concentrations as data
                denom_i = max(Yconc(:,:,sum(nrm_i) > 1));
        end

    % If there are multiple denominator values, check if they are the same,
    % and throw an error if they are not
    elseif length(denom_i) > 1
        if ~all(denom_i == denom_i(1))
            eid = 'Values:tooManyDifferentInputs';
            msg = ['Too many denominator values provided for the' ...
                'requested normalization.'];
            error(eid, msg)
        else
            % Select just one value if they are all the same
            denom_i = denom_i(1);
        end
    end

    % If the denominator is still empty (e.g., when using the antibody data
    % as a basis and looking at the negative cell lines), replace it with
    % NaN for standardization
    if isempty(denom_i); denom_i = NaN; end

    % Normalize all of the model output and the selected simulations to
    % compare against the experimental data
    Y       = Y ./ denom_i;
    Yconc   = Yconc ./ denom_i;

    % Find difference between normalized model output and data
    nz = size(exp_i, 3);    % # of combos in normalization (# of pages)
    idxz1 = idxz + nz - 1;  % End of indices for current normalization
    diff(:,:,idxz:idxz1) = Yconc - exp_i;
    diffID(:,:,idxz:idxz1) = repmat(IDconc, 1, nrep);
    norm(:,:,idxz:idxz1) = Y;
    normID(:,:,idxz:idxz1) = ID;
    denom(idxz:idxz1) = denom_i;

    % Get names of species and cell lines used
    combo_names = unique([out_i.Species, out_i.Cell], "rows", "stable");
    colnames(idxz:(idxz + nz - 1)) = join(combo_names, ",");

    % Set counter for next iteration
    idxz = idxz + nz;
end

% -------------------- FORMAT OUTPUT --------------------
% Convert outputs into table form if requested
if option == "table"

    % Create column of replicates to match with the differences
    reps = unique(data.Rep);
    reps = repelem(reps, nconcE, 1);

    % Create column of time points to match with the normalized outputs
    times = unique(out.Time, "stable");
    times = repelem(times, nconcO, 1);

    % Reshape the out matrices into two-dimensional matrices where the rows
    % correspond to different replicates/times and initial concentrations
    % and the columns correspond to different antibodies
    % The columns will be stacked so the concentrations are changing
    % fastest, then the replicates/times
    diff = reshape(diff, [], ncombo);
    diffID = reshape(diffID, [], ncombo);
    norm = reshape(norm, [], ncombo);
    normID = reshape(normID, [], ncombo);

    % Make a table from the concentrations, replicates, and differences,
    % and convert the table into long form
    diff = array2table([repmat(conc_exp, nrep, 1), reps, diff], ...
        'VariableNames', ["Conc"; "Rep"; colnames]);
    diff = stack(diff, 3:width(diff), 'IndexVariableName', "Combo", ...
        'NewDataVariableName', "Value");
    diff = addvars(diff, reshape(diffID', [], 1), 'Before', 1, ...
        'NewVariableNames', "ID");

    % Make a table from the concentrations, time points, and normalized
    % values, and convert the table into long form
    norm = array2table([repmat(conc_out, ntime, 1), times, norm], ...
        'VariableNames', ["Conc"; "Time"; colnames]);
    norm = stack(norm, 3:width(norm), 'IndexVariableName', "Combo", ...
        'NewDataVariableName', "Value");
    norm = addvars(norm, reshape(normID', [], 1), 'Before', 1, ...
        'NewVariableNames', "ID");

    % Make a table from the denominators used for the normalization
    denom = array2table(denom, 'VariableNames', "Value");
    denom = addvars(denom, colnames, 'Before', 1, ...
        'NewVariableNames', "Combo");

    % Convert the antibody-cell line combinations into separate columns
    combo_names = split(string(diff.Combo), ",", 2);
    diff = [diff(:,1), diff(:,2), array2table(combo_names, ...
        'VariableNames', ["Species", "Cell"]), diff(:,3), diff(:,5)];
    combo_names = split(string(norm.Combo), ",", 2);
    norm = [norm(:,1), norm(:,2), array2table(combo_names, ...
        'VariableNames', ["Species", "Cell"]), norm(:,3), norm(:,5)];
    combo_names = split(string(denom.Combo), ",", 2);
    denom = [array2table(combo_names, ...
        'VariableNames', ["Species", "Cell"]), denom(:,2)];

    % Sort the rows to match the input tables
    diff = sortrows(diff, ["Species", "Cell", "Conc", "Rep"]);
    norm = sortrows(norm, ["Species", "Cell", "Conc", "Time"]);
    denom = sortrows(denom, ["Species", "Cell"]);
end

end
