function binding_driver(importfile, options)

% BINDING_DRIVER	Antibody-receptor binding model driver file.
%
%	Driver file for the determination of the antibody-receptor binding
%   rates from the flow cytometry equilibrium binding data and simulations
%   with the optimal rates
%
%	USAGE:
%		BINDING_DRIVER(IMPORTFILE, OPTIONS)
%
%	INPUT:
%		IMPORTFILE = import driver variables saved from the output of
%		SETUP_VARIABLES
%
%       OPTIONS = paired name-value arguments with options for running the
%       driver file, including how many workers to use, the memory to
%       reserve, wall time, etc.
%
%	NOTES:
%		Organized as a function to be able to easily run multiple
%		simulations with different options.
%
%	See also BINDING_SIM, BINDING_OPT, BINDING_COST, SETUP_VARIABLES.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    importfile (1,:) string = strings(0)
    options.Sim string {mustBeMember(options.Sim, ["optimization", ...
        "binding-curve", "time", "concentration", "monovalent", ...
        "compare-ab", "compare-recep", "local", "global", ...
        "compare-opt"])} = "time"
    options.Norm string {mustBeMember(options.Norm, ...
        ["bs1-data", "ab-data", "bs1-max", "ab-max", "none"])} = "none"
    options.Save (1,:) logical = true
    options.Log (1,:) logical = false
    options.VarOnly (1,:) logical = false
    options.Messages (1,:) double = 25
    options.NSplit (1,1) double = 100
    options.Benchmark (1,:) logical = false
    options.BenchmarkFraction (1,1) double = 0.01
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User-Specified Variables

% Import the driver variables from previous simulation or generated file
if isempty(importfile)
    % Generate the driver variables from preset options
    switch options.Sim
        case "optimization"
            varfile = setup_variables(...
                "Name", options.Sim, ...
                "Case", "opt", ...
                "Antibodies", "binding", ...
                "Periods", "experimental", ...
                "Norm", options.Norm ...
                );
        case "binding-curve"
            varfile = setup_variables(...
                "Name", options.Sim, ...
                "Case", "sim", ...
                "Antibodies", "all", ...
                "Periods", "experimental", ...
                "Conc", "high", ...
                "Params", "optimal", ...
                "Outputs", "experimental", ...
                "Norm", "bs1-data" ...
                );
        case "time"
            varfile = setup_variables(...
                "Name", options.Sim, ...
                "Case", "sim", ...
                "Antibodies", "binding", ...
                "Periods", "long-washout", ...
                "Conc", "lowres", ...
                "Params", "optimal", ...
                "Outputs", "base" ...
                );
        case "concentration"
            varfile = setup_variables(...
                "Name", options.Sim, ...
                "Case", "sim", ...
                "Antibodies", "compare-double", ...
                "Periods", "long-initial", ...
                "Conc", "highres", ...
                "Receptors", "total", ...
                "Params", "optimal", ...
                "Outputs", "endpoint" ...
                );
        case "monovalent"
            varfile = setup_variables(...
                "Name", options.Sim, ...
                "Case", "sim", ...
                "Antibodies", "compare-double", ...
                "Periods", "long-initial", ...
                "Conc", "highres", ...
                "Receptors", "total", ...
                "Params", "monovalent", ...
                "Outputs", "endpoint" ...
                );
        case "compare-ab"
            varfile = setup_variables(...
                "Name", options.Sim, ...
                "Case", "sim", ...
                "Antibodies", "compare-double", ...
                "Periods", "long-initial", ...
                "Conc", "highres", ...
                "Receptors", "both-lowres", ...
                "Params", "optimal", ...
                "Outputs", "endpoint" ...
                );
        case "compare-recep"
            varfile = setup_variables(...
                "Name", options.Sim, ...
                "Case", "sim", ...
                "Antibodies", "compare-double", ...
                "Periods", "long-initial", ...
                "Conc", "lowres", ...
                "Receptors", "both", ...
                "Params", "optimal", ...
                "Outputs", "endpoint" ...
                );
        case "local"
            varfile = setup_variables(...
                "Name", options.Sim, ...
                "Case", "local", ...
                "Antibodies", "bs1", ...
                "Periods", "initial", ...
                "Conc", "local", ...
                "Receptors", "local", ...
                "Params", "optimal", ...
                "Outputs", "local" ...
                );
        case "global"
            varfile = setup_variables(...
                "Name", options.Sim, ...
                "Case", "global", ...
                "Antibodies", "bs1", ...
                "Periods", "long-initial", ...
                "Conc", "lowres", ...
                "Receptors", "lowres", ...
                "Params", "optimal", ...
                "Outputs", "key" ...
                );
        case "compare-opt"
            varfile = setup_variables(...
                "Name", options.Sim, ...
                "Case", "sim", ...
                "Antibodies", "negs", ...
                "Periods", "experimental", ...
                "Conc", "standard", ...
                "Receptors", "none", ...
                "Params", "import", ...
                "Outputs", "experimental", ...
                "Norm", options.Norm ...
                );
    end

    % Load in the variables saved by the setup function
    var = load(varfile, '*');
else
    % Load in the variables from the previous simulation
    varfile = importfile + "_vars.mat";
    var = load(varfile, '*');

    % Save the variables with the filename syntax for current simulation
    d = string(datetime("now", "Format", "yyyy-MM-dd"));
    varfile = var.folder + "/" + ...
        join([d, var.name, "vars"], "_") + ".mat";
    save(varfile, '-struct', 'var')
end

% Set the output file name for remaining outputs
d = string(datetime("now", "Format", "yyyy-MM-dd"));
outfile = join([d, var.name], "_");

% If the function was run exclusively to save input variables, stop here
if options.VarOnly
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup and Flags

% Start code timer to measure full run time of driver file, used as an
% argument to the code_time() helper function
tic_driver = tic;

% If saving standard output to a log file, get the file ID
if options.Log
    % Create the directory for log files if it doesn't exist
    logpath = "output/logs/";
    if ~exist(logpath, "dir")
        mkdir(logpath);
    end

    % Create the log file and get its file ID
    logfile = "output/logs/" + outfile + ".log";
    logID = fopen(logfile, 'w');
else
    % FileID of 1 means standard out = will write to screen
    logID = 1;
end

% Close previous parallel pool if needed, then start parallel pool
delete(gcp('nocreate'))    % Closes previous parallel pool if needed
poolobj = parpool(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulations

% Select specific simulation to run
switch var.sim

    % Run simulations for varying system conditions; includes options to
    % perform univariate global sensitivity and to calculate normalized
    % model output
    case {"sim", "global", "local"}

        % --------------- INPUTS ---------------
        % Create molecule and parameter structures
        [m, p, cond] = setup_case(var.setup);

        % Set antibody concentrations
        concAb = logspace(var.minconc, var.maxconc, ceil(var.nconc))';

        % If also using the exact experimental concentrations
        if var.expconc == 1
            conc_exp = unique(vertcat(var.data{:}).Conc);
            concAb = unique([concAb; conc_exp]);
        end

        % Add a column for the second time period
        if ismember("washout", var.period_types)
            concAb = [concAb, zeros(length(concAb), 1)];
        % Adding 0 concentration of all antibodies because the second
        % period is just present here to use a different time step value
        elseif ismember("addition", var.period_types) && ...
                var.options.Periods == "long-initial"
            concAb = [concAb, zeros(length(concAb), 1)];
        end

        % Generate initial concentrations for each receptor depending on
        % what is being varied - the ratio of IL6R to IL8R or the total
        % receptor concentration, or use the default experimental receptor
        % concentration
        concR = receptor_conc(var.options.Receptors, ...
            var.minrecep, var.maxrecep, var.nrecep, var.ratio);

        % Combine antibodies and receptors into initial concentration table
        yin = antibody_conc(var.ab, var.recep, concAb, concR, ...
            var.period_types, var.period_times, var.setup);

        % --------------- PARAMETERS ---------------
        % Set the binding parameter values based on the case
        switch var.sim
            % For regular simulations, use the parameter values from the
            % simulation inputs
            case "sim"
                % Make parameter table
                params = setup_input(1, var.knames, var.kvalues);

            % For the global univariate sensitivity analysis, vary the
            % parameter values based on the ranges specified in the inputs
            case "global"
                % Set up anonymous functions to generate ranges of
                % parameter values
                kexp = @(x) floor(log10(x));
                logvalues = @(x,m,n) ...
                    logspace(kexp(x) - m, kexp(x) + m, n) * x/(10^kexp(x));

                % Set parameter values for simulations - each parameter is
                % varied individually and used with the baseline values for
                % the other parameters
                all_values = cell(length(var.knames), 1);
                for i = 1:length(var.knames)
                    vary_i = logvalues(...
                        var.kvalues(i), var.nmag, var.nkvalues)';
                    all_values{i} = repmat(var.kvalues, var.nkvalues, 1);
                    all_values{i}(:,i) = vary_i;
                end
                all_values = vertcat(all_values{:});

                % Add the baseline values as the first simulation
                all_values = [var.kvalues; all_values];

                % Remove duplicate parameter sets
                all_values = unique(all_values, 'rows', 'stable');

                % Make parameter table
                params = setup_input(1, var.knames, all_values);

            % For the local univariate sensitivity analysis, vary the
            % parameter values by a set percentage
            case "local"

                % Set parameter values for simulations - each parameter is
                % varied individually and used with the baseline values for
                % the other parameters
                all_values = repmat(var.kvalues, length(var.kvalues), 1);
                all_values = all_values + ...
                    all_values .* eye(length(var.kvalues)) * 0.1;

                % Add the baseline values as the first simulation
                all_values = [var.kvalues; all_values];

                % Make parameter table
                params = setup_input(1, var.knames, all_values);
        end

        % Create combined ID table
        id = setup_id(yin, params);

        % --------------- OUTPUTS ---------------
        % Set up output calculation table
        output = setup_output(var.interest, var.value, var.add, var.time);

        % --------------- SIMULATIONS ---------------
        % Run simulations for all of the input variables

        % If benchmarking the current simulation setup, select only a
        % random fraction of simulation IDs to run
        if options.Benchmark
            [id, yin, params] = sample_inputs(...
                id, yin, params, options.BenchmarkFraction);
        end

        % Number of partitions to send to parallel workers
        n_split = options.NSplit;

        % Run model simulations in parallel
        out = run_simulations(n_split, options.Messages, logID, ...
            id, yin, params, output, m, p, var.odeopts);

        % --------------- CALCULATIONS AND SAVING ---------------
        % Calculate normalized model output if requested, and set the
        % variables to be saved

        % If no normalization basis was given, just save the output
        if var.options.Norm == "none"

            % Split the merged columns in the output table
            output = splitvars(output);

            % Which variables to save as .csv files
            savevars = struct("id", id, "out", out, "yin", yin, ...
                "params", params, "output", output);

        % If a normalization basis was given, normalize the output
        elseif var.sim == "sim"
            % Set the time for the normalization from the experimental
            % conditions
            norm_time = (cond.primary + cond.secondary) * 3600;

            % Calculate normalized model output using provided variables
            [norm, ~, denom] = calculate_norm(norm_time, ...
                var.data, var.data_fn, var.basis, var.calc, ...
                id, yin, params, output, out, logID);

            % Split the merged columns in the output table
            output = splitvars(output);

            % Sort the normalized values in the same way as the output
            norm = sortrows(norm, ["Basis", "Calc", "ID", "Time"]);

            % Which variables to save as .csv files
            savevars = struct("id", id, "out", out, "yin", yin, ...
                "params", params, "output", output, "norm", norm, ...
                "denom", denom);

        % For the univariate sensitivity analysis, calculate the cost when
        % using the varied parameter values
        elseif var.sim == "global"
            % Set the time for the normalization from the experimental
            % conditions
            norm_time = (cond.primary + cond.secondary) * 3600;

            % Calculate cost of the normalized model output
            [~, cost] = calculate_norm(norm_time, var.data, var.data_fn,...
                var.basis, var.calc, id, yin, params, output, out, logID);

            % --------------- CONVERSION AND SAVING ---------------
            % Split the merged columns in the output table
            output = splitvars(output);

            % Which variables to save as .csv files
            savevars = struct("id", id, "out", out, "cost", cost, ...
                "yin", yin, "params", params, "output", output);
        end

    % Optimize binding rate constants to flow cytometry equilibrium data
    case "opt"

        % -------------------- OPTIMIZATION --------------------
        % Start timer for optimization code
        tic_opt = tic;

        % Initialize structures and experimental conditions
        [m, p] = setup_case(var.setup);

        % Initialize outputs
        nopt = size(var.initial, 1);
        optimal = zeros(nopt, var.nparams);
        knames = strings(nopt, var.nparams);
        final_cost = zeros(nopt, 1);
        flag = zeros(nopt, 1);
        iter = zeros(nopt, 1);

        % Create the directory for log files if it doesn't exist
        diarypath = "output/logs/" + outfile;
        if ~exist(diarypath, "dir")
            mkdir(diarypath);
        end

        % Initialize the progress message function, create a data queue
        % object to receive messages from parallel workers, and display
        % progress messages throughout code execution
        display_progress(nopt, options.Messages, tic_opt, logID)
        queue = parallel.pool.DataQueue;
        afterEach(queue, @display_progress)

        % Specify variables to pass to parfor loop to avoid passing entire
        % variable structure
        [initial, kopt_names, lb, ub, kin, kin_names, data, ...
            basis, calc, ab, recep, setup] = deal(var.initial, ...
            var.kopt_names, var.lb, var.ub, var.kin, var.kin_names, ...
            var.data, var.basis, var.calc, var.ab, var.recep, var.setup);

        % Run optimization
        parfor i = 1:size(initial,1)
            % Set the diary file for logging the optimization progress
            istr = pad(string(i), 3, "left", "0");
            diaryfile = diarypath + "/" + outfile + "_" + istr;

            [optimal(i,:), knames(i,:), final_cost(i), flag(i), ...
                iter(i)] = binding_opt(initial{i,:}, kopt_names, ...
                lb, ub, kin, kin_names, data, basis, calc, ab, recep, ...
                setup, m, p, "Display", "iter-detailed", ...
                "Diary", diaryfile);

            % Ping the data queue after each optimization is complete to
            % update the progress table
            send(queue, [])
        end

        % Display time elapsed for optimization
        elapsed = code_time(tic_opt);
        fprintf(logID, "Optimization completed in %s.\n\n", elapsed);

        % --------------- CONVERSION AND SAVING ---------------
        % Try to convert output to table form; save as MAT files if error
        try
            % Convert output to table form for saving
            id = (1:nopt)';
            initial = addvars(initial, id, 'Before', 1, ...
                'NewVariableNames', "ID");
            optimal = array2table([id, optimal], ...
                'VariableNames', ["ID", knames(1,:)]);
            results = array2table([id, final_cost, flag, iter], ...
                'VariableNames', ["ID", "Cost", "Flag", "Iterations"]);
        catch ME
            % Save workspace if the variables were not saved correctly
            mat_fn = "docs/logs/" + var.sim + "_" + var.suffix + ".mat";
            save(mat_fn)

            % Throw the error after saving
            rethrow(ME)
        end

        % Which variables to save as .csv files
        savevars = struct("initial", initial, "optimal", optimal, ...
            "results", results);

    otherwise
        % Throw an error if the specified case does not exist
        eid = 'Case:doesNotExist';
        msg = 'The requested simulation does not exist.';
        error(eid, msg)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Output

% -------------------- SAVE OUTPUT --------------------
if options.Save
    % Try to save the outputs into CSV files; save as MAT files if error
    try
        % Start timer for writing to output files
        tic_save = tic;

        % Set path to output files
        filename = var.folder + "/" + outfile;

        % Write each saved variable to a separate .csv file
        varnames = names(savevars);
        for i = 1:length(varnames)
            varname = regexprep(varnames(i), "_", "-");
            filename_i = join([filename, varname], "_") + ".csv";

            % Only write the variable to the file if it is non-empty
            if ~isempty(savevars.(varnames(i)))
                writetable(savevars.(varnames(i)), filename_i)
            else
                fprintf(logID, "`%s` is empty and was not saved.\n\n", ...
                    varnames(i));
            end
        end

        % Display time elapsed for writing output
        elapsed = code_time(tic_save);
        fprintf(logID, "Writing output files completed in %s.\n\n", ...
            elapsed);

    catch ME
        % Save workspace if the variables were not saved correctly
        mat_fn = "output/logs/" + outfile + ".mat";
        save(mat_fn)

        % Throw the error after saving
        rethrow(ME)
    end
end

% Print job details and shut down parallel pool after script completion
delete(poolobj);

% Display total time elapsed for driver file
elapsed = code_time(tic_driver);
fprintf(logID, "Driver completed in %s.\n\n", elapsed);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions

% -------------------- RECEPTOR CONCENTRATION --------------------
% Vary receptor concentrations by total receptor or by ratio of 6R to 8R
function concR = receptor_conc(vary, minrecep, maxrecep, nrecep, ratio)

% Generate initial concentrations for each receptor based on how they are
% being varied for the simulations
switch vary
    % Not varying the receptor concentrations
    case "none"
        % If the receptors aren't being varied, setting the total receptor
        % and ratio to empty vectors makes this function return an empty
        % vector for concR; then, the antibody_conc() function will use the
        % experimental conditions to determine the receptor concentration
        concR = [];

    % Varying the individual receptors independently
    case {"both", "both-lowres"}
        % Set the concentrations for each receptor
        IL6R = logspace(minrecep, maxrecep, ceil(nrecep));
        IL8R = logspace(minrecep, maxrecep, ceil(nrecep));

        % Generate all combinations of the receptors
        concR = combine(IL6R, IL8R);
        concR = [concR.IL6R, concR.IL8R];

    % Vary each receptor concentration by a set percentage for the local
    % sensitivity analysis
    case "local"
        % Set the total concentration and vary each receptor by a
        % percentage of the initial value
        IL6R = [10^minrecep, 10^minrecep * 1.1]/2;
        IL8R = [10^minrecep, 10^minrecep * 1.1]/2;

        % Generate all combinations of the receptors
        concR = combine(IL6R, IL8R);
        concR = [concR.IL6R, concR.IL8R];

% Generate initial concentrations for each receptor depending on
% what is being varied - the ratio of IL6R to IL8R or the total
% receptor concentration
    otherwise
        % Set the total concentrations and ratios depending on which is
        % being varied
        totalR = logspace(minrecep, maxrecep, ceil(nrecep));

% Initialize receptor vectors
IL6R = zeros(length(totalR), length(ratio));
IL8R = zeros(length(totalR), length(ratio));

% Calculate the concentrations of each receptor based on the ratio
% and the total
for i = 1:length(totalR)
    IL6R(i,:) = (ratio./(ratio + 1)) * totalR(i);
    IL8R(i,:) = totalR(i) - IL6R(i,:);
end

% Concatenate the receptor concentrations into a single matrix
concR = [reshape(IL6R, [], 1), reshape(IL8R, [], 1)];
end

end

% -------------------- ANTIBODY CONCENTRATION --------------------
% Setup concentration tables for each antibody
function yin = antibody_conc(ab, recep, concAb, concR, types, times, setup)

% Initialize initial concentration table
yin = cell(length(ab), 1);
id = 1;

% Generate initial concentrations for each antibody
for i = 1:length(ab)

    % Select current antibodies as a string - this will correctly
    % handle if the antibody variable is a cell array or string
    % array
    abi = string(ab{i});

    % Ensure that the antibodies and receptors have consistent formatting
    [abi, recepi] = validate_ab(abi, recep{i});

    % Select the cell line for the current simulation
    uniq_recep = unique(horzcat(recepi{:}));

    % Get receptor number if it wasn't already given
    if isempty(concR)
        % "nostruct" option prevents the function from generating m and
        % p to save time
        [~, ~, cond] = setup_case(setup, uniq_recep, "nostruct");
        concRi = cond.recep;
    else
        % If the receptor numbers were given, assign the correct values
        % based on which cell line is being used
        if length(uniq_recep) == 2
            concRi = concR;
        elseif uniq_recep == "6R"
            concRi = concR(:,1);
        elseif uniq_recep == "8R"
            concRi = concR(:,2);
        end
    end

    % Add concentrations for second antibody if present
    if length(abi) > 1
        concAbi = repmat({concAb/2}, 1, 2);
    else
        concAbi = concAb;
    end

    % Use setup_input to generate the initial concentrations table
    % with the antibody and receptor names and the time periods
    yin{i} = setup_input(id, abi, recepi, ...
        types, times, abi, concAbi, recepi, concRi);

    % Increment ID
    id = max(yin{i}.ConcID) + 1;
end

% Concatenate the initial concentration table
yin = vertcat(yin{:});
end

% -------------------- SELECT SAMPLE OF INPUTS --------------------
% Benchmark simulations with a fraction of the total simulations
function [id, yin, params] = sample_inputs(id, yin, params, fraction)

% Select a random sample of simulation IDs
rng(317)
sample = randsample(id.ID, round(length(id.ID)*fraction));

% Filter the input tables to only the selected ID values
id = id(ismember(id.ID, sample), :);
yin = yin(ismember(yin.ConcID, unique(id.ConcID)), :);
params = params(...
    ismember(params.ParamID, unique(id.ParamID)), :);
end

% -------------------- RUN SIMULATIONS --------------------
% Run simulations in parfor loop
function out = run_simulations(n_split, n_message, logID, ...
    id, yin, params, output, m, p, opts)

% Start timer for simulation code
tic_sim = tic;

% Split inputs for parfor loop
[id_i, yin_i, params_i] = split_inputs(n_split, id, yin, params);

% Initialize output variable and counter
out = cell(n_split,1);

% Initialize the progress message function, create a data queue
% object to receive messages from parallel workers, and display
% progress messages throughout code execution
display_progress(n_split, n_message, tic_sim, logID)
queue = parallel.pool.DataQueue;
afterEach(queue, @display_progress)

% Run simulations
parfor i = 1:n_split
    out{i} = binding_sim(id_i{i}, yin_i{i}, params_i{i}, ...
        output, m, p, opts);

    % Ping the data queue after each simulation is complete to
    % update the progress table
    send(queue, [])
end

% Concatenate output
out = vertcat(out{:});

% Display time elapsed for simulations
elapsed = code_time(tic_sim);
fprintf(logID, "Simulations completed in %s.\n\n", elapsed);
end

% -------------------- NORMALIZE OUTPUT --------------------
% Calculate normalized output and cost
function [norm, cost, denom] = calculate_norm(norm_time, data, data_fn, ...
    basis, calc, id, yin, params, output, out, logID)

% Start timer for normalization
tic_norm = tic;

% Select only the outputs for the concentrations of antibody-receptor
% complexes that can be used to determine the normalization denominators
idx_norm = output{...
    output.Interest == "receptors" & output.Opts.Add == 1 & ...
    (output.Calc == "conc" | ...
    (output.Calc == "conct" & output.Opts.Time == norm_time)), ...
    "OutputID"};

% If there are no outputs that can be used to determine the normalization
% denominators, print a message and return to the calling function
if isempty(idx_norm)
    [norm, cost, denom] = deal(table('Size', [0 0]));
    fprintf(logID, ...
        "No suitable output calculations to use for normalization.\n\n");
    return
end

% Filter the output table
out_norm = out(out.OutputID == idx_norm, :);

% Select only the outputs for the concentrations of antibody-receptor
% complexes that can normalized with the pre-calculated normalization
% denominators; e.g., the total bound antibody is used as the basis for
% normalization, so then the binary and ternary complexes can be normalized
idx_other = output{...
    (ismember(output.Interest, ["binary", "ternary", "receptors"])) & ...
    output.Opts.Add == 1 & ...
    (output.Calc == "conc" | output.Calc == "conct"), ...
    "OutputID"};

% Remove outputs that will be used to calculate the denominators to avoid
% duplicate results
idx_other = idx_other(~ismember(idx_other, idx_norm));
% The output values will be filtered separately for each output type

% norm_output requires that the only input that is varied is the
% concentration, so each set of rate constants need to be
% normalized separately

% Select just the simulation ID and parameter ID pairs for indexing
param_idx = unique(id(:, ["ID", "ParamID"]), "rows");

% Initialize outer loop output
norms = combine(basis, calc);
norm = cell(height(norms), 1);
cost = cell(height(norms), 1);
denom = cell(height(norms), 1);

% For each normalization basis and calculation type requested
for i = 1:height(norms)

    % Set basis and calculation type
    basis_i = norms.basis(i);
    calc_i = norms.calc(i);

    % If only one data set was given use it; otherwise, select by norm type
    if length(data) == 1
        data_idx = 1;
    else
        data_idx = regexpl(data_fn, basis_i + "-Norm");
    end

    % Select corresponding data
    data_i = data(data_idx);

    % Initialize output
    n_norm = max(params.ParamID);
    norm_i = cell(n_norm, 1);
    cost_i = cell(n_norm, 1);
    denom_i = cell(n_norm, 1);

    % Calculate normalized values for each parameter set
    for j = 1:n_norm

        % Find input IDs that correspond to current parameter set
        idx = param_idx{param_idx.ParamID == j, "ID"};

        % Select the outputs corresponding to current parameter set
        out_norm_j = out_norm(ismember(out_norm.ID, idx), :);

        % Calculate normalized values using first data set, and calculate
        % diff compared to all data sets
        diff = cell(length(data_i), 1);
        for k = 1:length(data_i)
            if k == 1
                [norm_i{j}, diff{k}, denom_i{j}] = norm_output(...
                    out_norm_j, id, yin, data_i{k}, basis_i, calc_i, ...
                    norm_time, "table");

                % Add a column with the output ID to the normalized values
                norm_i{j} = addvars(norm_i{j}, ...
                    repmat(idx_norm, height(norm_i{j}), 1), ...
                    'NewVariableNames', "OutputID", 'After', 1);
                denom_i{j} = addvars(denom_i{j}, ...
                    repmat(idx_norm, height(denom_i{j}), 1), ...
                    'NewVariableNames', "OutputID", 'After', 1);

                % Use the output denominator to calculate the normalized
                % values the additional data, iterating through each of the
                % other outputs

                % Initialize output
                norm_other = cell(length(idx_other), 1);
                denom_other = cell(length(idx_other), 1);

                for a = 1:length(idx_other)
                    % Select the outputs corresponding to current parameter
                    % set and output type
                    out_other_a = out(...
                        ismember(out.ID, idx) & ...
                        ismember(out.OutputID, idx_other(a)), :);

                    [norm_other{a}, ~, denom_other{a}] = norm_output(...
                        out_other_a, id, yin, data_i{k}, ...
                        basis_i, calc_i, norm_time, "table", denom_i{j});

                    % Add a column with the output ID to the normalized
                    % values
                    norm_other{a} = addvars(norm_other{a}, ...
                        repmat(idx_other(a), height(norm_other{a}), 1), ...
                        'NewVariableNames', "OutputID", 'After', 1);
                    denom_other{a} = addvars(denom_other{a}, ...
                        repmat(idx_other(a), height(denom_other{a}), 1),...
                        'NewVariableNames', "OutputID", 'After', 1);
                end
            else
                [~, diff{k}] = norm_output(out_norm_j, id, yin, ...
                    data_i{k}, basis_i, calc_i, norm_time, "table");
            end
        end

        % Concatenate the norm tables for different outputs together
        norm_other = vertcat(norm_other{:});
        norm_i{j} = vertcat(norm_i{j}, norm_other);
        denom_other = vertcat(denom_other{:});
        denom_i{j} = vertcat(denom_i{j}, denom_other);

        % Concatenate differences and calculate cost from sum of squares
        diff = vertcat(diff{:}).Value;
        cost_i{j} = table(j, sum(diff .* diff), ...
            'VariableNames', ["ParamID", "Value"]);

        % Add ID column to denominator output
        denom_i{j} = addvars(denom_i{j}, ...
            repmat(j, height(denom_i{j}), 1), ...
            'NewVariableNames', "ParamID", 'Before', 1);
    end

    % Concatenate output
    norm_i = vertcat(norm_i{:});
    cost{i} = vertcat(cost_i{:});
    denom_i = vertcat(denom_i{:});

    % Add columns for the normalization basis and calculation type
    norm{i} = addvars(norm_i, ...
        repmat(basis_i, height(norm_i), 1), ...
        repmat(calc_i, height(norm_i), 1), ...
        'NewVariableNames', ["Basis", "Calc"]);
    denom{i} = addvars(denom_i, ...
        repmat(basis_i, height(denom_i), 1), ...
        repmat(calc_i, height(denom_i), 1), ...
        'NewVariableNames', ["Basis", "Calc"]);
end

% Concatenate output
norm = vertcat(norm{:});
cost = vertcat(cost{:});
denom = vertcat(denom{:});

% Display time elapsed for normalization
elapsed = code_time(tic_norm);
fprintf(logID, "Normalization completed in %s.\n\n", elapsed);
end

% -------------------- DISPLAY PROGRESS --------------------
% Create a formatted progress table to update as the code progresses
function display_progress(n_total, n_message, start_tic, fileID)

% Define persistent variables that will be stored in memory outside of the
% function - the previous value will be available the next time the
% function is called
persistent N N_digit N_done per_message tic0 ID

% The function is only called with three arguments when the table is being
% initialized
if nargin == 4

    % Define the base values for the persistent variables
    N = n_total;
    N_digit = ceil(log10(abs(N) * (1 + eps(N))));
    N_done = 0;
    per_message = round(N/n_message);
    tic0 = start_tic;
    ID = fileID;

    % Determine the current system time and display the start message
    [~, current] = code_time(tic0);
    fprintf(ID, ...
        "\n==========  Total Iterations: %d     " + ...
        "Start Time: %s  ==========\n\n", ...
        N, current);

else
    % Increment the counter of the number of times this function has been
    % called
    % Since the function is called every time an iteration is completed,
    % this is also the count of the number of done iterations
    N_done = N_done + 1;

    % Display a new row of the progress table every time a certain number
    % of iterations is compleeted (based on the number of messages
    % requested)
    if mod(N_done, per_message) == 0

        % Determine the amount of time elapsed since the start tic and the
        % current system time for the progress table
        [elapsed, current] = code_time(tic0);

        % The first message displayed needs to include the table header
        if N_done == per_message

            % Table header row and dividing line
            fprintf(ID, ...
                "     Completed%*s|" + ...
                "   Elapsed    |" + ...
                "   Time\n", ...
                N_digit+2, "");
            fprintf(ID, ...
                "   %s\n", strjoin(repmat("-", N_digit+42, 1), ""));

            % First data row in the progress table
            fprintf(ID, ...
                "     %*d ( %2.0f%% )   |" + ...
                "   %s   |" + ...
                "   %s\n", ...
                N_digit, N_done, N_done/N * 100, ...
                elapsed, current);

            % The final message displayed has one fewer space to accommodate
            % the extra digit in 100%
        elseif N_done == N

            % Final data row in the progress table
            fprintf(ID, ...
                "     %*d ( %2.0f%% )  |" + ...
                "   %s   |" + ...
                "   %s\n", ...
                N_digit, N_done, N_done/N * 100, ...
                elapsed, current);

            % Dividing line under the progress table
            fprintf(ID, ...
                "   %s\n\n", strjoin(repmat("-", N_digit+42, 1), ""));

            % All other table rows
        else
            % Standard data row in the progress table
            fprintf(ID, ...
                "     %*d ( %2.0f%% )   |" + ...
                "   %s   |" + ...
                "   %s\n", ...
                N_digit, N_done, N_done/N * 100, ...
                elapsed, current);
        end
    end
end
end
