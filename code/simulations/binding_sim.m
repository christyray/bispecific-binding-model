function out = binding_sim(id, yin, params, output, m, p, opts)

% BINDING_SIM	Perform model simulations for several input values.
%
%	Runs the model simulations for each set of input values and calculates
%	the specified output values for molecules of interest. Calls
%	BINDING_MAIN for the ODE solver.
%
%	USAGE:
%		OUT = BINDING_SIM(ID, YIN, PARAMS, OUTPUT, M, P, OPTS)
%
%	INPUT:
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
%       PARAMS = table with columns for the simulation ID, parameter name,
%       and parameter value for each parameter in the simulation; generated
%       by SETUP_INPUT
%
%       OUTPUT = table with all combinations of antibody/receptor combos,
%       molecules of interest, calculation types, and options for
%       OUTPUT_CALC; created by SETUP_OUTPUT
%
%       M = structure containing the molecules included in the model and
%       their assigned numbers, correct naming for species is antibody name
%       or antibody receptor complex name (with components separated by
%       underscores); defined by SETUP_CASE
%
%       P = structure containing parameters for the model, correct naming
%       for parameters is "kon_" or "koff_", followed by the antibody name
%       or antibody receptor complex name (with components separated by
%       underscores), ending with the receptor being bound to; e.g.,
%       kon_BS1_6R_8R; defined by SETUP_CASE
%
%       OPTS = optional cell array containing options for the ODE solver in
%       BINDING_MAIN; defaults to the default options in BINDING_MAIN
%
%	OUTPUT:
%		OUT = long form table with calculated simulation output; contains
%       columns for input and output ID (corresponding to INPUT and
%       OUTPUT), which calculation was performed, and the time points,
%       species name, and value for each output point
%
%	NOTES:
%		OUT can be directly exported to a .csv file to be used in R for
%		plotting or can be analyzed further in MATLAB. Specific outputs
%		corresponding to a single input and output calculation can be found
%		from the combination of the ID and OutputID columns.
%
%	See also BINDING_MAIN, SETUP_INPUT, SETUP_ID, SETUP_OUTPUT, SETUP_CASE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    id table
    yin table
    params table
    output table
    m struct
    p struct
    opts cell = cell(0)
end

% Throw an error if there is no value for p.alpha
if p.alpha == 0
    eid = 'Parameter:notSet';
    msg = 'The alpha parameter has not been set.';
    throwAsCaller(MException(eid,msg))
end

% If the inputs are all empty, return to the calling function with an empty
% table for the outputs
% This can happen when there are fewer than 100 simulations but the n_split
% is still set to 100 for the parfor loop
if all([isempty(id), isempty(yin), isempty(params)])
    varname = ["ID", "OutputID", "Calc", "Time", "Species", "Value"];
    vartype = ["double", "double", "string", "double", "string", "double"];
    out = table('Size', [0 length(varname)], 'VariableTypes', vartype, ...
        'VariableNames', varname);
    return
elseif any([isempty(id), isempty(yin), isempty(params)])
    eid = 'Input:empty';
    msg = 'One or more of the input tables are empty.';
    error(eid,msg)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup Simulations

% Number of simulations to perform
sims = unique(id.ID);           % For parfor loops with selected IDs
calcs = max(output.OutputID);
n = length(sims) * calcs;

% Number of species in system
nmol = length(names(m));

% Store initial alpha value for reinitialization
alpha = p.alpha;

% Load in the bound and single complexes to select molecules of interest
load complexes.mat bound singles

% Initialize counter
ni = 1;

% Initialize output
out = cell(n, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Simulations

% Loop through each set of input values
for i = 1:length(sims)

    % Select IDs for current initial concentations and parameter values
    sim = sims(i);
    concID = id{id.ID == sim, "ConcID"};
    paramID = id{id.ID == sim, "ParamID"};

    % Select current data
    in_i = yin(yin.ConcID == concID, :);
    ab = strsplit(in_i.Ab(1), "_");         % Split into separate Ab
    recep = strsplit(in_i.Recep(1), "_")';  % Split by different Ab first
    recep = mat2cell(recep, ones(length(recep), 1));    % Convert to cell
    recep = cellfun(@(x) strsplit(x, "-"), recep, 'UniformOutput', false);
    % Split antibodies for a single receptor

    % Assign time period information, may have duplicate names but
    % different times
    times = unique([in_i.Start, in_i.End]);
    ytype = unique([in_i.Period, in_i.Start], "rows", "stable");
    ytype = ytype(:,1);

    % Filter the input table to just the input concentrations
    in_conc = unique(in_i(:, ["Period", "Species", "Conc"]), ...
        "rows", "stable");

    % Initialize input concentration array
    ntimes = length(ytype);
    yin_i = [zeros(1,nmol); inf(ntimes-1, nmol)];

    % Assign input concentrations
    for j = 1:height(in_conc)
        % Which time period the value is for
        period = in_conc.Period(j) == ytype;

        % Assign the value to the input concentration matrix
        yin_i(period, m.(in_conc.Species(j))) = in_conc.Conc(j);
    end

    % Filter the input table to just the input parameters
    in_params = params(params.ParamID == paramID, :);

    % Replace the "Ab" in the parameter names with the specific antibody
    ab_rep = append("(", join(ab, "|"), ")");
    in_params.Param = regexprep(in_params.Param, "Ab", ab_rep);

    % Re-initialize p to erase any previous values
    p = struct_init(p, 0);
    p.alpha = alpha;

    % Assign parameter values
    for j = 1:height(in_params)
        p = field_assign(p, in_params.Param(j), in_params.Value(j));
    end

    % Solve the ODE solver
    [T, Y, net, bal] = binding_main(times, yin_i, ytype, m, p, opts{:});

    % Calculate outputs
    for j = 1:calcs

        % Select current data
        out_j = output(j,:);

        % Total bound receptor is a special case and calculated separately
        if out_j.Interest ~= "complexes"
            % Determine specific molecules of interest from selection
            interest = select_molecules(ab, recep, out_j.Interest, ...
                bound, singles);

            % Calculate output
            values = calc_output(T, Y, net, bal, m, p, ...
                singles, interest, out_j.Calc, out_j.Opts);
        else
            % calc_complexes combines select_molecules and calc_output for
            % the special case of total bound receptor
            values = calc_complexes(T, Y, net, bal, m, p, ab, recep, ...
                bound, singles, out_j.Calc, out_j.Opts);
        end

        % Generate an ID column for the output table
        id_col = table(sim, j, out_j.Calc, ...
            'VariableNames', ["ID", "OutputID", "Calc"]);
        id_col = repmat(id_col, height(values), 1);

        % Concatenate output
        out{ni} = [id_col, values];
        ni = ni + 1;
    end
end

% Convert output cell array to a long form table
out = vertcat(out{:});
