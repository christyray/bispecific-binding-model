function varfile = setup_variables(options)

% SETUP_VARIABLES	Generate simulation variables for the driver file.
%
%	Generate all of the input variables for BINDING_DRIVER, including which
%	simulation case to run, which antibodies and receptors to simulate, the
%	normalization to use and the experimental data files to import, the
%	outputs to calculate, the parameter values to simulate, and parameters
%	to optimize and their initial guesses.
%
%	USAGE:
%		VARFILE = SETUP_VARIABLES(OPTIONS)
%
%	INPUT:
%		OPTIONS = named arguments for which preset values to use for
%		certain variables; each named variable defaults to "other" which
%		allows for manual specification of each variable or "none" which
%		sets no value for the variable
%
%	OUTPUT:
%		VARFILE = path to file where the generated variables are saved;
%		follows the same naming syntax as the simulation results
%
%	NOTES:
%		The OPTIONS include common presets for many of the simulation
%		variables. If no value is given for a particular argument, it will
%		default to "other", which is the manual specification for that
%		variable. This allows for options other than the presets to easily
%		be used.
%
%	See also BINDING_DRIVER, BINDING_SIM.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    options.Name string = "sim"
    options.Case string {mustBeMember(options.Case, ...
        ["sim", "opt", "global", "local", "other"])} = "other"
    options.Antibodies string {mustBeMember(options.Antibodies, ...
        ["binding", "negs", "compare", "all", "bs1", "toci-h2", ...
        "compare-double", "other"])} = "other"
    options.Periods string {mustBeMember(options.Periods, ...
        ["experimental", "long-washout", "long-initial", "initial", ...
        "other"])} = "other"
    options.Conc string {mustBeMember(options.Conc, ...
        ["standard", "lowres", "highres", "high", "small", "local", ...
        "other", "none"])} = "none"
    options.Receptors string {mustBeMember(options.Receptors, ...
        ["total", "lowres", "ratio", "both", "both-lowres", "local", ...
        "other", "none"])} = "none"
    options.Params string {mustBeMember(options.Params, ...
        ["optimal", "import", "monovalent", "other", "none"])} = "none"
    options.Outputs string {mustBeMember(options.Outputs, ...
        ["experimental", "time", "key", "base", "endpoint", "local", ...
        "other", "none"])} = "none"
    options.Norm string {mustBeMember(options.Norm, ...
        ["bs1-data", "ab-data", "bs1-max", "ab-max", "all", "other", ...
        "none"])} = "none"
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Options for Driver

% Suppress warning about variables not being used
%#ok<*NASGU>

%% -------------------- FILE PARAMETERS --------------------

% Specific case to run
% Options: ["sim", "opt", "global", "local", "bivariate"]
if options.Case == "other"
    % Option to manually specify the simulation case to run
    options.Case = "opt";
end
sim = options.Case;

% Set name for output files
if options.Name == "optimization" || options.Name == "compare-opt"
    name = join([options.Name, options.Norm], "-");
else
    name = options.Name;
end

% Folder to save output in
folder = "output";

%% -------------------- SIMULATION INPUTS --------------------

% Initialize the ODE solver options to add onto as necessary
odeopts = cell(0);

% Antibodies and receptor partners to include in optimization and
% simulation
switch options.Antibodies
    % Antibodies and their possible binding partners
    case "binding"
        ab = ["Toci", "Toci", ...
            "H2", "H2", ...
            "BS1", "BS1", "BS1"];
        recep = {["6R", "6R"], ["6R", "6R", "8R"], ...
            ["8R", "8R"], ["8R", "8R", "6R"], ...
            "6R", "8R", ["6R", "8R"]};
    % Antibodies and possible binding partners, include non-binding cells
    case "negs"
        ab = ["Toci", "Toci", "Toci",...
            "H2", "H2", "H2", ...
            "BS1", "BS1", "BS1"];
        recep = {["6R", "6R"], "8R", ["6R", "6R", "8R"], ...
            "6R", ["8R", "8R"], ["8R", "8R", "6R"], ...
            "6R", "8R", ["6R", "8R"]};
    % Include combination of Toci and 10H2 for comparison
    case "compare"
        ab = {"Toci", "H2", "BS1", ["Toci", "H2"]};
        recep = {["6R", "6R", "8R"], ["8R", "8R", "6R"], ...
            ["6R", "8R"], ...
            ["6R", "6R"; "8R", "8R"]};
    % Antibodies in all cell lines, including combination of Toci and 10H2
    case "all"
        ab = {"Toci", "Toci", "Toci", ...
            "H2", "H2", "H2", ...
            "BS1", "BS1", "BS1", ...
            ["Toci", "H2"], ["Toci", "H2"], ["Toci", "H2"]};
        recep = {["6R", "6R"], ["8R", "8R"], ["6R", "6R", "8R"], ...
            ["6R", "6R"], ["8R", "8R"], ["8R", "8R", "6R"], ...
            "6R", "8R", ["6R", "8R"], ...
            ["6R", "6R"; "6R", "6R"], ["8R", "8R"; "8R", "8R"], ...
            ["6R", "6R"; "8R", "8R"]};
    % Simulation with only BS1 in the double-positive cell line
    case "bs1"
        ab = "BS1";
        recep = {["6R", "8R"]};
    % Simulation with only tocilizumab and 10H2 in the double-positive cell
    % line
    case "toci-h2"
        ab = {["Toci", "H2"]};
        recep = {["6R", "6R"; "8R", "8R"]};
    % Simulation with only BS1 and tocilizumab + 10H2 in the
    % double-positive cell line
    case "compare-double"
        ab = {"BS1", ["Toci" "H2"]};
        recep = {["6R", "8R"], ["6R", "6R"; "8R", "8R"]};
    % Option to manually specify which antibodies and receptors to use
    case "other"
        ab = "BS1";
        recep = {["6R", "8R"]};
end

% Experiment performed
setup = "flow";

% Which time points and time periods to use
switch options.Periods
    % Use the time points and period types from the experimental setup
    case "experimental"
        % Use setup_case to determine the experimental conditions
        [~, ~, cond] = setup_case(setup, recep{1}, "nostruct");
        period_types = cond.ytype;
        period_times = ...
            [0 cond.primary cond.primary + cond.secondary] * 3600;
    % Use the period types from the experimental setup with an extended
    % washout period
    case "long-washout"
        % Use setup_case to determine the experimental conditions
        [~, ~, cond] = setup_case(setup, recep{1}, "nostruct");
        period_types = cond.ytype;
        period_times = [0 cond.primary cond.primary * 5] * 3600;

        % Extend the time step for the washout period to limit the number
        % of output points
        odeopts = [odeopts, {"tstep", [30, 300]}];
    % Run the simulation for the initial period only, longer than the
    % experimental time periods
    case "long-initial"
        % Single period, maximum time of 24 hours; using the "addition"
        % period to use different time steps for the beginning and end of
        % the simulation
        period_types = ["initial", "addition"];
        period_times = [0 2 24] * 3600;

        % Extend the time step to limit the number of output points
        odeopts = [odeopts, {"tstep", [30, 600]}];
    % Run the simulation for the initial period only, with no extended time
    % or washout
    case "initial"
        period_types = "initial";
        period_times = [0 7200];
    % Option to manually specify time points and period types
    case "other"
        period_types = "initial";
        period_times = [0 7200];
end

% Initial antibody concentrations to use
% stepconc is the number of concentrations to use between each order of
% magnitude
% minconc and maxconc set the range of magnitudes for the concentrations
% expconc is the option to include the exact experimental concentrations
switch options.Conc
    % Use the range of concentrations used in the experimental data, with
    % an extra order of magnitude on either side
    case "standard"
        % Set the options for generating the range of concentrations
        stepconc = 10;
        minconc = -3;
        maxconc = 4;
        expconc = true;

    % Use the same range of concentrations but with fewer values, for
    % simulations that have a large number of variables
    case "lowres"
        % Set the options for generating the range of concentrations
        stepconc = 1;
        minconc = -3;
        maxconc = 4;
        expconc = false;

    % Use the same range of concentrations but with more values, for
    % simulations where the antibody concentration is the main variable
    case "highres"
        % Set the options for generating the range of concentrations
        stepconc = 30;
        minconc = -3;
        maxconc = 4;
        expconc = false;

    % Extend the concentration range to include higher concentrations to
    % model the entire range of the binding curve
    case "high"
        % Set the options for generating the range of concentrations
        stepconc = 10;
        minconc = -3;
        maxconc = 7;
        expconc = true;

    % Use a smaller range of concentration values for simulations that are
    % primarily varying other inputs
    case "small"
        % Set the options for generating the range of concentrations
        stepconc = 1;
        minconc = -1;
        maxconc = 3;
        expconc = false;

    % Vary the concentrations by a set percentage for the local sensitivity
    % analysis
    case "local"
        % Set the options for generating the percentage change in
        % concentration
        stepconc = 1;
        minconc = 1;
        maxconc = log10(10^minconc * 1.1);
        expconc = false;

    % Option to manually specify initial antibody concentrations
    case "other"
        stepconc = 10;
        minconc = -3;
        maxconc = 4;
        expconc = false;
end

% Determine number of initial antibody concentrations from the options
if options.Conc == "none"
    nconc = []; minconc = []; maxconc = []; stepconc = [];
else
    nconc = (maxconc - minconc) * stepconc + 1; 
end

% How to vary the receptor concentrations - total receptors, ratio of one
% receptor to the other, "both" to vary each receptor independently,
% "local" to vary the concentrations by a set percentage for the local
% sensitivity analysis, or "none" to use the experimental receptor
% concentrations
switch options.Receptors
    % Vary total receptor concentration
    case "total"
        % Set the options for generating the range of concentrations
        steprecep = 20;
        minrecep = 2;
        maxrecep = 7;
        ratio = 1;

    % Vary total receptor concentration using a limited number of values
    case "lowres"
        % Set the options for generating the range of concentrations
        steprecep = 1;
        minrecep = 2;
        maxrecep = 7;
        ratio = [2 1 1/2];

    % Vary ratio of IL6R to IL8R
    case "ratio"
        % Set the options for generating the range of concentrations
        steprecep = 0.5;
        minrecep = 3;
        maxrecep = 7;
        ratio = [10, 5, 2, 1.5, 1, 1/1.5, 1/2, 1/5, 1/10];

    % Varying the individual receptors independently
    case "both"
        % Set the options for generating the range of concentrations
        steprecep = 30;
        minrecep = 2;
        maxrecep = 7;
        ratio = 1;

    % Varying the individual receptors independently using a limited number
    % of values
    case "both-lowres"
        % Set the options for generating the range of concentrations
        steprecep = 1;
        minrecep = 2;
        maxrecep = 7;
        ratio = 1;

    % Varying receptor concentrations for local sensitivity analysis
    case "local"
        steprecep = 1;
        minrecep = 5;
        maxrecep = 5;
        ratio = 1;

    % Use the experimental receptor concentrations
    case "none"
        % Set everything to 0 because these values will not be used
        steprecep = 0;
        minrecep = 0;
        maxrecep = 0;
        ratio = 0;

    % Option to manually specify receptor concentrations
    case "other"
        steprecep = 28;
        minrecep = 2;
        maxrecep = 7;
        ratio = 1;
end

% Determine number of initial receptor concentrations from the options
if options.Receptors == "none"
    nrecep = []; minrecep = []; maxrecep = []; steprecep = [];
else
    nrecep = (maxrecep - minrecep) * steprecep + 1;
end

% Number of parameter values to simulate
nkvalues = 17;  % Number of different parameter values when varied
nmag = 2;       % Orders of magnitude to vary parameter values over

% Baseline parameter values to simulate
switch options.Params
    % Use the parameter values from the optimization of kon and koff
    case "optimal"
        % Baseline values
        kvalues = [5.91965e-6, 9.02608e-6, 8.10825e-8, 1.23632e-7, ...
            5.60921e-5, 6.37590e-5];
        knames = ["kon_Ab_6R$", "kon_Ab_8R$", "kon_Ab.*_[68]R_6R$", ...
            "kon_Ab.*_[68]R_8R$", "koff_Ab.*_6R$", "koff_Ab.*_8R$"];

    % Import optimized parameter values from a saved file
    case "import"
        % Import the optimized k values
        opt_fn = "output/optimal_parameters.csv";
        optimal = readtable(opt_fn);

        % Filter the columns to just the normalizations matching the
        % simulation options
        if ismember(options.Norm, ...
                ["bs1-data", "bs1-max", "ab-data", "ab-max"])
            norm_filter = strsplit(options.Norm, "-");
            optimal = optimal(strcmp(optimal.Basis, norm_filter(1)) & ...
                strcmp(optimal.Calc, norm_filter(2)), :);
        end

        % Select the columns containing k values
        col_names = names(optimal);
        cols = regexpl(col_names, "^(Param|ko).*");
        kvalues = optimal{:, cols};

        % Determine the names corresponding to the kvalues
        knames = col_names(cols);
        knames = regexprep(knames, "^Param_", "");
        knames = replace_names(knames);

    % Simulations where only the first binding step is allowed
    case "monovalent"
        % Baseline values
        kvalues = [5.91965e-6, 9.02608e-6, 8.10825e-8, 1.23632e-7, ...
            5.60921e-5, 6.37590e-5];
        knames = ["kon_Ab_6R$", "kon_Ab_8R$", "kon_Ab.*_[68]R_6R$", ...
            "kon_Ab.*_[68]R_8R$", "koff_Ab.*_6R$", "koff_Ab.*_8R$"];

        % Values for only the first binding step
        kvalues = [kvalues; ...
            kvalues(1), kvalues(2), 0, 0, kvalues(5), kvalues(6)];

    % Option to manually specify baseline parameter values
    case "other"
        kvalues = [5.91965e-6, 9.02608e-6, 4.915e-12, 7.4943e-12, ...
            5.60921e-5, 6.37590e-5];
        knames = ["kon_Ab_6R$", "kon_Ab_8R$", "kon_Ab.*_[68]R_6R$", ...
            "kon_Ab.*_[68]R_8R$", "koff_Ab.*_6R$", "koff_Ab.*_8R$"];
end

% Which parameter values to vary in the bivariate analysis
vary = [1, 2];  % Indices used to select from knames and kvalues

%% -------------------- OUTPUT CALCULATIONS --------------------

% Specify which output calculations should be performed for each simulation
switch options.Outputs
    % Output at the end time of the experimental data
    case "experimental"
        interest = ["receptors", "binary", "ternary", "complexes"];
        value = "conct";
        add = 1;
        time = 8100;
    % Output the concentrations at each time point
    case "time"
        interest = ["receptors", "binary", "ternary"];
        value = "conc";
        add = 1;
        time = 0;
    % Main outputs for sensitivity analysis
    case "key"
        interest = "present";
        value = ["end", "peak"];
        add = 0;
        time = 0;
    % Main outputs for the baseline simulations
    case "base"
        interest = "present";
        value = "conc";
        add = 0;
        time = 0;
    % Output all concentrations at final time point of simulation
    case "endpoint"
        interest = "present";
        value = "end";
        add = 0;
        time = 0;
    % Main outputs for local sensitivity analysis
    case "local"
        interest = ["present", "complexes"];
        value = ["conct", "peak", "auc"];
        add = 0;
        time = 7200;
    % Option to manually specify what outputs to calculate
    case "other"
        interest = "receptors";
        value = "end";
        add = [0 1];
        time = 0;
end

% Specify which normalization basis and calculation type to use and which
% data file(s) to use for simulations
switch options.Norm
    % Normalize to BS1 binding using concentrations from the data
    case "bs1-data"
        basis = "BS1";
        calc = "data";
        data_fn = "Flow-MFI_BS1-Norm";
    % Normalize each antibody to itself using concentrations from the data
    case "ab-data"
        basis = "Ab";
        calc = "data";
        data_fn = "Flow-MFI_Ab-Norm";
    % Normalize to BS1 binding using the max concentration
    case "bs1-max"
        basis = "BS1";
        calc = "max";
        data_fn = "Flow-MFI_BS1-Norm";
    % Normalize each antibody to itself the max concentration
    case "ab-max"
        basis = "Ab";
        calc = "max";
        data_fn = "Flow-MFI_Ab-Norm";
    % Use each normalization calculation for the results
    case "all"
        basis = ["Ab", "BS1"];
        calc = ["data", "max"];
        data_fn = ["Flow-MFI_BS1-Norm", "Flow-MFI_Ab-Norm"];
    % If not performing the normalization, do not set basis and calc
    case "none"
        basis = strings(0);
        calc = strings(0);
        data_fn = strings(0);
    % Option to manually specify which normalization options to use
    case "other"
        basis = "BS1";
        calc = "data";
        data_fn = "Flow-MFI_BS1-Norm";
end

%% -------------------- IMPORT DATA --------------------

% Import data
if isempty(data_fn)
    data = [];
else
    data = cell(length(data_fn), 1);
    for i = 1:length(data_fn); data{i} = import_data(data_fn(i)); end
end

%% -------------------- OPTIMIZATION PARAMETERS --------------------

% Set names and initial guesses for optimized parameters
if sim == "opt"
    % Set names of parameters being optimized
    kopt_names = ["kon_Ab_6R$", "kon_Ab_8R$", ...
        "kon_Ab.*_[68]R_6R$", "koff_Ab.*_6R$", "koff_Ab.*_8R$"];
    initial_names = ["initial_6R", "initial_8R", ...
        "initial_6R_prime", "initial_off6R", "initial_off8R"];

    % Set initial guesses and bounds for parameters being fit
    ninitial = 300;

    % Set minimum power and range for the guesses for each parameter
    minpower = repmat([-9, -9, -13, -7, -7], ninitial, 1);
    rangepower = repmat([7, 7, 8, 5, 5], ninitial, 1);

    % Use Latin Hypercube sampling to generate the initial guesses for
    % each parameter, use rng to set consistent guesses
    rng(138)
    initial = lhsdesign(ninitial, size(minpower, 2));
    initial = (initial .* rangepower) + minpower;
    initial = 10.^initial;

    % Set lower and upper bounds for optimization
    lb = round_signif(min(initial), 1, "floor") * 1e-2;
    ub = round_signif(max(initial), 1, "ceil") * 1e2;

    % Convert initial guesses into a table and name the columns based
    % on the parameter names
    initial = array2table(initial, 'VariableNames', initial_names);

    % Set values and names for parameters that are not being optimized
    kin = [];
    kin_names = strings(0);

    % Determine total number of input parameters for optimization
    nparams = width(initial) + length(kin) + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Driver Variables

% Save driver variables and options using the same file name syntax as the
% simulation results
d = string(datetime("now", "Format", "yyyy-MM-dd"));
varfile = folder + "/" + join([d, name, "vars"], "_") + ".mat";

% Save the variables necessary for the specified case in a MAT file and
% return the file name for loading

% Variables necessary for all cases
vars = ["options", "sim", "name", "folder", "ab", "recep", "setup"];

% Variables only necessary for specific cases
switch sim
    case "sim"
        vars = [vars, "period_types", "period_times", ...
            "nconc", "minconc", "maxconc", "expconc", ...
            "nrecep", "minrecep", "maxrecep", "ratio", ...
            "kvalues", "knames", "odeopts", ...
            "interest", "value", "add", "time", ...
            "data", "data_fn", "basis", "calc"];
    case "opt"
        vars = [vars, "initial", "kopt_names", "lb", "ub", ...
            "kin", "kin_names", "nparams", "data", "data_fn", ...
            "basis", "calc"];
    case {"local", "global"}
        vars = [vars, "period_types", "period_times", ...
            "nconc", "minconc", "maxconc", "expconc", ...
            "nrecep", "minrecep", "maxrecep", "ratio", ...
            "kvalues", "knames", "nkvalues", "nmag", ...
            "odeopts", "interest", "value", "add", "time", ...
            "data", "data_fn", "basis", "calc"];
    otherwise
        % Throw an error if the specified case does not exist
        eid = 'Case:doesNotExist';
        msg = 'The requested simulation does not exist.';
        error(eid, msg)
end

% Brackets and colon are necessary to process each element of the vars
% vector as a separate text string - save requires individual text scalars
% for each variable to be saved
save(varfile, vars{:})

end
