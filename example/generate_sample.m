function generate_sample()

% GENERATE_SAMPLE	Create a MAT file of sample inputs for model functions.
%
%	Generates a MAT file with all of the variables necessary to run
%	BINDING_SIM, BINDING_MAIN, BINDING_COST, BINDING_OPT, and NORM_OUTPUT. 
%   Used for sharing the project code and making it easier to illustrate 
%   how the model works.
%
%	USAGE:
%		GENERATE_SAMPLE()
%
%	NOTES:
%		The sample inputs are saved to "docs/sample.mat". Note that this
%		function does take about 10 to 20 seconds to run because it
%       performs 55 model simulations to generate the `out` table.
%
%       Use the command `load example/sample.mat` in the Command Window to
%       load the saved variables into the workspace.
%
%	See also BINDING_SIM, BINDING_COST, BINDING_OPT, NORM_OUTPUT.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Inputs for binding_sim()

% binding_sim() runs the model simulations for a provided set of initial
% conditions and parameter values and calculates different outputs for
% specified molecules of interest

% USAGE: OUT = BINDING_SIM(ID, YIN, PARAMS, OUTPUT, M, P, OPTS)

%% --------------- SIMULATION CONDITIONS ---------------

% Input to setup_case - determines which experimental conditions are used
setup = "flow";

% Antibodies and receptors present in the system - used for creating the
% initial conditions and filling in the parameter names
ab = ["Toci", "H2", "BS1", "BS1", "BS1"];
recep = {["6R", "6R"], ["8R", "8R"], "6R", "8R", ["6R", "8R"]};

% Create the molecule and parameter structures - uses the experimental
% conditions to determine the value of the unit conversion parameter
% `alpha`, generates the molecule and parameter names based on the
% `antibodies.txt` file
[m, p] = setup_case(setup);

% Set antibody concentrations for simulations - the rows correspond to
% different simulations and the columns correspond to different time
% periods in each simulation

% For the sample inputs, I am using the same concentrations that are
% present in the experimental data. The second period is a washout, so the
% unbound antibody concentration is set to zero
conc_exp = [-10.7712, -10.2941, -9.81697, -9.33985, -8.86273, -8.38561, ...
    -7.90849, -7.43136, -6.95424, -6.47712, -6];
conc_exp = 10 .^ conc_exp * 1e9;
concAb = [conc_exp', zeros(length(conc_exp), 1)];    % Units = nM

% Set time period types and length of each time period - the different
% available types are "initial" for the first period where everything is
% input into the system, "addition" for adding more of some species, and
% "washout" for removing some of some species

% The time points give the starting time, the boundaries between time
% periods, and the final time (so there will always be # of periods + 1
% time points)
times = [0 2 4] * 3600;     % Units = s
periods = ["initial", "washout"];

%% --------------- INITIAL CONCENTRATIONS ---------------

% Initialize cell array to hold concentrations tables
yin = cell(length(ab), 1);
id = 1;

% The function to generate the initial concentrations table,
% `setup_input()`, requires only one set of antibodies and receptors be
% passed in at a time (this is because it simulations can include multiple
% antibodies).

% Here, I iterate through the different antibodies to simulate and generate
% a separate table for each. The starting ID number is incremented for each
% iteration so the tables can be concatenated together.
for i = 1:length(ab)

    % Select current antibodies as a string - this will correctly
    % handle if the antibody variable is a cell array or string
    % array
    abi = string(ab{i});

    % Determine the concentration of each receptor depending on the cell
    % type - I am optimizing to experiments in IL6R+, IL8R+, and
    % IL6R+/IL8R+ cells, and each cell line expresses a different
    % concentration of receptors; receptor concentrations = # / cell

    % The "nostruct" option prevents the function from generating the m and
    % p structures to save time. The `cond` structure contains information
    % about the experimental conditions, including the number of cells, the
    % volume, and the receptor concentration
    [~, ~, cond] = setup_case(setup, recep{i}, "nostruct");

    % If two antibodies are present, the concentrations for both antibodies
    % need to be passed into the `setup_input` function as matrices in a
    % cell array (that was the easiest way to distinguish that the values
    % correspond to separate species). Here, I am holding the total
    % antibody concentration constant
    if length(abi) > 1
        concAbi = repmat({concAb/2}, 1, 2);
    else
        concAbi = concAb;
    end

    % `setup_input` generates a table with columns for the simulation ID,
    % the antibodies and receptors present in the simulation, the time
    % period type and starting and ending times, and each species and its
    % concentration at the start of each time period. Any species that are
    % not specified are assumed to be unchanged (e.g., the receptors do not
    % have values for the "washout" period because they are not washed out)
    yin{i} = setup_input(id, abi, recep{i}, ...
        periods, times, abi, concAbi, recep{i}, cond.recep);

    % Increment starting number for the simulation ID
    id = max(yin{i}.ConcID) + 1;
end

% Concatenate the separate initial concentration tables for each antibody
% into a single table
yin = vertcat(yin{:});

%% --------------- PARAMETERS ---------------

% Specify the regular expression patterns that will select the correct
% parameter fields from the `p` structure. The "Ab" in each pattern will be
% replaced with the specific antibody present in each simulation, and the
% pattern will be used to match the right fields in `p` and assign the
% corresponding parameter values
knames = ["kon_Ab_6R$", "kon_Ab_8R$", "kon_Ab.*_[68]R_6R$", ...
    "kon_Ab.*_[68]R_8R$", "koff_Ab"];
% The patterns work out to the on rates for binding to IL6R, binding to
% IL8R, binding to IL6R after already being bound to a receptor, binding to
% IL8R after already being bound to a receptor, and then the off rates for
% all binding steps

% Specify the parameter values corresponding to the parameter names
kvalues = [6.0248e-6, 8.9807e-6, 8.6016e-8, 1.2822e-7, 6.2109e-5];

% When `setup_input` is only passed three arguments, it generates the
% parameter table with columns for the simulation ID number, the parameter
% name pattern, and the parameter value
params = setup_input(1, knames, kvalues);

% The IDs from the initial concentrations table and from the parameters
% table are combined into a single table that gives a unique ID value for
% each simulation. This structure allows the simulations to vary which
% antibodies and receptors are present, what time periods are used, the
% length of the time periods, the concentrations of each species in each
% time period, and the parameter values used, all in one set of simulations
id = setup_id(yin, params);

% Columns in the BINDING_SIM input tables:
%
%   ---------- ID ----------
%   ID = a unique number identifying each separate simulation; each ID
%   value represents a unique combination of initial concentrations and
%   parameter values for each simulation
%
%   CONCID = the ID value from the `yin` table; identifies which initial
%   concentration set is used for the simulation
%
%   PARAMID = the ID value from the `params` table; identifies which
%   parameter set is used for the simulation
%
%   ---------- YIN ----------
%   CONCID = the ID value for each set of initial concentrations
%
%   AB = a string containing which antibody or antibodies are present in
%   the system during the simulation; different antibodies are separated
%   with underscores
%
%   RECEP = a string containing which antibodies are present in the system
%   during the simulation; multiple of the same antibody means the antibody
%   can bind multivalently to that receptor; receptors for the same
%   antibody are separated with dashes, receptors for different antibodies
%   are separated with underscores
%
%   PERIOD = the name of the time period that the concentrations correspond
%   to; for "initial" and "washout", the concentration of the species is
%   set to that value, for "addition" the value is added to the existing
%   concentration of the species (relevant for the mole balance)
%
%   START = the starting time of the current period in seconds
%
%   END = the ending time of the current period in seconds
%
%   SPECIES = the name of the molecule in the `m` structure that the
%   concentration corresponds to
%
%   CONC = the concentration of the molecule given in the SPECIES column
%   for the current time period; the antibody concentrations are given in
%   nM and the receptor concentrations are given in # / cell
%
%   ---------- PARAMS ----------
%   PARAMID = the ID value for each set of parameter values
%
%   PARAM = the name or a regular expression pattern identifying the
%   parameter in the `p` structure that the value corresponds to
%
%   VALUE = the value of the parameter given in the PARAM column

%% --------------- OUTPUTS ---------------

% Specify which molecules of interest to calculate outputs for; "inputs" is
% all unbound antibodies and receptors, "receptors" is all
% antibody-receptor complexes, "binary" is all complexes of an antibody
% with one receptor, "ternary" is all complexes of an antibody with two
% receptors, "complexes" is the total bound receptor - equal to binary +
% 2*ternary, "present" is all molecules present in the system (so only the
% antibodies and receptors that had initial concentrations assigned), and
% "all" is all molecules in the model
interest = ["inputs", "receptors", "binary", "ternary", "complexes"];
% `select_molecules()` will return the molecules of interest based on this
% argument

% Specfiy which output values to return; "conc" is the concentration at
% each time point, "conct" is the concentration at specified time point(s)
% (given by `time`), "peak" is the maximum concentration, "end" is the
% concentration at the final time point, "auc" is the area under the
% concentration-time curve, "dydt" is the derivative at each time point,
% "net" is the net amount of each species added into the system at each
% time point, and "bal" is the mole balance at each time point
value = ["conc", "conct", "peak", "auc"];

% `add` is a binary 1 or 0 for whether the calculation should return the
% sum of all molecules of interest or report each molecule separately; this
% is only used where it is reasonable to do so (e.g., `add` will not be set
% to 1 when the molecules of interest have different units, like unbound
% antibodies and receptors)
add = [0 1];

% `time` is the time point to return the concentration when the "conct"
% calculation is used; it is only relevant for "conct" and is ignored
% otherwise
time = 7200;

% `setup_output` generates a table with all combinations of the output
% calculation arguments; this table specifies which output calculations
% will be performed on the output from each simulation
output = setup_output(interest, value, add, time);

% Columns in the combined output table:
%
%   OUTPUTID = a unique number identifying each separate output calculation
%   to perform on the simulation output
%
%   INTEREST = the selector string for which molecules of interest to use
%   in the output calculation
%
%   CALC = the selector string for which output calculation to perform
%
%   ADD = the binary 1 or 0 for if the molecules of interest should be
%   added together or reported separately
%
%   TIME = the time point to return the concentration from for the "conct"
%   calculation type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Inputs for binding_main()

% binding_main() solves the system of ODEs for a specific set of initial
% conditions, parameters, and time periods and is called by binding_sim()
% to run a simulation for each set of inputs

% USAGE: [T, Y, NET, BAL] = BINDING_MAIN(TIMES, YIN, YTYPE, M, P, OPTS)

% m and p were already defined above for binding_sim

%% --------------- SINGLE SIMULATION ---------------

% binding_main() performs a simulation for a single set of initial
% concentrations for each species and parameters for each equation, so only
% the inputs for an individual simulation need to be selected

% Time periods are given as a vector where the first value is the start
% time of the simulation, the last value is the end time of the simulation,
% and any values in between are the end time of one period and the start
% time of the next period
times = [0 7200 14400];
% This indicates two periods, one from 0 to 7200 seconds and the other from
% 7200 to 14400 seconds

% The initial concentrations are given in a cell array where the each cell
% contains the concentrations for a different time period and each element
% is the concentration for a different species, in the same order as the
% species structure
% NaN or Inf indicate that the species starts at its ending concentration
% from the previous period (i.e., it does not change between periods)
y0 = {zeros(length(names(m)), 1), nan(length(names(m)), 1)};
y0{1}(m.IL6R) = 3.16e5;
y0{1}(m.IL8R) = 6.18e5;
y0{1}(m.BS1)  = 1000;
y0{2}(m.BS1)  = 0;
% Named y0 here to not overwrite the yin table needed for binding_sim()

% The period types are similarly given in a cell array or vector where each
% element is a string with the type of initial concentrations given:
% "initial" for the first time period where the starting concentrations are
% given, "washout" where concentrations are decreased to the values
% provided (typically 0), and "addition" where the specified concentration
% is added to the amount already present in the system
ytype = ["initial", "washout"];

% Specific parameter values have to be set for the individual simulation
p0 = p;
p0.kon_BS1_6R = kvalues(1);
p0.kon_BS1_8R = kvalues(2);
p0.kon_BS1_8R_6R = kvalues(3);
p0.kon_BS1_6R_8R = kvalues(4);
p0.koff_BS1_6R = kvalues(5);
p0.koff_BS1_8R = kvalues(5);
p0.koff_BS1_8R_6R = kvalues(5);
p0.koff_BS1_6R_8R = kvalues(5);
% Named p0 here to not overwrite the p structure needed for binding_sim()

% Additional options can be given for the ODE solver, including the
% absolute and relative tolerances and the time steps to use for each
% period, but those are not required inputs to the function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Inputs for binding_cost()

% binding_cost() calculates the difference between the experimental data
% and the normalized model output for a set of simulations and is used as
% the cost function for `lsqnonlin()` for parameter optimization

% USAGE: DIFF = BINDING_COST(KOPT, KOPT_NAMES, KIN, KIN_NAMES, ...
%   DATA, BASIS, CALC, AB, RECEP, SETUP, M, P)

% ab, recep, setup, m, and p were already defined above for binding_sim

%% --------------- EXPERIMENTAL DATA ---------------

% Multiple data files can be used for the same optimization if they
% represent the same experimental setup. In this case, I am using a single
% data file with data from two different replicates of the same experiment
file = "Flow-MFI_BS1-Norm";

% The `import_data` function reads in the data CSV file, converts the data
% to a tidy, long-form table, and adds a column for if the data point was
% used in the normalization of the data
data = import_data(file);

% Columns in the data table:
%
%   CONC = the initial concentation of the antibody in nM
%
%   SPECIES = the antibody used in the experiment (tocilizumab, 10H2, BS1,
%   or BS2)
%
%   CELL = the cell line used in the experiment (IL6R+, IL8R+, IL6R+/IL8R+,
%   or IL6R-/IL8R-)
%
%   REP = the replicate number; used for concatenating the data tables
%   together
%
%   VALUE = the normalized mean fluorescent intensity of the bound antibody
%
%   NORM = a binary 1 or 0 indicating whether the data in that row was used
%   in normalizing the antibody; used for selecting the model output points
%   that correspond to the same concentrations used to normalize the
%   experimental data

%% --------------- NORMALIZATION OPTIONS ---------------

% `norm_output` can calculate the normalization using different methods and
% with different bases

% Set which species to use as the basis for the normalization. "Ab"
% normalizes each antibody-cell line combination separately (e.g.,
% tocilizumab in IL6R+ cells is normalized against itself). "BS1"
% normalizes each antibody to the BS1 data for that cell line (e.g.,
% tocilizumab in IL6R+ cells is normalized against BS1 in IL6R+ cells). The
% value here should match the normalization done in the imported data set
basis = "BS1";

% Set how to perform the normalization. "data" will average the binding
% from the same concentrations in the model output that were averaged
% together to normalize the experimental data (based on the "Norm" column
% in the data table). "max" will select the maximum binding from the range
% of concentrations that match the experimental data concentrations
calc = "data";

%% --------------- PARAMETER VALUES ---------------

% The model parameter values are split into `kopt` and `kin` to distinguish
% between parameters that are currently being optimized (i.e., the values
% passed in from `lsqnonlin`) and the parameters that are required for the
% system but do not change during the optimization

% Although `binding_cost` is only returning the difference between the
% model output and the experimental data and is not actually doing the
% optimization, the parameters still need to be separated because
% `binding_cost` is the function given to `lsqnonlin`. `lsqnonlin` requires
% that the cost function is set up to take the values being optimized as a
% single argument

% Set the names/regular expression patterns and the values for the
% parameters being optimized; same idea as above in the `binding_sim`
% section where the pattern will be used to match the right fields in `p`
% and assign the corresponding parameter values
kopt_names = ["kon_Ab_6R$", "kon_Ab_8R$", "kon_Ab.*_[68]R_6R$"];
kopt = [1e-6 1e-6 1e-7];

% Set the names/regular expression patterns and the values for the
% parameters that are static and are not changed during the optimization
kin_names = "koff_Ab";
kin = 1e-4;

% The rate constant for binding to IL8R after already being bound to a
% receptor is not set here because the final rate constant in this system
% is constrained by the fact that the binding reactions occur in a cycle.
% Reactions in a cycle have a free energy change of 0, and this constrains
% the final parameter to be in a specific ratio with the other parameters.
% `binding_cost` will automatically calculate the missing parameter based
% on this constraint and the values of the provided parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Inputs for binding_opt()

% binding_opt() calls lsqnonlin() to optimize the binding rate constant
% values based on the provided experimental data. lsqnonlin() calls
% binding_cost() as the cost function for the normalization

% USAGE: [OPTIMAL, KNAMES, FINAL_COST, FLAG] = BINDING_OPT(...
%   KOPT, KOPT_NAMES, LB, UB, KIN, KIN_NAMES, DATA, BASIS, ...
%   CALC, AB, RECEP, SETUP, M, P)

% kopt, kopt_names, kin, kin_names, data, basis, calc, ab, recep, setup, m,
% and p were already defined above for binding_cost and binding_sim

%% --------------- OPTIMIZATION OPTIONS ---------------

% The lower and upper bounds for the simulation are defined for `lsqnonlin`

% There must be a single lower bound and upper bound for each value of kopt
lb = [1e-12 1e-12 1e-14];
ub = [1e-2 1e-2 1e-2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Inputs for norm_output()

% norm_output() normalizes the model output based on the experimental data
% and the options provided and returns the normalized output and the
% difference between the normalized output and the experimental data. It is
% called within binding_cost() to perform the normalization for the
% parameter optimization

% USAGE: [NORM, DIFF] = NORM_OUTPUT(OUT, INPUT, DATA, BASIS, CALC, ...
%   TIME, OPTION)

% input, data, basis, and calc were already defined above for binding_cost
% and binding_sim

%% --------------- RUN SIMULATION ---------------

% `norm_output` normalizes the model output after the simulation has been
% performed, so it requires the model output as an input argument

% `binding_sim` runs the model simulations based on the input initial
% conditions and parameter values
out = binding_sim(id, yin, params, output, m, p);

% The normalization requires the concentrations at each time point for the
% antibody-receptor complexes. The other output calculations are not used
% in the normalization, so they can be filtered out

% Determine the output ID that corresponds to the concentrations at each
% time point for the antibody-receptor complexes
id_out = output{output.Interest == "receptors" & ...
    output.Calc == "conc" & output.Opts.Add == 1, ...
    "OutputID"};

% Select the rows of the model output from the correct output ID
out = out(out.OutputID == id_out, :);

%% --------------- NORMALIZATION OPTIONS ---------------

% Specify which time point to use for the normalization; this should be the
% same time point that was used for the experimental data
time = 8100;
% In this experiment, the antibodies were incubated for 2 hours on the
% cells, and then the unbound antibody was washed off and the secondary
% antibody was incubated for 15 minutes. The measurements were taken at 2
% hours and 15 minutes, so the model output should also be normalized at 2
% hours and 15 minutes.

% Specify if the difference and the normalized should be returned as
% matrices or tables. The tables are easier to understand and work with -
% the matrix option is provided primarily for `binding_cost` since it only
% needs the actual values for the difference
option = "table";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save MAT File

% Save the sample inputs as a MAT file for easy access
save example/sample.mat id yin params output m p times y0 ytype p0 ...
    kopt kopt_names kin kin_names data basis calc ab recep setup ...
    lb ub out time option

end
