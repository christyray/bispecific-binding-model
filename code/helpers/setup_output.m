function output = setup_output(interest, calc, add, time)

% SETUP_OUTPUT	Generate a table of outputs to collect from BINDING_SIM.
%
%   Combines the interest, calculation, add, and time objects into a table
%   formatted correctly for use in the BINDING_SIM function.
%
%   USAGE:
%       OUTPUT = SETUP_OUTPUT(INTEREST, CALC, ADD, TIME)
%
%   INPUT:
%       INTEREST = a string vector with the molecules of interest to pass
%       to SELECT_MOLECULES; options are ["inputs", "antibodies",
%       "receptors", "binary", "ternary", "complexes", "present", "all"]
%
%       CALC = a string vector of the output calculations to use for
%       CALC_OUTPUT; options are ["conc", "conct", "end", "peak", "auc",
%       "net", "bal", "dydt"]
%
%       ADD = a logical vector 1 or 0 for if the output in CALC_OUTPUT
%       should be summed or if each individual molecule of interest should
%       be reported; default is 0
%
%       TIME = a vector or cell array of time points to calculate the
%       concentration at when the "conct" calculation type is used; default
%       is 0
%
%   OUTPUT:
%       OUTPUT = table with all combinations of antibody/receptor combos,
%       molecules of interest, calculation types, and options for
%       OUTPUT_CALC
%
%   NOTES:
%       This function generates the output table to be used as an input to
%       the BINDING_SIM function.
%
%   See also BINDING_SIM, CALC_OUTPUT, SELECT_MOLECULES, COMBINE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    interest (:,1) string
    calc (:,1) string
    add (:,1) double = 0
    time (:,1) {mustBeA(time, ["double", "cell"])} = 0
end

% Convert the molecules of interest and the calculation type to lowercase
interest = lower(interest);
calc = lower(calc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Output Table

% Combine all of the variables for the output into one table
output = combine(interest, calc, add, time);

% Rename the table columns for standardization
output = rename(output, ...
    new=["Interest", "Calc", "Add", "Time"]);

% Select only the rows with relevant calculations
% The "Time" option is only relevant for the "conct" calculation
[~,idx] = unique(output(:, ["Interest", "Calc", "Add"]), "rows", "stable");
rows = ismember((1:height(output))', idx) | output.Calc == "conct";
output = output(rows, :);

% "inputs", "present", and "all" should not be added up because the units
% are not consistent
condition = any(output.Interest == ["inputs", "present", "all"], 2) ...
    & output.Add == 1;
output{condition, "Add"} = 0;

% "complexes" must be added because it is specifically the total bound
% receptor (binary + 2 * ternary)
condition = output.Interest == "complexes" & output.Add == 0;
output{condition, "Add"} = 1;

% Remove any duplicated rows after restricting the calculation options
output = unique(output, "stable");

% Merge the options column for CALC_OUTPUT
output = mergevars(output, ["Add", "Time"], ...
    'MergeAsTable', true, 'NewVariableName', "Opts");

% Add an ID column
output = addvars(output, (1:height(output))', 'Before', 1, ...
    'NewVariableNames', "OutputID");

end
