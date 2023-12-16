function output = calc_complexes(T, Y, net, bal, m, p, ab, recep, ...
    bound, singles, calc, opts)

% CALC_COMPLEXES	Calculate simulation output for total bound receptor.
%
%	Combination of SELECT_MOLECULES and CALC_OUTPUT specifically for the
%	total bound receptor. Total bound receptor needs to be handled
%	separately because it is a special case of [binary] + 2 * [ternary].
%
%	USAGE:
%		OUTPUT = CALC_COMPLEXES(T, Y, NET, BAL, M, P, AB, RECEP, ...
%           BOUND, SINGLES, CALC, OPTS)
%
%	INPUT:
%       T = vector of time values corresponding to the concentrations in Y,
%       given in seconds
%
%       Y = matrix of concentration values output from the ODEs; each
%       column corresponds to an equation in order, each row is a time
%       point
%
%       NET = matrix of net input concentrations output from the ODEs; each
%       column corresponds to an input molecule, each row is a time point
%
%       BAL = matrix of mole balance values output from the ODEs; each
%       column corresponds to an input molecule, each row is a time point
%
%       M = structure containing the molecules included in the model and
%       their assigned numbers, input as a structure, correct naming for
%       species is antibody name or antibody receptor complex name (with
%       components separated by underscores)
%
%       P = structure containing parameters for the model, input as a
%       structure; correct naming for parameters is "kon_" or "koff_",
%       followed by the antibody name or antibody receptor complex name
%       (with components separated by underscores), ending with the
%       receptor being bound to; e.g., kon_BS1_6R_8R
%
%       BOUND = antibody-receptor combinations corresponding to
%       molecules in the m structure; e.g., only includes BS1_6R_8R and not
%       BS1_8R_6R because they are the same molecule
%
%       SINGLES = all of the individual antibodies and receptors that
%       contribute to the bound complexes in BOUND
%
%       CALC = a string specifying which output calculation to use; options
%       are ["conc", "conct", "end", "peak", "auc", "net", "bal", "dydt"]
%
%       OPTS = a structure or table containing the options for the
%       calculations; has two fields: "Add", a logical for whether to sum
%       the output or output each individual molecule, and "Time", a vector
%       or cell array of time points to use for "conct"
%
%	OUTPUT:
%		OUTPUT = a long form table with three columns, one for the time
%		points for the output, one with the species names used in the
%		output calculation, and one with the output values
%
%	NOTES:
%		This is a special case of SELECT_MOLECULES and CALC_OUTPUT
%		specifically for the total bound receptor. For Tocilizumab, 10H2,
%		and BS1, the total bound receptor is calculated by [binary] + 2 *
%		[ternary]. Because of how CALC_OUTPUT handles the searching of
%		molecules and the adding of outputs, it was more straightforward to
%		separate these outputs into a separate wrapper function.
%
%	See also SELECT_MOLECULES, CALC_OUTPUT, BINDING_SIM.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    T (:,1) double
    Y double
    net double = 0
    bal double = 0
    m struct = struct()
    p struct = struct()
    ab (:,1) string = strings(0)
    recep {mustBeA(recep, ["string", "cell"])} = strings(0)
    bound (:,1) string = strings(0)
    singles (:,1) string = strings(0)
    calc (1,1) string {mustBeMember(calc, ...
        ["conc", "conct", "end", "peak", "auc", ...
        "net", "bal", "dydt"])} = "conc"
    opts {mustBeA(opts, ["table", "struct"])} = struct()
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select Molecules of Interest

% Determine binary and ternary complexes present in system
binary = select_molecules(ab, recep, "binary", bound, singles);
ternary = select_molecules(ab, recep, "ternary", bound, singles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Model Output

% The "complexes" calculation requires the values to be added together
opts.Add = 1;

% Get the required direct outputs from the calculation function based on
% the requested calculation
if calc ~= "peak"

    % Get values for binary and ternary complexes separately
    binary_out = calc_output(T, Y, net, bal, m, p, singles, binary, ...
        calc, opts);
    ternary_out = calc_output(T, Y, net, bal, m, p, singles, ternary, ...
        calc, opts);

    % Add the values together and re-create the output table
    % If there is no ternary species, remove the trailing "+2*"
    values = binary_out.Value + ternary_out.Value * 2;
    names = paste([binary_out.Species, ternary_out.Species], "+2*");
    names = regexprep(names, "\+2\*$", "");
    output = table(binary_out.Time, names, values, ...
        'VariableNames', ["Time", "Species", "Value"]);

% Peak concentration is a special case because the peak values cannot just
% be added together
elseif calc == "peak"

    % Get concentrations for binary and ternary complexes separately
    % Get values for binary and ternary complexes separately
    binary_out = calc_output(T, Y, net, bal, m, p, singles, binary, ...
        "conc", opts);
    ternary_out = calc_output(T, Y, net, bal, m, p, singles, ternary, ...
        "conc", opts);

    % Add the values together and determine the time of the peak value
    values = binary_out.Value + ternary_out.Value * 2;
    [peak, idx] = max(values);
    peakT = binary_out.Time(idx);

    % Re-create the output table
    % If there is no ternary species, remove the trailing "+2*"
    names = paste([binary_out.Species(1), ternary_out.Species(1)], "+2*");
    names = regexprep(names, "\+2\*$", "");
    output = table(peakT, names, peak, ...
        'VariableNames', ["Time", "Species", "Value"]);
end
