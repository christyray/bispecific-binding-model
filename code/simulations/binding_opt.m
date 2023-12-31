function [optimal, knames, final_cost, flag, iter] = binding_opt(...
    kopt, kopt_names, lb, ub, kin, kin_names, data, basis, calc, ...
    ab, recep, setup, m, p, options)

% BINDING_OPT	Optimization function for binding model parameters.
%
%	Optimization setup for the determination of the antibody-receptor
%   binding rates from the flow cytometry equilibrium binding data; wrapped
%   as a function so it can be used with different initial guesses.
%
%	USAGE:
%		[OPTIMAL, KNAMES, FINAL_COST, FLAG, ITER] = BINDING_OPT(...
%           KOPT, KOPT_NAMES, LB, UB, KIN, KIN_NAMES, DATA, BASIS, ...
%           CALC, AB, RECEP, SETUP, M, P, OPTIONS)
%
%	INPUT:
%       KOPT = initial guesses for the parameter(s) to be optimized
%
%       KOPT_NAMES = string vector with the regular expression patterns or
%       names for each parameter to be optimized
%
%       LB = lower bounds for the optimization; must be same length as KOPT
%
%       UB = upper bounds for the optimization; must be same length as KOPT
%
%       KIN = vector with values for the parameters not being optimized
%
%       KIN_NAMES = string vector with the regular expression patterns or
%       names for each parameter not being optimized
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
%       AB = a string vector with the antibodies to use in the complexes
%
%       RECEP = a cell array of string vectors with the receptor(s) to use
%       in the complexes; the rows of the array should correspond to the
%       antibodies in AB
%
%		SETUP = experimental design used, must be a string
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
%		OPTIONS = named arguments for what type of display to use from the
%		optimization iterations and whether to write the display to a diary
%		file
%
%	OUTPUT:
%		OPTIMAL = optimized parameter values generated by LSQNONLIN;
%		includes the static parameters and the parameter calculated from
%		the thermodynamic cycle constraint (if applicable)
%
%       KNAMES = regular expression patterns or names corresponding to the
%       values in OPTIMAL
%
%       FINAL_COST = sum of squares of the difference between the
%       experimental data and the model output for the optimized parameter
%       values
%
%       FLAG = flag output by LSQNONLIN describing why the optimization
%       stopped
%
%       ITER = number of iterations performed by LSQNONLIN in the
%       optimization
%
%	NOTES:
%		When multiple optimizations are performed in a single loop, KNAMES
%		will be the same for all optimizations.
%
%	See also BINDING_COST, BINDING_SIM, LSQNONLIN.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    kopt (1,:) double
    kopt_names (1,:) {mustBeA(kopt_names, "string")}
    lb (1,:) double
    ub (1,:) double
    kin (1,:) double
    kin_names (1,:) {mustBeA(kin_names, "string")}
    data {mustBeA(data, ["table", "cell"])}
    basis string {mustBeMember(basis, ["Ab", "BS1"])}
    calc string {mustBeMember(calc, ["data", "max"])}
    ab string
    recep {mustBeA(recep, ["string", "cell"])}
    setup string
    m struct
    p struct
    options.Display string {mustBeMember(options.Display, ...
        ["iter", "iter-detailed", "final", "final-detailed", ...
        "none"])} = "none"
    options.Diary string = strings(0)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optimization

% Start timer for optimization code
tic_opt_iter = tic;

% Initialize variables shared with the output function - used for
% collecting the values after each iteration
valID = 1;

% Variables with shared scope:
% valID = ID of current output text file, set by MATLAB when the file is
% opened and used to ensure that the data for each optimization writes to
% the correct place
% kopt_names = used as the column names for the table of current parameter
% values
% options = provides the name to use for the diary file

% The cost function file requires multiple inputs to allow it to run
% inside the loop. Because lsqnonlin requires the cost function to only
% take the parameters to be optimized as input, this anonymous
% function gives the cost function the other parameters
cost = @(kopt) binding_cost(kopt, kopt_names, kin, kin_names, ...
    data, basis, calc, ab, recep, setup, m, p);

% Optimization options
lsqoptions = optimoptions('lsqnonlin', 'Display', options.Display, ...
    'OutputFcn',@outfun);
% Setting the 'OutputFcn' makes lsqnonlin call my `outfun()` after each
% iteration, which will print the current parameter values in a separate
% log file for each optimization (especially useful when running on the
% cluster)

% If a diary file was provided, save the iteration display to the diary
if ~isempty(options.Diary)
    diaryfile = options.Diary + "-iter.log";
    diary(diaryfile);
end

% Run lsqnonlin for the anonymous cost function and find the cost for the
% optimized values
[optimal, final_cost, ~, flag, optimout] = lsqnonlin(cost, kopt, ...
    lb, ub, lsqoptions);

% Turn off the diary if it was initiated
diary off

% Combine the optimized parameters and the static parameters
optimal = [optimal, kin];
knames = [kopt_names, kin_names];

% Calculate the input parameter that was not provided using the
% thermodynamic cycle ratio
[kvalue, kname] = thermo_ratio(optimal, knames);
optimal = [optimal, kvalue];
knames = [knames, kname];

% Select just the number of iterations from the optimout structure
iter = optimout.iterations;

% If a diary file was provided, save the elapsed time to the diary
if ~isempty(options.Diary)
    elapsed = code_time(tic_opt_iter);
    diaryID = fopen(diaryfile, 'a');
    fprintf(diaryID, "Single optimization completed in %s.\n\n", elapsed);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions

    % Collect the parameter values used for each iteration
    % Nested because it has a shared scope with the binding_opt function
    function stop = outfun(x, optimValues, state)

        % Stop allows the output function to add additional stopping
        % criteria for lsqnonlin, but that isn't needed here
        stop = false;

        % State allows for different behavior depending on what the
        % optimization solver is currently doing
        switch state
            % After the solver is initialized
            case 'init'
                % If a diary file was provided, save to the diary
                % Default value of valID = 1, set at the top of the file
                if ~isempty(options.Diary)
                    valfile = options.Diary + "-xval.log";
                    valID = fopen(valfile, 'w');
                end

                % Print table header
                header = pad(replace_names(kopt_names), 15, 'left', " ");
                header = [pad("Iteration", 10, 'left', " "), header];
                header = strjoin(header, "");
                fprintf(valID, '%s\n', header);

            % After each iteration is complete
            case 'iter'
                % Print table row
                row = [optimValues.iteration, x];
                formatSpec = "%10d" + ...
                    strjoin(repmat("%15.4e", length(x), 1), "") + "\n";
                fprintf(valID, formatSpec, row);
        end
    end
end
