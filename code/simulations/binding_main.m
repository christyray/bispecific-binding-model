function [T, Y, net, bal] = binding_main(times, yin, ytype, m, p, opts)

% BINDING_MAIN	Solves the system of ODEs in binding_eqns.
%
%	Solves the system of ODEs given in the binding_eqns file based on
%   the parameters given. Inputs are the parameters and the length of
%   time; outputs are the time and concentration matrices from the ODE
%   solver.
%
%	USAGE:
%		[T, Y, NET, BAL] = BINDING_MAIN(TIMES, YIN, YTYPE, M, P, OPTS)
%
%	INPUT:
%       TIMES = vector of times to use for the beginning and end of the
%       system, given in seconds; if multiple periods, give the beginning,
%       end, and every change point
%
%       YIN = initial values for the model, input as vectors in a cell
%       array where each vector corresponds to the initial conditions for a
%       specific time period; the order of values in the vector should
%       correspond to the order of species in the molecule structure
%
%       YTYPE = classification for each set of initial values, input as a
%       string vector; options are "initial", "addition", and "washout";
%       tells the model how to add or subtract from the previous conditions
%       to get the values in the system at the start of the time period
%
%       M = structure containing the molecules included in the model and
%       their assigned numbers, input as a structure; correct naming for
%       species is antibody name or antibody receptor complex name (with
%       components separated by underscores)
%
%       P = structure containing parameters for the model, input as a
%       structure; correct naming for parameters is "kon_" or "koff_",
%       followed by the antibody name or antibody receptor complex name
%       (with components separated by underscores), ending with the
%       receptor being bound to; e.g., kon_BS1_6R_8R
%
%       TSTEP = optional and named, contained within OPTS; value for the
%       time step to be used in the ODE solver output, given in seconds;
%       default is 30
%
%       ABSTOL = optional and named, contained within OPTS; value for the
%       absolute tolerance for the ODE solver options
%
%       RELTOL = optional and named, contained within OPTS; value for the
%       relative tolerance for the ODE solver options
%
%       ABS_BAL = optional and named, contained within OPTS; value for the
%       absolute tolerance to be used for the mole balance; default is
%       1e-12
%
%       REL_BAL = optional and named, contained within OPTS; value for the
%       relative tolerance to be used for the mole balance; default is 2e-5
%
%	OUTPUT:
%       T = vector of time values corresponding to the concentrations in Y,
%       given in seconds
%
%       Y = matrix of concentration values output from the ODEs; each
%       column corresponds to an equation in order, each row is a time
%       point
%
%       NET = matrix with the amount added or removed from the system for
%       each species at each time point; can be output if debugging
%
%       BAL = matrix with the mole balance for each species at each
%		time point; can be output if debugging
%
%	NOTES:
%		BINDING_MAIN sets up the initial conditions for the time period
%		based on the concentrations in the system at the end of the
%		previous time period and the type of the next time period.
%
%       For "initial", it uses the initial values as provided. For
%       "addition", the values provided will be added to the ending
%       concentrations, and any species that are not changing should be
%       given as NaN. For "washout", the the initial conditions will be set
%       to the values provided, and any species that are not changing
%       should be given as NaN.
%
%       BAL does not need to be output by default but is provided as an
%       output for debugging purposes.
%
%	See also BINDING_EQNS, MOLE_BALANCE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    times double
    yin {mustBeA(yin, ["double", "cell"])}
    ytype (1,:) string
    m struct
    p struct
    opts.tstep double = 30
    opts.AbsTol double = 1e-8
    opts.RelTol double = 1e-6
    opts.abs_bal double = 1e-12
    opts.rel_bal double = 2e-5
end

% If only one time step is given, repeat it for each time period
if length(opts.tstep) == 1
    opts.tstep = repmat(opts.tstep, 1, length(ytype));
end

% -------------------- ERRORS --------------------
% Number of time periods, time steps, initial condition types, and initial
% condition vectors must be equal
if length(ytype) == 1 && class(yin) == "cell" && any(size(yin) > 1)
    eid = 'Size:notEqual';
    msg = 'Too many initial conditions provided; only one type given.';
    throwAsCaller(MException(eid,msg))
elseif length(ytype) == 1 && ~any(size(yin) == 1)
    eid = 'Size:notEqual';
    msg = 'Too many initial conditions provided; only one type given.';
    throwAsCaller(MException(eid,msg))
elseif length(ytype) > 1 && ~any(size(yin) == length(ytype))
    eid = 'Size:notEqual';
    msg = 'Initial condition and type length are not equal.';
    throwAsCaller(MException(eid,msg))
elseif length(ytype) ~= (length(times) - 1)
    eid = 'Size:notEqual';
    msg = 'Number of time periods is not equal to number of period types.';
    throwAsCaller(MException(eid,msg))
elseif length(opts.tstep) ~= length(ytype)
    eid = 'Size:notEqual';
    msg = 'Number of time steps is not equal to number of time periods.';
    throwAsCaller(MException(eid,msg))
end

% -------------------- CONVERSION --------------------
% Convert initial conditions to cell array for correct indexing
if length(ytype) == 1 && class(yin) == "double"
    yin = {yin};
elseif class(yin) == "double"

    % Reshape so the rows correspond to different time periods
    if size(yin, 1) ~= length(ytype)
        yin = yin';
    end

    % Convert to a cell array with each time period in a different cell
    rowdim = ones(1, size(yin, 1));
    yin = mat2cell(yin, rowdim);
end

% Replace any Inf values with NaN - allows other functions to work with Inf
% to make unique() work (Inf == Inf but NaN ~= NaN)
for i = 1:length(yin)
    yin{i}(isinf(yin{i})) = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ODE Solver

% Convert ytype to lowercase for standardization
ytype = lower(ytype);

% Set the ODE solver options
neqns = length(fieldnames(m));
options = odeset('AbsTol', opts.AbsTol, 'RelTol', opts.RelTol, ...
    'NonNegative', 1:neqns, 'InitialStep', 1e-2);

% Initialize result vectors
T = [];
Y = [];

% Initialize array for net amount added to system for mole balance
% Has a column for each species, same as Y
net = zeros(0, neqns);

% Initialize any conditions not specified to zero
y0 = zeros(1, neqns);

% Solve for each time period
for j = 1:(length(times)-1)

    % Select time span for current period
    tspan = [times(j) times(j+1)];

    % Adjusts the tspan vector for points based on the input time step
    tspan = tspan(1):opts.tstep(j):tspan(2);

    % Select the initial conditions for the current period
    switch ytype(j)
        case "initial"      % For the first period

            % Set the given initial conditions
            y0(1:length(yin{j}))    = yin{j};

        case "addition"     % When more is added to the system

            % For each species that changes, add the new concentration to
            % the previous
            % When giving conditions, use NaN for any species not changed
            idx = ~isnan(yin{j});
            y0(idx) = y0(idx) + yin{j}(idx);

        case "washout"      % When some is removed from the system

            % For each species that changes, set the concentration to the
            % given value
            % When giving conditions, use NaN for any species not changed
            idx = ~isnan(yin{j});
            y0(idx) = yin{j}(idx);
    end

    % Run the solver
    [T1,Y1] = ode15s(@(t,y) binding_eqns(t,y,m,p),tspan,y0,options);

    % Determine the total amount added or removed in the time period
    if j == 1
        % For the first time period, the amount added is the initial
        % condition
        value = y0;
    else
        value = y0 - Y(end,:) + net(end,:);
        % In system at beginning of this time period - in system at end
        % of previous time period + total amount added before now
        % Total change = In system now - what was in system before
        % washout/addition
        % Total for balance = Total change + previous total for balance
    end

    % Create an array with the total change for each species
    change = repmat(value, length(T1), 1);

    % Concatenate the change during this time period to the full net array
    net = cat(1, net, change);

    % Concatenate the results to the full results vector
    T = cat(1,T,T1);
    Y = cat(1,Y,Y1);

    % Set the next y0 to be equal to the end conditions
    y0 = Y(end,:);

end

% Calculate mole balance, will display an error message and table if the
% moles are not balanced
bal = mole_balance(Y,net,m,p,opts.abs_bal,opts.rel_bal);
end
