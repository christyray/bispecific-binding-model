function [kout, koutname] = thermo_ratio(kvalues, knames)

% THERMO_RATIO	Calculate a parameter value based on the thermo cycle.
%
%	Calculates the value of one of the binding rate constants for a
%	bivalent antibody based on the other rate constants and the
%	thermodynamic cycle relationship.
%
%	USAGE:
%		[KOUT, KOUTNAME] = THERMO_RATIO(KVALUES, KNAMES)
%
%	INPUT:
%       KVALUES = a matrix of on and off rate constants for the binding
%       model; rows should correspond to different values to calculate, and
%       columns should correspond to the species in KNAMES
%
%       KNAMES = a string vector of regular expression patterns
%       correponding to the parameter values given in KVALUES
%
%
%	OUTPUT:
%		KOUT = the rate constant value that was not provided in the input
%		values
%
%       KOUTNAME = the regular expression pattern corresponding to the name
%       of the KOUT parameter
%
%	NOTES:
%		Uses the thermodynamic cycle relationship between the k_on and
%       k_off values.
%
%       KOUT and KOUTNAME can be concatenated to the existing KVALUES and
%       KNAMES to complete the parameter definitions for the binding model
%       simulation.
%
%	See also BINDING_COST, BINDING_OPT, BINDING_DRIVER.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    kvalues (:,:) double
    knames (1,:) string
end

% Confirm that the number of values given is equal to the number of names
if ~any(size(kvalues) == length(knames))
    eid = 'Size:notEqual';
    msg = "The number of values and names given must be equal.";
    throwAsCaller(MException(eid, msg))
end

% Reshape kvalues so the columns correspond to the knames
if size(kvalues, 2) ~= length(knames)
    kvalues = kvalues';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign Values

% Find value indices using the provided name patterns
param_names = ["kon_Ab_6R"; "kon_Ab_8R"; ...
    "kon_Ab_8R_6R"; "kon_Ab_6R_8R"; ...
    "koff_Ab_6R"; "koff_Ab_8R"; ...
    "koff_Ab_8R_6R"; "koff_Ab_6R_8R"];
params = zeros(size(kvalues, 1), length(param_names));

% For each input parameter name
for i = 1:length(knames)
    idx = regexpl(param_names, knames(i));
    params(:,idx) = repmat(kvalues(:,i), 1, sum(idx));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cycle Calculation

% Performs the calculation based on which input parameter is absent
missing = find(params(1,:) == 0);

% Throw an error if too many parameters are missing
if length(missing) > 1
    eid = 'Arguments:notEnoughInputs';
    msg = "Not enough input parameters provided.";
    throwAsCaller(MException(eid, msg))

% Throw an error if all parameters are provided, but they are do not meet
% the thermodynamic cycle constraint
elseif isempty(missing)

    % Calculate the equilibrium constants
    K6R         = params(:,5) ./ params(:,1);
    K8R         = params(:,6) ./ params(:,2);
    K6Rprime    = params(:,7) ./ params(:,3);
    K8Rprime    = params(:,8) ./ params(:,4);

    % Calculate the ratio of the equilibrium constants
    ratio       = K6R ./ K6Rprime .* K8Rprime ./ K8R;

    % Throw an error if the ratio is not approximately one
    if ~ismembertol(ratio, 1, 1e-8)
        eid = 'Calculation:constantsIncorrectRatio';
        msg = "The provided constants do not satisfy the constraint.";
        throwAsCaller(MException(eid, msg))

    % Set the output values to empty because no other calculation is needed
    % and exit the function
    else
        kout = [];
        koutname = strings(0);
        return
    end
end

% Determine calculation based on missing parameter
switch missing
    case 1  % kon_6R
        % Calculate the possible equilibrium constants
        K8R         = params(:,6) ./ params(:,2);
        K6Rprime    = params(:,7) ./ params(:,3);
        K8Rprime    = params(:,8) ./ params(:,4);

        % Calculate missing equilibrium constant
        K6R         = K8R ./ K8Rprime .* K6Rprime;

        % Calculate missing value
        kout        = params(:,5) ./ K6R;
        koutname    = "kon_Ab_6R$";

    case 2  % kon_8R
        % Calculate the possible equilibrium constants
        K6R         = params(:,5) ./ params(:,1);
        K6Rprime    = params(:,7) ./ params(:,3);
        K8Rprime    = params(:,8) ./ params(:,4);

        % Calculate missing equilibrium constant
        K8R         = K6R ./ K6Rprime .* K8Rprime;

        % Calculate missing value
        kout        = params(:,6) ./ K8R;
        koutname    = "kon_Ab_8R$";

    case 3  % kon_6R_prime
        % Calculate the possible equilibrium constants
        K6R         = params(:,5) ./ params(:,1);
        K8R         = params(:,6) ./ params(:,2);
        K8Rprime    = params(:,8) ./ params(:,4);

        % Calculate missing equilibrium constant
        K6Rprime    = K6R ./ (K8R ./ K8Rprime);

        % Calculate missing value
        kout        = params(:,7) ./ K6Rprime;
        koutname    = "kon_Ab.*_[68]R_6R$";

    case 4  % kon_8R_prime
        % Calculate the possible equilibrium constants
        K6R         = params(:,5) ./ params(:,1);
        K8R         = params(:,6) ./ params(:,2);
        K6Rprime    = params(:,7) ./ params(:,3);

        % Calculate missing equilibrium constant
        K8Rprime    = K8R ./ (K6R ./ K6Rprime);

        % Calculate missing value
        kout        = params(:,8) ./ K8Rprime;
        koutname    = "kon_Ab.*_[68]R_8R$";

    case 5  % koff_6R
        % Calculate the possible equilibrium constants
        K8R         = params(:,6) ./ params(:,2);
        K6Rprime    = params(:,7) ./ params(:,3);
        K8Rprime    = params(:,8) ./ params(:,4);

        % Calculate missing equilibrium constant
        K6R         = K8R ./ K8Rprime .* K6Rprime;

        % Calculate missing value
        kout        = params(:,1) .* K6R;
        koutname    = "koff_Ab_6R$";

    case 6  % koff_8R
        % Calculate the possible equilibrium constants
        K6R         = params(:,5) ./ params(:,1);
        K6Rprime    = params(:,7) ./ params(:,3);
        K8Rprime    = params(:,8) ./ params(:,4);

        % Calculate missing equilibrium constant
        K8R         = K6R ./ K6Rprime .* K8Rprime;

        % Calculate missing value
        kout        = params(:,2) .* K8R;
        koutname    = "koff_Ab_8R$";

    case 7  % koff_6R_prime
        % Calculate the possible equilibrium constants
        K6R         = params(:,5) ./ params(:,1);
        K8R         = params(:,6) ./ params(:,2);
        K8Rprime    = params(:,8) ./ params(:,4);

        % Calculate missing equilibrium constant
        K6Rprime    = K6R ./ (K8R ./ K8Rprime);

        % Calculate missing value
        kout        = params(:,3) .* K6Rprime;
        koutname    = "koff_Ab.*_[68]R_6R$";

    case 8  % koff_8R_prime
        % Calculate the possible equilibrium constants
        K6R         = params(:,5) ./ params(:,1);
        K8R         = params(:,6) ./ params(:,2);
        K6Rprime    = params(:,7) ./ params(:,3);

        % Calculate missing equilibrium constant
        K8Rprime    = K8R ./ (K6R ./ K6Rprime);

        % Calculate missing value
        kout        = params(:,4) .* K8Rprime;
        koutname    = "koff_Ab.*_[68]R_8R$";
end
