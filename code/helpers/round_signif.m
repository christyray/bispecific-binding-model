function rounded = round_signif(n, digits, option)

% ROUND_SIGNIF	Round to a number of significant digits with options.
%
%	Round input to a specified number of significant figures using any of
%	MATLAB's rounding functions: "round", "ceil", "floor", or "fix."
%
%	USAGE:
%		ROUNDED = ROUND_SIGNIF(N, DIGITS, OPTION)
%
%	INPUT:
%		N = number or vector of numbers to round
%
%       DIGITS = number of significant digits to round to; default is 1
%
%       OPTION = rounding function to use; options are "round", "floor",
%       "ceil", or "fix"; defaults to "round"
%
%	OUTPUT:
%		ROUNDED = rounded input number(s)
%
%	NOTES:
%		Includes a catch for floating point error, but may still give
%		erroneous output at large numbers of digits due to floating point
%		error.
%
%	See also ROUND, CEIL, FLOOR, FIX.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    n {mustBeNumeric, mustBeVector}
    digits (1,1) {mustBeA(digits, "double")} = 1
    option string {mustBeMember(option, ...
        ["round", "floor", "ceil", "fix"])} = "round"
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Round Inputs

% Order of magnitude of the input values
power = ceil(log10(abs(n) .* (1 + eps(n))));

% Resolution to round to
resolution = 10.^(power - digits);

% Pre-rounded divided number to avoid floating point error
rounded = round(n./resolution, digits + 2, 'significant');

% Use the specified rounding function
switch option
    case "round"
        rounded = round(rounded) .* resolution;
    case "floor"
        rounded = floor(rounded) .* resolution;
    case "ceil"
        rounded = ceil(rounded) .* resolution;
    case "fix"
        rounded = fix(rounded) .* resolution;
end

% Locate any values that cannot be handled in the inputs and replace the
% respective outputs

% Infinity and negative infinity
rounded(isinf(n)) = n(isinf(n));

% Zeros
rounded(n == 0) = 0;
