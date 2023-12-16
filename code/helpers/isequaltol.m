function iseq = isequaltol(A, B, tol, options)

% ISEQUALTOL	Determine array equality within a tolerance.
%
%	Element by element comparison of two arrays or tables A and B within a
%	specified tolerance. Returns a logical 1 (TRUE) where the values are
%	equal within the tolerance and a logical 0 (FALSE) where the values are
%	not equal within the tolerance.
%
%	USAGE:
%		ISEQ = ISEQUALTOL(A, B, TOL, OPTIONS)
%
%	INPUT:
%		A = first array or table for comparison
%
%       B = second array or table for comparison; each dimension must
%       either be the same as A or 1
%
%       TOL = number for the tolerance for the comparison; default is 1e-12
%
%       TYPE = which type of tolerance to use, "relative" or "absolute";
%       default is "relative"
%
%	OUTPUT:
%		ISEQ = logical array with the same dimensions as the inputs with
%		TRUE where the values are equal within the tolerance and FALSE
%		where the values are not equal within the tolerance.
%
%	NOTES:
%       Designed for comparing two numerical arrays that may fail an
%       equality test due to floating-point error. Counterpart to
%       ISMEMBERTOL to specifically compare matching elements instead of
%       searching for any match.
%
%	See also ISEQUAL, ISMEMBER, ISMEMBERTOL.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    A {mustBeA(A, ["double", "table"])}
    B {mustBeA(B, ["double", "table"])}
    tol (1,1) double = 1e-12;
    options.Type string {mustBeMember(options.Type, ...
        ["relative", "absolute", "rel", "abs"])} = "relative"
end

% Throw an error if the input arrays are not compatible sizes
if ~all(size(A) == size(B) | size(A) == 1 | size(B) == 1)
    eid = 'Size:notEqual';
    msg = 'Provided arrays are not compatible sizes.';
    error(eid,msg)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Standardize Inputs

% If both inputs were tables and both tables use the same variable names
% but the names are different places, sort the table columns to match the
% first input array
if (class(A) == "table" && class(B) == "table" ...
        && ~all(names(A) == names(B)) ...
        && all(ismember(names(A), names(B))))
    B = movevars(B, names(A), 'Before', 1);
end

% Convert inputs from tables to arrays if needed
if class(A) == "table"; A = table2array(A); end
if class(B) == "table"; B = table2array(B); end

% If A or B are vectors, convert them to arrays of matching size
sz = max(size(A), size(B));
A = A .* ones(sz);
B = B .* ones(sz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare Equality with Tolerance

% Determine equality depending on the type of tolerance being used
switch options.Type
    case {"relative", "rel"}
        % max(A, B) returns the larger value from each position in A and B
        % Using the larger value between the two matrices to scale the
        % relative tolerance for the error calculation
        reltol = @(x,y) abs(x - y) <= max(abs(x), abs(y)) * tol;

        % Separate calculations for real and imaginary components
        iseq = reltol(real(A), real(B)) & reltol(imag(A), imag(B));

    case {"absolute", "abs"}
        % Absolute tolerance does not require scaling
        abstol = @(x,y) abs(x - y) <= tol;

        % Separate calculations for real and imaginary components
        iseq = abstol(real(A), real(B)) & abstol(imag(A), imag(B));
end
