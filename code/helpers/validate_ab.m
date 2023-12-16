function [ab, recep] = validate_ab(ab, recep)

% VALIDATE_AB	Validate antibody and receptor inputs to functions.
%
%	Validates that the number of antibodies is equal to the number of
%	receptor sets provided. Converts receptors to a cell array for correct
%	indexing in functions that use antibody and receptor input.
%
%	USAGE:
%		[AB, RECEP] = VALIDATE_AB(AB, RECEP)
%
%	INPUT:
%		AB = a string vector with the antibodies to use in the complexes
%
%       RECEP = a string or cell array of string vectors with the
%       receptor(s) to use in the complexes; the rows of the array should
%       correspond to the antibodies in AB
%
%	OUTPUT:
%		AB = a column string vector with the antibodies to use in the
%		complexes
%
%       RECEP = a cell array of string vectors where the rows correspond to
%       the receptor partners of the rows of AB
%
%	NOTES:
%		Throws an error if the input sizes do not match. Used for argument
%		validation for functions that take input of antibodies and
%		receptors.
%
%	See also ASSEMBLY, CREATE_BINDING, CREATE_MOLECULES, CREATE_PARAMS.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    ab (:,1) string
    recep {mustBeA(recep, ["string", "cell"])}
end

% -------------------- ERRORS --------------------
% Antibodies and receptors must be same length
if length(ab) == 1 && class(recep) == "cell" && any(size(recep) > 1)
    eid = 'Size:notEqual';
    msg = 'Too many receptor partners provided; only one antibody given.';
    throwAsCaller(MException(eid,msg))
elseif length(ab) == 1 && ~any(size(recep) == 1)
    eid = 'Size:notEqual';
    msg = 'Too many receptor partners provided; only one antibody given.';
    throwAsCaller(MException(eid,msg))
elseif length(ab) > 1 && ~any(size(recep) == length(ab))
    eid = 'Size:notEqual';
    msg = 'Antibody and receptor length are not equal.';
    throwAsCaller(MException(eid,msg))
end

% -------------------- CONVERSION --------------------
% Convert receptors to cell array for correct indexing
if length(ab) == 1 && class(recep) == "string"
    recep = {recep};
elseif class(recep) == "string"

    % Reshape so the rows correspond to different antibodies
    if size(recep, 1) ~= length(ab)
        recep = recep';
    end

    % Convert to a cell array with each antibody in a different cell
    rowdim = ones(1, size(recep, 1));
    recep = mat2cell(recep, rowdim);
end

% Make receptors in each cell a row vector
recep = cellfun(@(x) reshape(x, 1, []), recep, 'UniformOutput', false);

% Sort receptors in each cell in order
recep = cellfun(@(x) sort(x), recep, 'UniformOutput', false);

% Make receptors a column cell vector
recep = reshape(recep, [], 1);

% Remove leading "IL" from receptor names
recep = cellfun(@(x) regexprep(x, "^IL", ""), ...
    recep, 'UniformOutput', false);
