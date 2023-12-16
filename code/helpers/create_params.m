function [m, p] = create_params(ab, recep)

% CREATE_PARAMS	Create the molecule and parameter structures.
%
%	Create the molecule and parameter structures for the mechanistic
%	binding model from input antibodies and their corresponding receptor
%	partners. Accepts input of antibodies and receptor partners, or
%	defaults to reading in the binding pairs from 'antibodies.txt'.
%
%	USAGE:
%		[M, P] = CREATE_PARAMS(AB, RECEP)
%
%	INPUT:
%		AB = a string vector with the antibodies to use in the complexes;
%		uses antibodies listed in 'antibodies.txt' by default
%
%       RECEP = a cell array of string vectors with the receptor(s) to use
%       in the complexes; the rows of the array should correspond to the
%       antibodies in AB; uses the receptor pairs listed in
%       'antibodies.txt' by default
%
%	OUTPUT:
%       M = structure containing the molecules included in the model and
%       their assigned numbers, correct naming for species is antibody name
%       or antibody receptor complex name (with components separated by
%       underscores)
%
%       P = structure containing parameters for the model, correct naming
%       for parameters is "kon_" or "koff_", followed by the antibody name
%       or antibody receptor complex name (with components separated by
%       underscores), ending with the receptor being bound to; e.g.,
%       kon_BS1_6R_8R; all parameters are assigned a value of 0
%
%	NOTES:
%		Outputs the m and p structures necessary for the system of ODEs
%		describing the mechanistic binding model. Uses the 'antibodies.txt'
%		file to determine the antibodies and receptors present in the
%		system if they are not provided, similar to the ASSEMBLY script.
%
%	See also CREATE_BINDING, CREATE_MOLECULES.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    ab string = strings(0)
    recep {mustBeA(recep, ["string", "cell"])} = strings(0)
end

% Validate that the antibodies and receptors are the same length and
% convert the receptors to a cell array for correct indexing
[ab, recep] = validate_ab(ab, recep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read Antibodies and Receptors

% If the antibodies and receptors were not passed as input
if isempty(ab)

    % Set options for reading from antibodies.txt
    text_options = delimitedTextImportOptions('Delimiter', '|', ...
        'Whitespace', '', 'CommentStyle', '#', ...
        'VariableNames', ["ab", "recep"], ...
        'VariableTypes', ["string", "string"]);

    % Read in antibodies and corresponding receptors
    [ab, recep] = readvars('code/simulations/antibodies.txt', ...
        text_options);

    % Convert the receptors to a cell array of string vectors
    recep = cellfun(@(x) convertCharsToStrings(split(x, ",", 2)), ...
        recep, 'UniformOutput', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Bound Complexes, Single Molecules, and Parameters

% Use create_molecules to generate the molecules in the system
[bound, singles, params] = create_molecules(ab, recep);

% Create molecule structure
mnames = [singles; bound];
mnum = (1:numel(mnames))';
m = cell2struct(num2cell(mnum), mnames);

% Create parameter structure
kon = append("kon_", params);
koff = append("koff_", params);
pnames = ["alpha"; kon; koff];
pnum = zeros(numel(pnames),1);
p = cell2struct(num2cell(pnum), pnames);

end
