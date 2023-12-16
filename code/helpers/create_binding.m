function [binding, stoich] = create_binding(ab, recep)

% CREATE_BINDING	Create array of all possible binding pairs.
%
%	Create an array of all binding partners for each given antibody and set
%   of receptors, and an array of the corresponding stoichiometry. Outputs
%   two arrays: (1) a string array where each row is a binding interaction
%   and each column is a molecule in the binding interaction, and (2) a
%   matrix with the on and off stoichiometry for each binding interaction.
%
%	USAGE:
%		[BINDING, STOICH] = CREATE_BINDING(AB, RECEP)
%
%	INPUT:
%		AB = a string vector with the antibodies to use in the complexes
%
%       RECEP = a cell array of string vectors with the receptor(s) to use
%       in the complexes; the rows of the array should correspond to the
%       antibodies in AB
%
%	OUTPUT:
%		BINDING = a string array where each row is a possible binding
%       interaction; the first column is the existing complex and the
%       second column is the next receptor being bound
%
%       STOICH = a matrix where each row corresponds to the binding
%       interactions in BINDING; the first column is the stoichiometry of
%       the association step and the second column is the stoichiometry of
%       the dissociation step
%
%	NOTES:
%		CREATE_BINDING outputs all possible binding interactions - not the
%		required parameters or the final molecules formed. The output from
%		it can be further used to generate those outputs using the
%       CREATE_MOLECULES function.
%
%       The calculation for STOICH is based on the number of each receptor
%       in RECEP.
%
%	See also COMBINE, PASTE, CREATE_MOLECULES.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    ab string
    recep {mustBeA(recep, ["string", "cell"])}
end

% Validate that the antibodies and receptors are the same length and
% convert the receptors to a cell array for correct indexing
[ab, recep] = validate_ab(ab, recep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Antibody Loop

% Initialize output arrays
binding = [];
stoich = [];

% Find all complexes for each antibody
for j = 1:length(ab)
    ab_j = ab(j);
    recep_j = recep{j};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization

% Sort the receptors to ensure consistent output
recep_j = sort(recep_j);

% Depth is the number of rounds of binding to combine
depth = length(recep_j);

% Total is the length of the output from the combinations
total = sum(depth.^(1:depth));

% Initialize combination output and counter
binding_j = strings(total, 2);
idx = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Possible Complexes

% For each number of molecules already in complex
for bound = 1:depth

    % If only one molecule in complex, just the molecule
    if bound == 1
        binding_i = table2array(combine(ab_j, recep_j));

    % If more than one molecule already in complex, start with the
    % previously bound complexes
    else
        binding_i = paste(binding_i);
        binding_i = table2array(combine(binding_i, recep_j));
    end

    % Sort the complexes based on the binding partner
    binding_i = sortrows(binding_i, 2);

    % Add complexes to combined output and increment counter
    endrow = idx + size(binding_i, 1) - 1;
    binding_j(idx:endrow, :) = binding_i;
    idx = endrow + 1;

    % Remove any existing complexes that are identical
    split_r = split(paste(binding_i), "_", 2);
    if contains(ab_j, "_")
        split_r = split_r(:,3:end);
    else
        split_r = split_r(:,2:end);
    end
    split_r_uniq = unique(sort(split_r,2), 'rows');
    binding_i = binding_i(ismember(split_r,split_r_uniq,'rows'),:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove Extraneous Complexes

% Remove extra empty rows from over-allocating
binding_j = binding_j(~all(binding_j == "", 2),:);

% Remove duplicate rows
binding_j = unique(binding_j, 'rows', 'stable');

% Get unique receptors
uniq_recep = unique(recep_j);

% Make array to hold receptor counts
recep_count = [uniq_recep', zeros(length(uniq_recep),1)];

% For each unique receptor, count number of occurrences and remove any
% complexes with too any of a particular receptor
for i = 1:length(uniq_recep)
    recep_i = uniq_recep(i);

    % Count occurrences of receptor in input
    total = sum(matches(recep_j, recep_i));
    recep_count(i,2) = string(total);

    % Count occurrences of receptor in complexes array
    n_recep = count(binding_j, recep_i);

    % Remove any complexes with too many of one receptor
    binding_j = binding_j(sum(n_recep, 2) <= total, :);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine Stoichiometry

% Initialize output
stoich_j = zeros(size(binding_j, 1), 2);

% Iterate over each binding interaction
for i = 1:size(binding_j, 1)
    binding_i = binding_j(i,1);
    recep_i = binding_j(i,2);

    % Find occurrences of receptor in input
    total = str2double(recep_count(recep_count(:,1) == recep_i,2));

    % Determine stoichiometry based on number of receptors in binding
    % interaction and total number of binding sites
    off = count(binding_i, recep_i) + 1;
    on = total - off + 1;

    stoich_j(i,:) = [on, off];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Concatenate Output

binding = cat(1, binding, binding_j);
stoich = cat(1, stoich, stoich_j);

end
