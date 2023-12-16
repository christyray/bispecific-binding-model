function [bal, bal_percent] = mole_balance(Y, net, m, p, abs_tol, rel_tol)

% MOLE_BALANCE	Calculate mole balance from ODE solver outputs.
%
%	Calculate the mole balance from the amount currently in the system and
%	the net amount added to the system. Output the balance for each species
%	in a matrix and print a warning if the moles are not balanced.
%
%	USAGE:
%		[BAL, BAL_PERCENT] = MOLE_BALANCE(Y, NET, M, P, ABS_TOL, REL_TOL)
%
%	INPUT:
%		Y = matrix with the concentrations of each species at each time
%		point; output from the ODE solver
%
%       NET = matrix with the net concentrations of each species added into
%       the system at each time point; calculated in the main function
%
%       M = structure containing the molecules included in the model and
%       their assigned numbers, input as a structure, correct naming for
%       species is antibody name or antibody receptor complex name (with
%       components separated by underscores)
%
%       P = structure containing parameters for the model, input as a
%       structure, correct naming for parameters is "kon_" or "koff_",
%       followed by the antibody name or antibody receptor complex name
%       (with components separated by underscores), ending with the
%       receptor being bound to; e.g., kon_BS1_6R_8R
%
%       ABS_TOL = absolute tolerance for the mole balance calculation
%       warning; mole balance values above this will prompt a warning to
%       display; default is 1e-12
%
%       REL_TOL = relative tolerance for the mole balance calculation
%       warning; calculated as the percentage of the net amount of the
%       molecule added to the system that is out of balance; percentages
%       above this will prompt a warning to display; default is 1e-6
%
%	OUTPUT:
%		BAL = matrix with the mole balance for each species at each
%		time point
%
%		BAL_PERCENT = matrix with mole balance as a percentage of the net
%       amount added for each species at each time point
%
%	NOTES:
%		It is not required to call MOLE_BALANCE with an output since it
%       will display a warning to the command window if the moles of the
%       components are not balanced.
%
%       Mole balance is successful if the balance is within either the
%       absolute tolerance or the relative tolerance; being within both
%       tolerances is not required.
%
%	See also BINDING_MAIN, COMPLEX_CREATE, MOLECULE_CREATE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    Y double
    net double
    m struct
    p struct
    abs_tol double = 1e-12
    rel_tol double = 1e-6
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mole Balance

% Load complexes from the saved .mat file
load complexes.mat bound singles

% Remove "IL" from the receptors for later regexp
singles = regexprep(singles, "^IL([0-9]R)$", "$1");

% Initialize output matrix
bal = zeros(size(Y,1), length(singles));
net_total = zeros(size(Y,1), length(singles));

% For each individual molecule
for i = 1:length(singles)
    molecule = singles(i);

    % Find number of occurrences of molecule in each complex
    % cellfun only works on cell arrays or arrays with length > 1
    if length(bound) == 1
        n = length(regexp(bound, molecule));
    else
        n = cellfun('length', regexp(bound, molecule));
    end
    matches = bound(n > 0);
    n = n(n > 0);

    % Get positions of matching complexes in input matrices
    idx_com = cellfun(@(name) m.(name), matches);

    % Calculate the mole balance for the input molecules
    if strcmp(molecule, "6R") || strcmp(molecule, "8R")

        % Get position of individual molecule in input matrices
        idx_mol = cellfun(@(name) m.(name), "IL" + molecule);

        % For input receptors, values are already # recep/cell
        bal_mol = Y(:,idx_mol) - net(:,idx_mol);
        net_mol = net(:, idx_mol);
    else
        % Get position of individual molecule in input matrices
        idx_mol = cellfun(@(name) m.(name), molecule);

        % For input ligands and molecules, convert nM to # recep/cell
        bal_mol = (Y(:,idx_mol) - net(:,idx_mol)) ./ p.alpha;
        net_mol = net(:,idx_mol) ./ p.alpha;
    end

    % Calculate the mole balance for the molecule in the complexes
    bal_com = (Y(:,idx_com) - net(:,idx_com)) .* n';
    net_com = net(:,idx_com) .* n';

    % Combine the complexes and the individual molecules
    bal(:,i) = bal_mol + sum(bal_com,2);
    net_total(:,i) = net_mol + sum(net_com,2);
end

% Calculate maximum error and percent error for each molecule
max_balance = max(abs(bal));
bal_percent = bal ./ net_total;
max_percent = max(abs(bal_percent));
max_percent(isinf(max_percent)) = NaN;
% Replace any Inf from very small divisors with NaN

% Print error message and table if the moles are not balanced
if any(max_balance > abs_tol & max_percent > rel_tol)
    warnID = 'balance:NotBalanced';
    warnMSG = "The moles of the components are not balanced";
    warning(warnID, warnMSG);
    balanceTable = table(singles,max_balance',max_percent', ...
        'VariableNames',["Species", "Balance", "Frac of Net"]);
    disp(balanceTable)
end

end
