function [m, p, cond] = setup_case(setup, receptor, option)

% SETUP_CASE	Define environment values based on experimental setup.
%
%	Generates the m and p structures necessary for model simulations, and
%	calculates the value for the unit conversion parameter p.alpha based on
%	experimental setup used. Also outputs a reference structure with other
%	information about the environment used for the experiment.
%
%	USAGE:
%		[M, P, COND] = SETUP_CASE(SETUP, RECEPTOR, OPTION)
%
%	INPUT:
%		SETUP = experimental design used, must be a string
%
%       RECEPTOR = string vector of receptors present in the system,
%       defaults to ["6R", "8R"]
%
%       OPTION = additional options for output; options are ["nostruct"]
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
%       COND = structure containing values for experimental conditions,
%       including cells, number of receptors per cell, etc.
%
%	NOTES:
%		Most important output from this function is the value for p.alpha -
%		necessary for correct unit conversion between nM and # recep/cell
%		in model simulations.
%
%	See also BINDING_SIM, CREATE_PARAMS, SETUP_INPUT.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    setup string
    receptor string = ["6R", "8R"]
    option string = ""
end

% Remove leading "IL" from receptor names
receptor = regexprep(receptor, "^IL", "");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign Environment Values

% Assign values based on experimental setup
switch setup

    % Flow cytometry experimental setup
    case "flow"
        cond.cells = 1e5;           % cells
        cond.workingVolume = 0.2;   % mL
        cond.primary = 2;           % hours
        cond.secondary = 0.25;      % hours
        cond.ytype = ["initial", "washout"];    % period types

        % Define default receptors value as a reference
        % Receptor numbers from cell measurements
        if any(strcmp(receptor, "6R"), 'all') && ...
                any(strcmp(receptor, "8R"), 'all')
            cond.recep = [3.16e5, 6.18e5];  % 6R+8R+ cells
        elseif any(strcmp(receptor, "6R"), 'all')
            cond.recep = 5.08e5;            % 6R+8R- cells
        elseif any(strcmp(receptor, "8R"), 'all')
            cond.recep = 1.3e6;             % 6R-8R+ cells
        end
        % 6R+ = 5.08e5; 8R+ = 1.3e6; 6R+8R+ = 3.16e5 and 6.18e5

        % Alpha parameter to use in ligand equations for correct units
        cond.alpha = recep2nM(1, cond.cells, cond.workingVolume);
        % Units = nM/(# receptors/cell)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Parameter Structures

% If flag was given to skip structure generation (for performance reasons)
if any(option == "nostruct")
    m = struct();
    p = struct();
else
    % Generates the structures based on complexes.mat
    [m, p] = create_params();
    p.alpha = cond.alpha;   % Assign alpha from experimental conditions
end
