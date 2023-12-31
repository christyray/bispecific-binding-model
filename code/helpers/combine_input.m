function input = combine_input(yin, p)

% COMBINE_INPUT	Combine concentration and parameter tables.
%
%	Combines the initial concentration and parameter tables generated by
%	SETUP_INPUT into a single table for use in BINDING_SIM. Adds an ID
%	column so each simulation has a unique ID.
%
%	USAGE:
%		INPUT = COMBINE_INPUT(YIN, P)
%
%	INPUT:
%		YIN = table with columns for the simulation ID, time period type,
%		starting and ending times for the time period, species name, and
%		initial concentration for that species for each time period;
%		generated by SETUP_INPUT
%
%       P = table with columns for the simulation ID, parameter name, and
%		parameter value for each parameter in the simulation; generated by
%       SETUP_INPUT
%
%	OUTPUT:
%		INPUT = table with columns for the simulation ID, time period type,
%		starting and ending times for the time period, species name,
%		initial concentration, parameter name, and parameter value for each
%		simulation to be performed; gives all possible combinations of the
%		input tables
%
%	NOTES:
%		Generates the table necessary for the INPUT argument of
%		BINDING_SIM. Uses the COMBINE function to generate all possible
%		combinations of the input tables.
%
%	See also COMBINE, SETUP_INPUT, BINDING_SIM.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    yin table
    p table
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Input Table

% Combine input tables together
input = combine(yin, p);

% Generate new ID column for combined input
% Select existing ID columns
id = [input.ConcID, input.ParamID];

% Create the new ID column based on the unique concentration values and
% parameter values
% The third output of unique is the indices that can be used with the
% unique values to generate the original matrix, means that each repeated
% pair of original IDs will have the same value for the new ID as intended
[~, ~, id] = unique(id, 'rows', 'stable');

% Concatenate the new ID onto the orignal table
input = addvars(input, id, 'Before', 1);
input = movevars(input, ["ConcID", "ParamID"], 'After', "id");

% Rename the columns to standardize the input table
input = rename(input, ...
    new=["ID", "ConcID", "ParamID", "Ab", "Recep", "Period", ...
    "Start", "End", "Species", "Conc", "Param", "Value"]);

% Sort the output table by the ID
input = sortrows(input, "ID");

end
