function out_names = replace_names(param_names, type)

% REPLACE_NAMES	Replace the parameter names with regular expressions.
%
%	Replaces the short-hand parameter names with the matching regular
%	expressions or replaces the regular expressions with the shorthand
%	names, depending on the input type.
%
%	USAGE:
%		OUT_NAMES = REPLACE_NAMES(PARAM_NAMES, TYPE)
%
%	INPUT:
%		PARAM_NAMES = string vector of names or regular expressions
%
%       TYPE = string option for which names were input; options are
%       "short" or "regex"; default is to determine type from the input
%       names
%
%	OUTPUT:
%		OUT_NAMES = string vector of names or regular expressions
%
%	NOTES:
%		Goal is to be able to automatically switch between the short-hand
%		names and the necessary regular expressions.
%
%	See also NAMES.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    param_names (1,:) string
    type string {mustBeMember(type, ["short", "regex"])} = strings(0)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Replace Names

% If the type was not given, determine the type from the input names
if isempty(type)
    if any(regexpl(param_names, "Ab"))
        type = "regex";
    else
        type = "short";
    end
end

% Define the possible short names and their corresponding regular
% expressions
short_names = ["kon6R", "kon8R", "kon6Rprime", "kon8Rprime", ...
    "koff", "koff6R", "koff8R"];
expressions = ["kon_Ab_6R$", "kon_Ab_8R$", "kon_Ab.*_[68]R_6R$", ...
    "kon_Ab.*_[68]R_8R$", "koff_Ab", "koff_Ab.*_6R$", "koff_Ab.*_8R$"];

% Set the names to match based on which names were input
switch type
    % Match the regular expressions to the short names
    case "regex"
        matches = expressions;
        outputs = short_names;

    % Match the short names to their regular expressions
    case "short"
        matches = short_names;
        outputs = expressions;
end

% Second output of ismember is the index of the match or 0 if there
% was not a match
[~, idx] = ismember(param_names, matches);

% Throw an error if any of the input names did not match, otherwise
% output the matched names
if any(idx == 0)
    missing = join(param_names(idx == 0), ", ");
    msg = "Parameter(s) missing from replacement list: %s";
    eid = "Inputs:noMatchesFound";
    error(eid, msg, missing);
else
    out_names = outputs(idx);
end
