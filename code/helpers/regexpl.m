function matches = regexpl(str, expression, varargin)

% REGEXPL	Return a logical vector of regular expression matches.
%
%	Uses REGEXP to find regular expression matches. Returns a logical
%	vector of matches instead of the location of the start of the match.
%
%	USAGE:
%		MATCHES = REGEXPL(STR, EXPRESSION, OPTION1, ..., OPTIONM)
%
%	INPUT:
%		STR = character vector, cell array, or string vector containing
%		input text to be matched
%
%       EXPRESSION = character vector, cell array, or string vector
%       containing regular expression pattern to use for matching; if an
%       array of multiple patterns is provided, they will be joined by OR
%
%       OPTION = character vector with search options, passed to the
%       REGEXP function
%
%	OUTPUT:
%		MATCHES = logical array of matches of EXPRESSION in STR
%
%	NOTES:
%		Emulates the behavior of grepl in R. Case sensitive (uses REGEXP
%		instead of REGEXPI) by default, and 'ignorecase' can be used to
%       get case insensitive results.
%
%       Multiple patterns are joined by OR ("|") so a value in STR that
%       matches any provided EXPRESSION will return TRUE.
%
%	See also REGEXP, REGEXPI.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    str string
    expression string
end

% Necessary because varargin represents a variable number of inputs
arguments (Repeating)
    varargin
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Regular Expression Matches

% If multiple patterns were provided, join them with "|"
expression = join(expression, "|");

% Use options if they were provided
if ~isempty(varargin)

    % Return location of positive matches
    matches = regexp(str, expression, varargin{:});
else

    % Return location of positive matches
    matches = regexp(str, expression);
end

% Convert matches to cell if it was not already a cell array
if class(matches) ~= "cell"
    matches = {matches};
end

% Convert positive matches to a logical array
matches = ~cellfun('isempty', matches);

end
