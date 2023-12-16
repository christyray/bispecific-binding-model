function str = make_unique(str)

% MAKE_UNIQUE	Make all values in a string vector unique.
%
%	Appends integers to the end of each non-unique string to make each
%	value in a string vector unique. Can be used to make column names for
%	tables from a vector with duplicate strings.
%
%	USAGE:
%		STR = MAKE_UNIQUE(STR)
%
%	INPUT:
%		STR = a string vector with any number of duplicated values
%
%	OUTPUT:
%		STR = a string vector with integers appended to the end of all
%       duplicated values
%
%	NOTES:
%		MATLAB requires that all column names in a table are non-duplcates.
%       This can be cumbersome when the columns are an intermediate step
%       (e.g., if the table will eventually be stacked by those columns),
%       so this function appends integers to the string values to make them
%       trivially unique.
%
%	See also JOIN, PASTE, COMBINE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    str (:,1) string
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make Unique Values

% Assign an index to each unique value
[~, ~, idxu] = unique(str, "stable");
% uniq(idxu) returns original vector

% Select all indices that represent duplicates in the original string
nuniq = 1:max(idxu);    % idxu has an index for each unique value
idx_repeats = nuniq(histcounts(idxu) > 1);
% Determines which values in idxu represent duplicate values in the
% original string (e.g, if the first, third, and fourth unique values in
% STR are duplicated, idx_repeats will be [1 3 4])
% histcounts returns the count of each value in idxu

% For each duplicated string
for i = 1:length(idx_repeats)

    % Find the indices of the duplicated string in the original string
    repeat_i = find(idxu == idx_repeats(i));

    % Determine how many zeros to pad the number with based on the number
    % of duplicate values (to ensure that the values still sort the same)
    ndup = length(repeat_i);
    ndigit = floor(log10(ndup)) + 1;

    % Make the vector of numbers to append to the strings
    id = pad(string(1:ndup), ndigit, 'left', "0")';

    % Replace the values in the original string
    str(repeat_i) = join([str(repeat_i), id], "");
end
