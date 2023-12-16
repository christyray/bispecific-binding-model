function molecules = select_molecules(ab, recep, interest, bound, singles)

% SELECT_MOLECULES	Select molecules of interest for simulation output.
%
%	Outputs names of molecules of interest for the simulation output
%	calculations based on the antibodies, receptors, and desired
%	calculation to perform.
%
%	USAGE:
%		MOLECULES = SELECT_MOLECULES(AB, RECEP, INTEREST, BOUND, SINGLES)
%
%	INPUT:
%		AB = a string vector with the antibodies to use in the complexes
%
%       RECEP = a cell array of string vectors with the receptor(s) to use
%       in the complexes; the rows of the array should correspond to the
%       antibodies in AB
%
%       INTEREST = a string vector with what type of molecules should be
%		returned; options are ["inputs", "antibodies", "receptors",
%       "binary", "ternary", "present", "all"]
%
%       BOUND = antibody-receptor combinations corresponding to
%       molecules in the m structure; e.g., only includes BS1_6R_8R and not
%       BS1_8R_6R because they are the same molecule
%
%       SINGLES = all of the individual antibodies and receptors that
%       contribute to the bound complexes in BOUND
%
%	OUTPUT:
%		MOLECULES = a string vector or cell array of the molecules
%		corresponding to the molecules of interest requested
%
%	NOTES:
%		"inputs" returns the unbound ligands, antibodies, and receptors.
%
%       "antibodies" returns the antibodies and antibody-receptor
%		complexes (can also be used with the ligands).
%
%       "receptors" returns all antibody-receptor complexes.
%
%       "binary" returns the binary antibody-receptor complexes.
%
%       "ternary" returns the ternary antibody-receptor complexes.
%
%       "present" returns the unbound ligands, antibodies, and receptors
%       and all of the antibody-receptor complexes present in the system.
%
%       "all" returns all of the molecules in the equations file,
%       regardless of the given AB and RECEP.
%
%       If multiple entries are given for the INTEREST argument, the
%       molecules of interest will be generated for all values and the
%       final list will be the unique molecules of interest from all
%       inputs.
%
%       BOUND and SINGLES are output from create_molecules. If no inputs
%       are given, the values from complexes.mat will be loaded.
%
%	See also CREATE_MOLECULES, BINDING_SIM.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    ab (:,1) string
    recep {mustBeA(recep, ["string", "cell"])}
    interest (:,1) string {mustBeMember(interest, ...
        ["inputs", "antibodies", "receptors", "binary", "ternary", ...
        "present", "all"])}
    bound (:,1) string = strings(0)
    singles (:,1) string = strings(0)
end

% Validate that the antibodies and receptors are the same length and
% convert the receptors to a cell array for correct indexing
[ab, recep] = validate_ab(ab, recep);

% Convert the molecules of interest to lower case for matching
interest = lower(interest);

% If complexes were not provided, load from saved .mat file
if isempty(bound) || isempty(singles)
    load complexes.mat bound singles
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select Molecules of Interest

% Replace "present" with the combination of "inputs" and "antibodies"
if any(strcmp(interest, "present"))
    interest = interest(~strcmp(interest, "present"));
    interest = [interest; "inputs"; "antibodies"];
end

% Initialize output array
molecules = cell(length(interest), 1);

% For each type of output requested
for i = 1:length(interest)

    % Generate list of molecules based on molecules requested
    switch interest(i)

        % Unbound ligands, antibodies, and receptors
        case "inputs"

            % Combine all receptors into a single column string vector
            recep_i = reshape(horzcat(recep{:}), [], 1);

            % Select only the specific molecules of interest from singles
            select = regexpl(singles, [ab; recep_i]);
            molecules{i} = singles(select);

        % Antibodies and antibody-receptor complexes
        case {"antibodies", "receptors"}

            % Make a search list of input antibodies with only the
            % receptor partners input into the function

            % Initialize search list
            select = strings(length(bound), 1);
            endrow = 1;

            % For each antibody
            for j = 1:length(ab)

                % For each possible number of bound receptors
                for k = 1:length(recep{j})

                    % Generate a list of all receptor pairs separated by _
                    recep_k = nchoosek(recep{j}, k);
                    if k > 1
                        recep_k = paste(recep_k, "_");
                    end

                    % Select complexes containing the input antibodies and
                    % receptor partners
                    select_k = "^" + ab(j) + "_" + unique(recep_k) + "$";
                    select(endrow:endrow+length(select_k)-1) = select_k;
                    endrow = endrow + length(select_k);
                end
            end

            % Select the complexes matching the search list
            select = regexpl(bound, select);

            % Include the antibodies in the output if requested
            if interest(i) == "antibodies"
                molecules{i} = [ab; bound(select)];
            else
                molecules{i} = bound(select);
            end

        % Binary antibody-receptor complexes
        case "binary"

            % Make a search list of input antibodies with only the
            % receptor partners input into the function

            % Initialize search list
            select = strings(length(bound), 1);
            endrow = 1;

            % For each antibody
            for j = 1:length(ab)

                % Select complexes containing the input antibodies and
                % a single receptor partner from the list
                select_k = "^" + ab(j) + "_" + unique(recep{j}) + "$";
                select(endrow:endrow+length(select_k)-1) = select_k;
                endrow = endrow + length(select_k);
            end

            % Select the complexes matching the search list
            select = regexpl(bound, select);
            molecules{i} = bound(select);

        % Ternary antibody-receptor complexes
        case "ternary"

            % Make a search list of input antibodies with only the
            % receptor partners input into the function

            % Initialize search list
            select = strings(length(bound), 1);
            endrow = 1;

            % For each antibody
            for j = 1:length(ab)

                % Generate a list of all receptor pairs separated by _
                recep_k = nchoosek(recep{j}, 2);
                recep_k = paste(recep_k, "_");

                % Select complexes containing the input antibodies and
                % a single receptor partner from the list
                select_k = "^" + ab(j) + "_" + unique(recep_k) + "$";
                select(endrow:endrow+length(select_k)-1) = select_k;
                endrow = endrow + length(select_k);
            end

            % Select the complexes matching the search list
            select = regexpl(bound, select);
            molecules{i} = bound(select);

        % All molecules present in the equations file
        case "all"
            molecules{i} = [singles; bound];
    end
end

% Combine all molecules into a single column string vector
molecules = reshape(vertcat(molecules{:}), [], 1);

% Remove all duplicates
molecules = unique(molecules, 'stable');

end
