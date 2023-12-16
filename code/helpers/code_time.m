function [elapsed, current] = code_time(start_tic)

% CODE_TIME Determine time elapsed and current system time.
%
%	Calculate the time elapsed since a timer was started and the current
%	system time, and return the values in a formatted string.
%
%	USAGE:
%		[ELAPSED, CURRENT] = CODE_TIME(START_TIC)
%
%	INPUT:
%		START_TIC = a timer initiated with the TIC function
%
%	OUTPUT:
%		ELAPSED = the time since the TIC function was called as a string in
%		HH:mm:ss format
%
%       CURRENT = the current system time as a string in HH:mm:ss format
%
%	NOTES:
%		Designed to be used in printing formatted time strings, such as in
%		code progress messages.
%
%	See also TIC, TOC, DATETIME.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function Arguments

arguments
    start_tic uint64 = uint64(0)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Format Time

% toc gives the time elapsed since tic was called; if toc is given a
% specific start tic, it calculates the time elapsed since that tic
elapsed = toc(start_tic);

% Format the elapsed time in HH:mm:ss to match the current time format
hh = floor(elapsed / 3600);
mm = floor((elapsed - hh * 3600)/60);
ss = round(elapsed - hh * 3600 - mm * 60);

elapsed = compose("%02d:%02d:%02d", hh, mm, ss);

% Determine the current system time and format it in HH:mm:ss
current = string(datetime("now", "Format", "HH:mm:ss"));

end