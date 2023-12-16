function num = recep2nM(num, cells, volume, params)

% RECEP2NM	Convert units from # receptors/cell to nM.
%
%	Converts units from # receptors/cell to nM based on the volume and the
%   cells present.
%
%	USAGE:
%		NUM = RECEP2NM(NUM, CELLS, VOLUME, UNITS="UNITS")
%
%	INPUT:
%		NUM = the number to be converted, in # receptors/cell
%
%       CELLS = the number of cells present in the system
%
%       VOLUME = the volume of the available space in mL or in the units
%       specified by the UNITS argument
%
%       UNITS = the unit for the VOLUME argument; defaults to "mL"
%
%	OUTPUT:
%		NUM = the input NUM converted to nM
%
%	NOTES:
%		This function and NM2RECEP can be used to calculate the alpha
%		parameter for the binding model unit conversion based on the
%		experimental conditions.
%
%	See also NM2RECEP, SETUP_CASE, BINDING_MAIN.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Argument Validation

arguments
    num double
    cells double
    volume double
    params.units = "mL"
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convert Volume to mL

% If the volume was given in a different unit, convert to mL
if ~strcmpi(params.units, "mL")
    if strcmpi(params.units, "L")
        volume = volume * 1000;
    elseif strcmpi(params.units, "uL")
        volume = volume / 1000;
    elseif strcmpi(params.units, "nL")
        volume = volume / 1e6;
    else
        error('unit:notFound', ...
            "The specified unit conversion is not included. " + ...
            "Please select a different unit.")
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convert Units

% nM = (# recep / cell) * (mol / # molecules) * (1 / mL) * (# cells) * ...
% (1e9 nmol / mol) * (1e3 mL / 1 L)
num = num / 6.022e23 / volume * cells * 1e12;

end
