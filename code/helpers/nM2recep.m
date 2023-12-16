function num = nM2recep(num, cells, volume, params)

% NM2RECEP	Convert units from nM to # receptors/cell.
%
%	Converts units from nM to # receptors/cell based on the volume and the
%   cells present.
%
%	USAGE:
%		NUM = NM2RECEP(NUM, CELLS, VOLUME, UNITS="UNITS")
%
%	INPUT:
%		NUM = the number to be converted, in # nM
%
%       CELLS = the number of cells present in the system
%
%       VOLUME = the volume of the available space in mL or in the units
%       specified by the UNITS argument
%
%       UNITS = the unit for the VOLUME argument; defaults to "mL"
%
%	OUTPUT:
%		NUM = the input NUM converted to # receptors/cell
%
%	NOTES:
%		This function and RECEP2NM can be used to calculate the alpha
%		parameter for the binding model unit conversion based on the
%		experimental conditions.
%
%	See also RECEP2NM, SETUP_CASE, BINDING_MAIN.

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

% # recep/cell = (nM) * (# molecules / mol) * (mL) * (1 / # cells) * ...
% (mol / 1e9 nmol) * (1 L / 1e3 mL)
num = num * 6.022e23 * volume / cells / 1e12;

end
