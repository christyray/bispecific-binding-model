function dydt = binding_eqns(~, y, m, p)

% BINDING_EQNS	System of ODEs for mechanistic binding model.
%
%   Sets up the system of ordinary differential equations for the binding
%   of antibodies and ligands to the IL6 and IL8 receptors. Used as an
%   input for an ODE solver.
%
%   USAGE:
%       DYDT = BINDING_EQNS(~, Y, M, P)
%
%   INPUT:
%       T = time, only present because it is necessary for the ODE solver,
%       represented by a tilde (~)
%
%       Y = concentration of each component (in the order of the species
%       list), input as a vector, input needs to be the initial
%       concentrations at the starting time
%
%       M = structure containing the molecules included in the model and
%       their assigned numbers, input as a structure, correct naming for
%       species is antibody name or antibody receptor complex name (with
%       components separated by underscores)
%
%       P = structure containing parameters for the model, input as a
%       structure, correct naming for parameters is 'kon_' or 'koff_',
%       followed by the antibody name or antibody receptor complex name
%       (with components separated by underscores), ending with the
%       receptor being bound to; e.g., kon_BS1_6R_8R
%
%   OUTPUT:
%       DYDT = matrix of concentration values at each time point; columns
%       are in the order that the equations are defined in (same as species
%       list)
%
%   NOTES:
%       Used exclusively as an input to an ODE solver. Sets up the system
%       of binding equations for the IL6R/IL8R system and all antibodies
%       being studied.
%
%   See also BINDING_MAIN, BINDING_SIM.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equation Snippets

v_Toci_6R = 2 * p.kon_Toci_6R * y(m.Toci) * y(m.IL6R) - 1 * p.koff_Toci_6R * y(m.Toci_6R);
v_Toci_6R_6R = 1 * p.kon_Toci_6R_6R * y(m.Toci_6R) * y(m.IL6R) - 2 * p.koff_Toci_6R_6R * y(m.Toci_6R_6R);
v_H2_8R = 2 * p.kon_H2_8R * y(m.H2) * y(m.IL8R) - 1 * p.koff_H2_8R * y(m.H2_8R);
v_H2_8R_8R = 1 * p.kon_H2_8R_8R * y(m.H2_8R) * y(m.IL8R) - 2 * p.koff_H2_8R_8R * y(m.H2_8R_8R);
v_BS1_6R = 1 * p.kon_BS1_6R * y(m.BS1) * y(m.IL6R) - 1 * p.koff_BS1_6R * y(m.BS1_6R);
v_BS1_8R = 1 * p.kon_BS1_8R * y(m.BS1) * y(m.IL8R) - 1 * p.koff_BS1_8R * y(m.BS1_8R);
v_BS1_8R_6R = 1 * p.kon_BS1_8R_6R * y(m.BS1_8R) * y(m.IL6R) - 1 * p.koff_BS1_8R_6R * y(m.BS1_6R_8R);
v_BS1_6R_8R = 1 * p.kon_BS1_6R_8R * y(m.BS1_6R) * y(m.IL8R) - 1 * p.koff_BS1_6R_8R * y(m.BS1_6R_8R);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Binding Equations

dydt = zeros(12,1);

dydt(m.Toci) = p.alpha * (- v_Toci_6R);
dydt(m.H2) = p.alpha * (- v_H2_8R);
dydt(m.BS1) = p.alpha * (- v_BS1_6R - v_BS1_8R);
dydt(m.IL6R) = - v_Toci_6R - v_Toci_6R_6R - v_BS1_6R - v_BS1_8R_6R;
dydt(m.IL8R) = - v_H2_8R - v_H2_8R_8R - v_BS1_8R - v_BS1_6R_8R;
dydt(m.Toci_6R) = + v_Toci_6R - v_Toci_6R_6R;
dydt(m.Toci_6R_6R) = + v_Toci_6R_6R;
dydt(m.H2_8R) = + v_H2_8R - v_H2_8R_8R;
dydt(m.H2_8R_8R) = + v_H2_8R_8R;
dydt(m.BS1_6R) = + v_BS1_6R - v_BS1_6R_8R;
dydt(m.BS1_8R) = + v_BS1_8R - v_BS1_8R_6R;
dydt(m.BS1_6R_8R) = + v_BS1_8R_6R + v_BS1_6R_8R;
