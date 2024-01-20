function [three_comp_model] = Three_Component_Dynamics_101122_652pm()

three_comp_model = sbiomodel('three_comp_model');

%compartments
cell = addcompartment(three_comp_model, 'cell', 1, 'ConstantCapacity', 0, 'CapacityUnit', 'liter');

A_cell = addspecies(cell, 'A');
B_cell = addspecies(cell, 'B');
C_cell = addspecies(cell, 'C');

set(A_cell, 'InitialAmount', 100);
set(A_cell, 'InitialAmountUnits', 'mole/liter');

set(B_cell, 'InitialAmount', 20);
set(B_cell, 'InitialAmountUnits', 'mole/liter');

set(C_cell, 'InitialAmount', 0);
set(C_cell, 'InitialAmountUnits', 'mole/liter');

k_vec = [0.1;0.05];

%Define the reactions 

R1_obj = addreaction(three_comp_model, 'cell.A -> cell.B');
k_1 = addparameter(three_comp_model, 'k_1', k_vec(1), 'ValueUnits', '1/hour');
Rate_R1 = addparameter(three_comp_model, 'Rate_R1', 'Constant', false);
Rate_R1 = addrule(three_comp_model, 'Rate_R1 = k_1*cell.A', 'repeatedAssignment');
set(R1_obj, 'ReactionRate', 'Rate_R1');



R2_obj = addreaction(three_comp_model, 'cell.B -> cell.C');
k_2 = addparameter(three_comp_model, 'k_2', k_vec(2), 'ValueUnits', '1/hour');
Rate_R2 = addparameter(three_comp_model, 'Rate_R2', 'Constant', false);
Rate_R2 = addrule(three_comp_model, 'Rate_R2 = k_2*cell.B', 'repeatedAssignment');
set(R2_obj, 'ReactionRate', 'Rate_R2');



