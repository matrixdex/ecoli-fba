function [maciek_model] = Toy_Model_Maciek_Dynamics_101122_7pm()

maciek_model = sbiomodel('maciek_model');

%Compartments

cell = addcompartment(maciek_model, 'cell', 1, 'ConstantCapacity', 0, 'CapacityUnit', 'liter');
media = addcompartment(maciek_model, 'media', 1, 'ConstantCapacity', 0, 'CapacityUnit', 'liter');

B_media = addspecies(media, 'B');
B_cell = addspecies(cell, 'B');
C_cell = addspecies(cell, 'C');
D_cell = addspecies(cell, 'D');
E_media = addspecies(media, 'E');
D_media = addspecies(media, 'D');

set(B_cell, 'InitialAmount', 10000);
set(B_cell, 'InitialAmountUnits', 'mole/liter');

set(C_cell, 'InitialAmount', 10000);
set(C_cell, 'InitialAmountUnits', 'mole/liter');

set(D_cell, 'InitialAmount', 10000);
set(D_cell, 'InitialAmountUnits', 'mole/liter');


set(B_media, 'InitialAmount', 10^10);
set(B_media, 'InitialAmountUnits', 'mole/liter');

load ws_k_vec
k_vec = ws_k_vec;

%R1 - B uptake
R1_obj = addreaction(maciek_model, 'media.B -> cell.B');
k_1 = addparameter(maciek_model, 'k_1', k_vec(1), 'ValueUnits', 'mole/hour', 'Constant', false);
Rate_R1 = addparameter(maciek_model, 'Rate_R1', 'Constant', false);
Rate_R1 = addrule(maciek_model, 'Rate_R1 = k_1', 'repeatedAssignment');
set(R1_obj, 'ReactionRate', 'Rate_R1');

%R2 - B to D

R2_obj = addreaction(maciek_model, 'cell.B -> cell.D');
k_2 = addparameter(maciek_model, 'k_2', k_vec(2), 'ValueUnits', '1/hour');
Rate_R2 = addparameter(maciek_model, 'Rate_R2', 'Constant', false);
Rate_R2 = addrule(maciek_model, 'Rate_R2 = k_2*cell.B', 'repeatedAssignment');
set(R2_obj, 'ReactionRate', 'Rate_R2');

%R3 - D to B

R3_obj = addreaction(maciek_model, 'cell.D -> cell.B');
k_3 = addparameter(maciek_model, 'k_3', k_vec(3), 'ValueUnits', '1/hour');
Rate_R3 = addparameter(maciek_model, 'Rate_R3', 'Constant', false);
Rate_R3 = addrule(maciek_model, 'Rate_R3 = k_3*cell.D', 'repeatedAssignment');
set(R3_obj, 'ReactionRate', 'Rate_R3');

%R4 - B-> C + E

R4_obj = addreaction(maciek_model, 'cell.B -> cell.C + media.E');
k_4 = addparameter(maciek_model, 'k_4', k_vec(4), 'ValueUnits', '1/hour');
Rate_R4 = addparameter(maciek_model, 'Rate_R4', 'Constant', false);
Rate_R4 = addrule(maciek_model, 'Rate_R4 = k_4*cell.B', 'repeatedAssignment');
set(R4_obj, 'ReactionRate', 'Rate_R4');

%R5 - B + C -> D + 2E

R5_obj = addreaction(maciek_model, 'cell.B + cell.C -> cell.D + 2 media.E');
k_5 = addparameter(maciek_model, 'k_5', k_vec(5), 'ValueUnits', '1/hour');
Rate_R5 = addparameter(maciek_model, 'Rate_R5', 'Constant', false);
Rate_R5 = addrule(maciek_model, 'Rate_R5 = k_5*cell.B*cell.C', 'repeatedAssignment');
set(R5_obj, 'ReactionRate', 'Rate_R5');

%R6 - D to F

R6_obj = addreaction(maciek_model, 'cell.D -> media.D');
k_6 = addparameter(maciek_model, 'k_6', k_vec(6), 'ValueUnits', '1/hour');
Rate_R6 = addparameter(maciek_model, 'Rate_R6', 'Constant', false);
Rate_R6 = addrule(maciek_model, 'Rate_R6 = k_6*cell.D', 'repeatedAssignment');
set(R6_obj, 'ReactionRate', 'Rate_R6');












