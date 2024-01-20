function [maciek_model] = toy_model_isotopomer_creation_100622()
 

maciek_model = sbiomodel('maciek_model')
 
%Compartments
cell = addcompartment(maciek_model,'cell',1,'ConstantCapacity',0,'CapacityUnit','liter') 
media = addcompartment(maciek_model,'media',1,'ConstantCapacity',0,'CapacityUnit','liter')

load metabolite_vec

metabolite_vec2 = {};
for i = 1:8
    curr_met = metabolite_vec{i};
    curr_met = replace(curr_met, '-', '_');
    met_string = strcat(curr_met, '_media');
    eval_string = strcat(met_string, ' = addspecies(media, ''', curr_met, ''' );');
    metabolite_vec2{i,1} = met_string;
    eval(eval_string);
end

for i = 9:28
    curr_met = metabolite_vec{i};
    curr_met = replace(curr_met, '-', '_');
    met_string = strcat(curr_met, '_cell');
    metabolite_vec2{i,1} = met_string;
    eval_string = strcat(met_string, ' = addspecies(cell, ''', curr_met, ''' );');
    eval(eval_string);
end

for i = 29:38
    curr_met = metabolite_vec{i};
    curr_met = replace(curr_met, '-', '_');
    met_string = strcat(curr_met, '_media');
    eval_string = strcat(met_string, ' = addspecies(media, ''', curr_met, ''' );');
    metabolite_vec2{i,1} = met_string;
    eval(eval_string);
end

load ws_k_vec 
k_vec = ws_k_vec

k_1 = addparameter(maciek_model, 'k_1', 100, 'ValueUnits', 'mole/hour', 'Constant', false);
k_2 = addparameter(maciek_model, 'k_2', k_vec(2), 'ValueUnits', '1/hour');
k_3 = addparameter(maciek_model, 'k_3', k_vec(3), 'ValueUnits', '1/hour');
k_4 = addparameter(maciek_model, 'k_4', k_vec(4), 'ValueUnits', '1/hour');
k_5 = addparameter(maciek_model, 'k_5', k_vec(5), 'ValueUnits', '1/mole*hour');
k_6 = addparameter(maciek_model, 'k_6', k_vec(6), 'ValueUnits', '1/hour');

metabolite_vec3 = {};
for i =1:length(metabolite_vec2)
    curr_met = metabolite_vec2{i};
    k1 = findstr(curr_met, '_media');
    k2 = findstr(curr_met, '_cell');
    
    new_met = {};
    if isempty(k1)~=1
        new_met = strcat('media.', curr_met(1:k1-1));
    end
    
    if isempty(k2)~=1
        new_met = strcat('cell.', curr_met(1:k2-1));
    end
    metabolite_vec3{i,1} = new_met;
end

load S_matrix S_matrix

rxn_string_vec = {};    
    
for i = 1:72
    substrates = find(S_matrix(:,i)<0);
    
    substrate_string = '';
    for j = 1:length(substrates)
        if j<length(substrates)
            substrate_string = strcat(substrate_string, strcat(metabolite_vec3{substrates(j)}, ' +', {' '}));
        else
            substrate_string = strcat(substrate_string, strcat(metabolite_vec3{substrates(j)}));  
        end
    end
    
    
    
    products = find(S_matrix(:,i)>0);
    
    product_string = '';
    for j = 1:length(products)
        if j<length(products)
            product_string = strcat(product_string, strcat(metabolite_vec3{products(j)}, ' +', {' '}));
        else
            product_string = strcat(product_string, strcat(metabolite_vec3{products(j)})); 
        end
    end
    
    
    rxn_string = strcat(substrate_string, ' ->', {' '}, product_string);
    rxn_string_vec{i,1} = char(rxn_string);
end





%R1  B uptake

uptake_1 = addparameter(maciek_model, 'uptake_1', 100, 'ValueUnits', 'mole/hour', 'Constant', false);
uptake_2 = addparameter(maciek_model, 'uptake_2', 0, 'ValueUnits', 'mole/hour', 'Constant', false);
uptake_3 = addparameter(maciek_model, 'uptake_3', 0, 'ValueUnits', 'mole/hour', 'Constant', false);
uptake_4 = addparameter(maciek_model, 'uptake_4', 0, 'ValueUnits', 'mole/hour', 'Constant', false);
uptake_5 = addparameter(maciek_model, 'uptake_5', 0, 'ValueUnits', 'mole/hour', 'Constant', false);
uptake_6 = addparameter(maciek_model, 'uptake_6', 0, 'ValueUnits', 'mole/hour', 'Constant', false);
uptake_7 = addparameter(maciek_model, 'uptake_7', 0, 'ValueUnits', 'mole/hour', 'Constant', false);
uptake_8 = addparameter(maciek_model, 'uptake_8', 0, 'ValueUnits', 'mole/hour', 'Constant', false);

for i = 1:8
    eval_string1 = '';
    eval_string1 = strcat(eval_string1, 'R', num2str(i), '_obj = addreaction(maciek_model,''');
    eval_string1 = strcat(eval_string1, rxn_string_vec{i}, "');");
    eval_string1 = char(eval_string1);
    
    eval_string2 = '';
    eval_string2 = strcat(eval_string2, 'Rate_R', num2str(i), ' = addparameter(maciek_model, ''');
    eval_string2 = strcat(eval_string2, 'Rate_R', num2str(i), "'");
    eval_string2 = strcat(eval_string2, ", 'Constant', false);")
    
    eval_string3 = '';
    eval_string3 = strcat(eval_string3, "Rate_R", num2str(i), " = addrule(maciek_model,");
    eval_string3 = strcat(eval_string3, "'Rate_R", num2str(i), " = uptake_", num2str(i),"', 'repeatedAssignment');");
    
    eval_string4 = '';
    eval_string4 = strcat(eval_string4, "set(R", num2str(i), "_obj, 'ReactionRate', 'Rate_R", num2str(i), "');");
    
    eval(eval_string1);
    eval(eval_string2);
    eval(eval_string3);
    eval(eval_string4);
    
end

%R2

for i = 9:16
    eval_string1 = '';
    eval_string1 = strcat(eval_string1, 'R', num2str(i), '_obj = addreaction(maciek_model,''');
    eval_string1 = strcat(eval_string1, rxn_string_vec{i}, "');");
    eval_string1 = char(eval_string1);
    
    eval_string2 = '';
    eval_string2 = strcat(eval_string2, 'Rate_R', num2str(i), ' = addparameter(maciek_model, ''');
    eval_string2 = strcat(eval_string2, 'Rate_R', num2str(i), "'");
    eval_string2 = strcat(eval_string2, ", 'Constant', false);");
    
    
    curr_reaction = rxn_string_vec{i};
    rxn_split = strsplit(curr_reaction, '->');
    substrate = strip(rxn_split{1});
    
    eval_string3 = '';
    eval_string3 = strcat(eval_string3, "Rate_R", num2str(i), " = addrule(maciek_model,");
    eval_string3 = strcat(eval_string3, "'Rate_R", num2str(i), " = k_2*", substrate, "', 'repeatedAssignment');");
    
    eval_string4 = '';
    eval_string4 = strcat(eval_string4, "set(R", num2str(i), "_obj, 'ReactionRate', 'Rate_R", num2str(i), "');");
    eval(eval_string1);
    eval(eval_string2);
    eval(eval_string3);
    eval(eval_string4);
    
end


%R3

for i = 17:24
    eval_string1 = '';
    eval_string1 = strcat(eval_string1, 'R', num2str(i), '_obj = addreaction(maciek_model,''');
    eval_string1 = strcat(eval_string1, rxn_string_vec{i}, "');");
    eval_string1 = char(eval_string1);
    
    eval_string2 = '';
    eval_string2 = strcat(eval_string2, 'Rate_R', num2str(i), ' = addparameter(maciek_model, ''');
    eval_string2 = strcat(eval_string2, 'Rate_R', num2str(i), "'");
    eval_string2 = strcat(eval_string2, ", 'Constant', false);");
    
    
    curr_reaction = rxn_string_vec{i};
    rxn_split = strsplit(curr_reaction, '->');
    substrate = strip(rxn_split{1});
    
    eval_string3 = '';
    eval_string3 = strcat(eval_string3, "Rate_R", num2str(i), " = addrule(maciek_model,");
    eval_string3 = strcat(eval_string3, "'Rate_R", num2str(i), " = k_3*", substrate, "', 'repeatedAssignment');");
    
    eval_string4 = '';
    eval_string4 = strcat(eval_string4, "set(R", num2str(i), "_obj, 'ReactionRate', 'Rate_R", num2str(i), "');");
    
    eval(eval_string1);
    eval(eval_string2);
    eval(eval_string3);
    eval(eval_string4);
    
end



%R4

for i = 25:32
    eval_string1 = '';
    eval_string1 = strcat(eval_string1, 'R', num2str(i), '_obj = addreaction(maciek_model,''');
    eval_string1 = strcat(eval_string1, rxn_string_vec{i}, "');");
    eval_string1 = char(eval_string1);
    
    eval_string2 = '';
    eval_string2 = strcat(eval_string2, 'Rate_R', num2str(i), ' = addparameter(maciek_model, ''');
    eval_string2 = strcat(eval_string2, 'Rate_R', num2str(i), "'");
    eval_string2 = strcat(eval_string2, ", 'Constant', false);");
    
    
    curr_reaction = rxn_string_vec{i};
    rxn_split = strsplit(curr_reaction, '->');
    substrate = strip(rxn_split{1});
    
    eval_string3 = '';
    eval_string3 = strcat(eval_string3, "Rate_R", num2str(i), " = addrule(maciek_model,");
    eval_string3 = strcat(eval_string3, "'Rate_R", num2str(i), " = k_4*", substrate, "', 'repeatedAssignment');");
    
    eval_string4 = '';
    eval_string4 = strcat(eval_string4, "set(R", num2str(i), "_obj, 'ReactionRate', 'Rate_R", num2str(i), "');");
    
    eval(eval_string1);
    eval(eval_string2);
    eval(eval_string3);
    eval(eval_string4);
    
end

%R5

for i = 33:64
    eval_string1 = '';
    eval_string1 = strcat(eval_string1, 'R', num2str(i), '_obj = addreaction(maciek_model,''');
    eval_string1 = strcat(eval_string1, rxn_string_vec{i}, "');");
    eval_string1 = char(eval_string1);
    
    eval_string2 = '';
    eval_string2 = strcat(eval_string2, 'Rate_R', num2str(i), ' = addparameter(maciek_model, ''');
    eval_string2 = strcat(eval_string2, 'Rate_R', num2str(i), "'");
    eval_string2 = strcat(eval_string2, ", 'Constant', false);");
    
    
    curr_reaction = rxn_string_vec{i};
    rxn_split = strsplit(curr_reaction, '->');
    substrates = strip(rxn_split{1});
    substrate_split = strsplit(substrates, '+');
    substrate_1 = strip(substrate_split{1});
    substrate_2 = strip(substrate_split{2});
    
    eval_string3 = '';
    eval_string3 = strcat(eval_string3, "Rate_R", num2str(i), " = addrule(maciek_model,");
    eval_string3 = strcat(eval_string3, "'Rate_R", num2str(i), " = k_5*", substrate_1,"*",substrate_2, "', 'repeatedAssignment');");
    
    eval_string4 = '';
    eval_string4 = strcat(eval_string4, "set(R", num2str(i), "_obj, 'ReactionRate', 'Rate_R", num2str(i), "');");
    
    eval(eval_string1);
    eval(eval_string2);
    eval(eval_string3);
    eval(eval_string4);
    
end

%R6
for i = 65:72
    eval_string1 = '';
    eval_string1 = strcat(eval_string1, 'R', num2str(i), '_obj = addreaction(maciek_model,''');
    eval_string1 = strcat(eval_string1, rxn_string_vec{i}, "');");
    eval_string1 = char(eval_string1);
    
    eval_string2 = '';
    eval_string2 = strcat(eval_string2, 'Rate_R', num2str(i), ' = addparameter(maciek_model, ''');
    eval_string2 = strcat(eval_string2, 'Rate_R', num2str(i), "'");
    eval_string2 = strcat(eval_string2, ", 'Constant', false);");
    
    
    curr_reaction = rxn_string_vec{i};
    rxn_split = strsplit(curr_reaction, '->');
    substrate = strip(rxn_split{1});
    
    eval_string3 = '';
    eval_string3 = strcat(eval_string3, "Rate_R", num2str(i), " = addrule(maciek_model,");
    eval_string3 = strcat(eval_string3, "'Rate_R", num2str(i), " = k_6*", substrate, "', 'repeatedAssignment');");
    
    eval_string4 = '';
    eval_string4 = strcat(eval_string4, "set(R", num2str(i), "_obj, 'ReactionRate', 'Rate_R", num2str(i), "');");
    
    eval(eval_string1);
    eval(eval_string2);
    eval(eval_string3);
    eval(eval_string4);
    
end

set(A_000_media, 'InitialAmount', 10^10);
set(A_000_media, 'InitialAmountUnits', 'mole/liter');

set(B_000_cell, 'InitialAmount', 10000);
set(B_000_cell, 'InitialAmountUnits', 'mole/liter');

set(C_00_cell, 'InitialAmount', 10000);
set(C_00_cell, 'InitialAmountUnits', 'mole/liter');

set(D_0_cell, 'InitialAmount', 10000);
set(D_0_cell, 'InitialAmountUnits', 'mole/liter');


e1 = addevent(maciek_model, 'time>10000', 'uptake_1 = 0');
e2 = addevent(maciek_model, 'time>10000', 'uptake_3 = 100');

