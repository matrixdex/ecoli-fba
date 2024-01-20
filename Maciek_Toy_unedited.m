all_reactions = [];
sp = [];

% V1
v1=[];
s1_matrix = get_isotopomer_matrix(3);

rxn_count = 1;
v1_substrates_products = [];
for i=1:size(s1_matrix,2)
    curr_s1 = s1_matrix(:,i);
    prod_1 = curr_s1;
    v1_substrates_products(rxn_count,1).s1 = curr_s1;
    v1_substrates_products(rxn_count,1).p1 = prod_1;
    
    rxn_count = rxn_count + 1;
end

all_reactions(1,1).substrates_products = v1_substrates_products;
all_reactions(1,1).substrate_names = {'A'};
all_reactions(1,1).product_names = {'B'};

% V2
v2=[];
s1_matrix = get_isotopomer_matrix(3);

rxn_count = 1;
v2_substrates_products = [];
for i=1:size(s1_matrix,2)
    curr_s1 = s1_matrix(:,i);
    prod_1 = curr_s1;
    v2_substrates_products(rxn_count,1).s1 = curr_s1;
    v2_substrates_products(rxn_count,1).p1 = prod_1;
    
    rxn_count = rxn_count + 1;
end

all_reactions(2,1).substrates_products = v2_substrates_products;
all_reactions(2,1).substrate_names = {'B'};
all_reactions(2,1).product_names = {'D'};


% V3
v3=[];
s1_matrix = get_isotopomer_matrix(3);

rxn_count = 1;
v3_substrates_products = [];
for i=1:size(s1_matrix,2)
    curr_s1 = s1_matrix(:,i);
    prod_1 = curr_s1;
    v3_substrates_products(rxn_count,1).s1 = curr_s1;
    v3_substrates_products(rxn_count,1).p1 = prod_1;
    
    rxn_count = rxn_count + 1;
end

all_reactions(3,1).substrates_products = v3_substrates_products;
all_reactions(3,1).substrate_names = {'D'};
all_reactions(3,1).product_names = {'B'};



% V4
v4=[];
s1_matrix = get_isotopomer_matrix(3);

rxn_count = 1;
v4_substrates_products = [];
for i=1:size(s1_matrix,2)
    curr_s1 = s1_matrix(:,i);
    prod_1 = [curr_s1(2);curr_s1(3)];
    prod_2 = [curr_s1(1)];
    v4_substrates_products(rxn_count,1).s1 = curr_s1;
    v4_substrates_products(rxn_count,1).p1 = prod_1;
    v4_substrates_products(rxn_count,1).p2 = prod_2;
    
    rxn_count = rxn_count + 1;
end

all_reactions(4,1).substrates_products = v4_substrates_products;
all_reactions(4,1).substrate_names = {'B'};
all_reactions(4,1).product_names = {'C'; 'E'};


% V5
v5=[];
s1_matrix = get_isotopomer_matrix(3);
s2_matrix = get_isotopomer_matrix(2);

rxn_count = 1;
v5_substrates_products = [];
for i = 1:size(s1_matrix,2)
    for j = 1:size(s2_matrix,2)
        curr_s1 = s1_matrix(:,i);
        curr_s2 = s2_matrix(:,j);
        
        prod_1 = [curr_s1(2);curr_s1(3);curr_s2(1)];
        prod_2 = [curr_s1(1)];
        prod_3 = [curr_s2(2)];
        
        v5_substrates_products(rxn_count,1).s1 = curr_s1;
        v5_substrates_products(rxn_count,1).s2 = curr_s2;
        v5_substrates_products(rxn_count,1).p1 = prod_1;
        v5_substrates_products(rxn_count,1).p2 = prod_2;
        v5_substrates_products(rxn_count,1).p3 = prod_3;
        
        rxn_count = rxn_count + 1;
    end
end

all_reactions(5,1).substrates_products = v5_substrates_products;
all_reactions(5,1).substrate_names = {'B';'C'};
all_reactions(5,1).product_names = {'D';'E';'E'};





% V6
v6=[];
s1_matrix = get_isotopomer_matrix(3);

rxn_count = 1;
v6_substrates_products = [];
for i=1:size(s1_matrix,2)
    curr_s1 = s1_matrix(:,i);
    prod_1 = curr_s1;
    v6_substrates_products(rxn_count,1).s1 = curr_s1;
    v6_substrates_products(rxn_count,1).p1 = prod_1;
    
    rxn_count = rxn_count + 1;
end

all_reactions(6,1).substrates_products = v6_substrates_products;
all_reactions(6,1).substrate_names = {'D'};
all_reactions(6,1).product_names = {'F'};



% All Metabolites

A_matrix = get_isotopomer_matrix(3);
A_chars = get_isotopomer_chars(A_matrix,'A');

B_matrix = get_isotopomer_matrix(3);
B_chars = get_isotopomer_chars(B_matrix, 'B');

C_matrix = get_isotopomer_matrix(2);
C_chars = get_isotopomer_chars(C_matrix, 'C');

D_matrix = get_isotopomer_matrix(3);
D_chars = get_isotopomer_chars(D_matrix, 'D');

E_matrix = get_isotopomer_matrix(1);
E_chars = get_isotopomer_chars(E_matrix, 'E')';

F_matrix = get_isotopomer_matrix(3);
F_chars = get_isotopomer_chars(F_matrix, 'F');

metabolite_vec = [A_chars;B_chars;C_chars;D_chars;E_chars';F_chars];


% All Reactions

reaction_vec = {};
count = 1;
for i=1:6
    SP = all_reactions(i).substrates_products;
    for j=1:length(SP)
        curr_string = strcat('v', num2str(i),'-',num2str(j));
        reaction_vec{count,1} = curr_string;
        count = count + 1;
    end
end
        

S_matrix = zeros(length(metabolite_vec),length(reaction_vec));

curr_reaction = 1;
for i =1:length(all_reactions)
    sp = all_reactions(i).substrates_products;
    num_substrates = length(all_reactions(i).substrate_names);
    num_products = length(all_reactions(i).product_names);
    
    for j=1:length(sp)
        for k=1:num_substrates
            eval(strcat('curr_substrate = all_reactions(i).substrates_products(j).s',num2str(k),';'));
            substrate_string = strcat(all_reactions(i).substrate_names{k}, '-');
            for l = 1:length(curr_substrate)
                substrate_string = strcat(substrate_string, num2str(curr_substrate(l)));
            end
            met_idx = find(ismember(metabolite_vec, substrate_string));
            S_matrix(met_idx,curr_reaction) = S_matrix(met_idx,curr_reaction)-1;
        end
        
        for k=1:num_products
            eval(strcat('curr_product = all_reactions(i).substrates_products(j).p',num2str(k),';'));
            product_string = strcat(all_reactions(i).product_names{k}, '-');
            for l = 1:length(curr_product)
                product_string = strcat(product_string, num2str(curr_product(l)));
            end
            met_idx = find(ismember(metabolite_vec, product_string));
            S_matrix(met_idx,curr_reaction) = S_matrix(met_idx,curr_reaction)+1;
        end
            
        curr_reaction = curr_reaction + 1;
    end
end
            
            
steady_state_vec = [];
for i=1:length(metabolite_vec)
    curr_met = metabolite_vec{i};
    if contains(curr_met,'B')==1 || contains(curr_met,'C') ==1 || contains(curr_met,'D')==1
        steady_state_vec(i,1) = 1;
    else
        steady_state_vec(i,1) = 0;
    end
end

S_matrix_SS = S_matrix(find(steady_state_vec==1),:);

    
    