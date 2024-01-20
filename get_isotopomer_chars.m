function isotopomer_charvec = get_isotopomer_chars(isotopomer_matrix, species_name)

count = 1;
for i = 1:size(isotopomer_matrix,2)
    isotopomer_char = '';
    for j=1:length(isotopomer_matrix(:,i))
        isotopomer_char(j) = num2str(isotopomer_matrix(j,i));
    end
    isotopomer_char = strcat(species_name,'-', isotopomer_char);
    
    isotopomer_charvec{count,1} = isotopomer_char;
    count = count + 1;
end
