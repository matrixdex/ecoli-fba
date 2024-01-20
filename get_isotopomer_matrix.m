function s_matrix = get_isotopomer_matrix(num_atoms)

for i=1:2^(num_atoms)
    bin_str = dec2bin(i-1);
    if length(bin_str)<num_atoms
        deficiency = num_atoms - length(bin_str);
        new_bin_str = '';
        for j=1:deficiency
            new_bin_str = strcat(new_bin_str, '0');
        end
        
        new_bin_str = strcat(new_bin_str,bin_str);
    else
        new_bin_str = bin_str;
    end
    
    s_vec = zeros(num_atoms,1);
    
    for k=1:length(s_vec)
        s_vec(k,1) = str2num(new_bin_str(k));
    end
    
    s_matrix(:,i) = s_vec;
end

        
    