%Cytoscape and adjacency matrix

N = length(metabolite_vec) + length(reaction_vec);
adj_matrix = zeros(N,N);
source_vec = zeros(N,1);
target_vec = zeros(N,1);

fid = fopen('Maciek_Network.txt', 'w+');

fprintf(fid, 'Source\tInteraction\tTarget\n');
count =1;
for i = 1:size(S_matrix,1)
    for j= 1:size(S_matrix,2)
        if S_matrix(i,j)<0
            fprintf(fid,'%s\tInteraction\t%s\n',metabolite_vec{i}, reaction_vec{j});
            source_vec(count,1) = i;
            target_vec(count,1) = length(metabolite_vec) + j;
            count = count + 1;
            adj_matrix(i,length(metabolite_vec)+j) = 1;
        end
        
        if S_matrix(i,j)>0
            fprintf(fid,'%s\tInteraction\t%s\n',reaction_vec{j}, metabolite_vec{i});
            source_vec(count,1) = length(metabolite_vec)+j;
            target_vec(count,1) = i;
            count = count+1;
            adj_matrix(length(metabolite_vec)+j, i) = 1;
        end
    end
end

G = digraph(source_vec,target_vec);
fclose(fid);

fid2 = fopen('Maciek_Node_Attributes.txt', 'w+')
fprintf(fid2, 'Unique_Node\tNode_Type\n');

for i=1:length(metabolite_vec)
    fprintf(fid2, '%s\tMetabolite\n', metabolite_vec{i});
end

for i=1:length(reaction_vec)
    fprintf(fid2, '%s\tReaction\n', reaction_vec{j});
end

fclose(fid2);