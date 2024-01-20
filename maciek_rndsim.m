

xt=zeros(100,21);
for i=1:100
    ws_k_vec = [normrnd(100,100/4);normrnd(0.00477,0.00477/4);normrnd(0.00477,0.00477/4);normrnd(0.00088,0.00088/4);normrnd(8.8*(10^-8),8.8*(10^-8)/4);normrnd(0.0075,0.0075/4)];
    save ws_k_vec ws_k_vec
    [maciek_model] = Toy_Model_Maciek_Dynamics_101122_7pm();
    [time, x, names] = Toy_Model_Maciek_Dynamics_execute_101122_730pm(maciek_model);
    xt(i,1:6) = ws_k_vec;
    for j=1:15
        xt(i,j+6)=max(x(:,j));
    end
end
