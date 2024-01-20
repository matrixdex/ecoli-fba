xt3=zeros(100,10);
k=[];

for i=1:100
    ws_k_vec = [normrnd(100,100/4);normrnd(0.00477,0.00477/4);normrnd(0.00477,0.00477/4);normrnd(0.00088,0.00088/4);normrnd(8.8*(10^-8),8.8*(10^-8)/4);normrnd(0.0075,0.0075/4)];
    save ws_k_vec ws_k_vec
    [maciek_model] = toy_model_isotopomer_creation_2();
    [time, x, names, F_distribution, M0, M1, M2, M3] = toy_model_isotopomer_execute_100622();
    xt3(i,1:6) = ws_k_vec;
    xt3(i,7)=M0;
    xt3(i,8)=M1;
    xt3(i,9)=M2;
    xt3(i,10)=M3;

    if (abs(M1-0.66) <= 0.05)
        if(abs(M2-0.33) <=0.05)
            k=x;
        end
    end

    
end


