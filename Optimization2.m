k = find(steady_state_vec==1);
Ceq = S_matrix(k,:);

R = 40;
con1 = zeros(1,R);
con1([5]) = 1;

con2 = zeros(1,R);
con2(1:8) = 1;

con3 = zeros(1,R);
con3(9:16) = 1;

con4 = zeros(1,R);
con4(17:24) = 1;

con5 = zeros(1,R);
con5(25:32) = 1;

con6 = zeros(1,R);
con6(33:40) = 1;




con_matrix = [con1;con2;con3;con4;con5;con6];
Ceq_new = [Ceq;con_matrix];
deq_new = zeros(size(Ceq_new,1),1);

deq_new(21,1) = 100;
deq_new(22,1) = 100;
deq_new(23,1) = 40;
deq_new(24,1) = 40;
deq_new(25,1) = 60;
deq_new(26,1) = 100;


x0 = zeros(40,0);
x0(5,1) = 100;
x0([9:16]) = 40/length(9:16);
x0([17:24]) = 40/length(17:24);
x0([25:32]) = 60/length([25:32]);
x0([33:40]) = 100/length([33:40]);


lb = zeros(R,1);
ub = 100*ones(R,1);

lb([18,20,22,24]) = 0;
ub([18,20,22,24]) = 0;

options = optimset('MaxIter',2000,'LargeScale','off');
x = lsqlin(Ceq_new,deq_new,[],[],[],[],lb,ub,x0,options);

c = zeros(R,1);
c(30) = 1;

x2 = linprog(-c,[],[],Ceq_new,deq_new,lb,ub);
x3 = linprog(c,[],[],Ceq_new, deq_new,lb,ub);
