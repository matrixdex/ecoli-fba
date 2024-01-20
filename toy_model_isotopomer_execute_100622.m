function [time, x, names, F_distribution, M0, M1, M2, M3] = toy_model_isotopomer_execute_100622()


maciek_model = toy_model_isotopomer_creation_100622();





csObj = getconfigset(maciek_model);
csObj.StopTime = 20000;
csObj.TimeUnits = 'hour';
% [time,x,names] = sbiosimulate(maciek_model, csObj, dObj1);
[time,x,names] = sbiosimulate(maciek_model, csObj);
set(csObj, 'SolverType', 'ode15s');

final = x(42,[56,64,72,80,112])';
F_distribution = x(end,112:119);
M0 = sum(F_distribution(1))/(sum(F_distribution));
M1 = sum(F_distribution([2,3,5])/sum(F_distribution));
M2 = sum(F_distribution([4,6,7])/sum(F_distribution));
M3 = sum(F_distribution([8])/sum(F_distribution));



