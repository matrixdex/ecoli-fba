function [time, x, names] = Toy_Model_Maciek_Dynamics_execute_101122_730pm(maciek_model)

csObj = getconfigset(maciek_model);
csObj.StopTime = 20000;
csObj.TimeUnits = 'hour';
set(csObj, 'SolverType', 'ode15s');

[time, x, names] = sbiosimulate(maciek_model, csObj);

