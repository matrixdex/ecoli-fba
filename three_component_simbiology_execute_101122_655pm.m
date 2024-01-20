function [time, x, names] = three_component_simbiology_execute_101122_655pm(model_obj)

csObj = getconfigset(model_obj);
csObj.StopTime = 60;
csObj.TimeUnits = 'hour';
set(csObj, 'SolverType', 'ode15s');

[time, x, names] = sbiosimulate(model_obj, csObj);