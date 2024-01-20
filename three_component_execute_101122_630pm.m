function [t,Y] = three_component_execute_101122_630pm()

A_0 = 100;
B_0 = 20;
C_0 = 0;


tspan = [0:0.01:60];
y0 = [A_0;B_0;C_0];
[t,Y] = ode45(@three_component_101122_627pm, tspan, y0);

plot(t, Y(:,1))