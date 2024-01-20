function dy = three_component_101122_627pm(t,Y)

load params
k1 = params(1);
k2 = params(2);

A = Y(1);
B = Y(2);
C = Y(3);

v(1,1) = k1*A;
v(2,1) = k2*B;

dy = zeros(3,1);
dy(1) = -v(1);
dy(2) = v(1) - v(2);
dy(3) = v(2);



