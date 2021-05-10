%% Calculate minimum gain
kp = 1;
kd = 0.0251;
c1 = 5.886;
c2 = 0.1064;

m = 0.003;

sys = tf([c2*kd c2*kp], [m 0 -c1]);
[r,k] = rlocus(sys);

r1_negative = find(r(1,:)<=0);
index_cross = r1_negative(1);
k_min = k(index_cross)
%plot(k, real(r(1,:)))

%% Calculate maximum gain
i0 = 0.196;
i_max = 10/2.41;
e_max = 0.01;
dedt_max = 1;

c_max = e_max*kp + kd*dedt_max;
k_max = (i_max - i0)/c_max

