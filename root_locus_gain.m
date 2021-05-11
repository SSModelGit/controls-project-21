%% Calculate minimum gain
% kp = 1;
% c1 = 3.36;
% c2 = 0.0117;
% m = 0.003;
% kd = 1.2*sqrt(m/c1);
xs = 0.01:0.0001:0.03;
ks_max = zeros(size(xs));
ks_min = ks_max;
kp = 1;
for i = 1:(size(xs,2))
    v_max = 7.5;
    x0 = xs(i);
    m = 0.003;
    R = 2.41;
    k = 1.75e-8;
    i0 = m*9.81/k*x0^4;
    v0 = i0*R;

    c1 = 4*k*i0/x0^5;
    c2 = k/x0^4;
    kd = 1.5*(m/c1)^0.5;

    sys = tf([c2*kd c2], [m 0 -c1]);
    [r,k] = rlocus(sys);
%     rlocus(sys)
    % Calculate minimum gain
    r1_negative = find(r(1,:)<0);
    index_cross = r1_negative(1);
    k_min = k(index_cross);
%     r_cross = r(1,index_cross);
%     plot(k, real(r(1,:)))

    % Calculate maximum gain
    i_max = v_max/R;
    e_max = 0.01;
    dedt_max = 1;

    c_max = e_max*kp + kd*dedt_max;
    k_max = (i_max - i0)/c_max;

%     k_maxed = find(k>k_max);
%     index_max = k_maxed(1);
%     r_maxed = r(1,index_max);
    ks_min(i) = k_min;
    ks_max(i) = k_max;
end
%%
v_max = 7.5;
x0 = 0.015;
m = 0.003;
R = 2.41;
k = 1.75e-8;
i0 = m*9.81/k*x0^4;
v0 = i0*R;

c1 = 4*k*i0/x0^5;
c2 = k/x0^4;
kd = 1.5*(m/c1)^0.5;

sys = tf([c2*kd c2], [m 0 -c1]);
[r,k] = rlocus(sys);

% Calculate minimum gain
r1_negative = find(r(1,:)<0);
index_cross = r1_negative(1);
k_min_15 = k(index_cross);

i_max = v_max/R;
e_max = 0.01;
dedt_max = 1;

c_max = e_max*kp + kd*dedt_max;
k_max = (i_max - i0)/c_max;

close all
figure
rlocus(sys)
figure
hold on
plot(k, r(1,:))
plot(k, 0*k)
xlim([0 k_max])
scatter(k_min_15, r(1, index_cross))
legend('Pole Real Value', '0', 'First Negative r, k = 24.2')
xlabel 'Controler Gain [-]'
ylabel 'Distance from Imaginary Axis [s^-1]'
title 'Minimum Gain for x0 = 0.015 m'
hold off

%%
close all
figure
hold on
plot(xs,ks_min)
plot(xs, ks_max)
xlabel 'x0 (m)'
ylabel 'Controler Gain [-]'
title 'Minimum and Maximum Controler Gain'
legend('Minimum Gain', 'Estimated Maximum Gain', 'location', 'northwest')

hold off