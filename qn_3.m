# Part (a)
G = 240 * tf([1, 1], [1, 10, 24, 0])
figure(1)
rlocus(G)
Mp = 5;
ts = 0.02;
z = sqrt((log(Mp/100.0)*log(Mp/100.0))/(pi*pi + (log(Mp/100.0)*log(Mp/100.0))));
wn = 4.0/(z*ts);
H_expected = tf([wn*wn], [1, 2*z*wn, wn*wn]);
[num_desired,den_desired] = tfdata(H_expected, 'v');
poles_expected = roots(den_desired)

# Part(b)
curr_zeroes = roots([1, 1]);
curr_poles = roots([1, 10, 24, 0]);

phi = -pi + (-1 * arg(poles_expected(1) - curr_zeroes(1)) + arg(poles_expected(1) - curr_poles(1)) + arg(poles_expected(1) - curr_poles(2)) + arg(poles_expected(1) - curr_poles(3)));
theta = pi - arg(poles_expected(1));
gamma = 0.5*(pi - theta - phi);

zc = wn * (sin(gamma) / sin(theta + phi) )
pc = wn * (sin(gamma + phi) / (sin(theta + gamma + phi)));

# zc = 1 / T and pc = 1 / (alpha * T)
T = 1 / zc
alpha = zc / pc

Gc_lead = tf([1, zc], [1, pc]);
G_new = Gc_lead * G;

# The bode plot is used to find magnitude of transfer function at a point as 
# octave doesn't allow direct substitution to transfer function.
[mag,phase] = bode( G_new, imag(poles_expected(1)) );
Kc = 1/abs(mag*( cos(phase)+1i*sin(phase) ))
G_new = Kc * G_new

# G_for_kv is s*Gc*G, after which s = 0 is substituted 
G_for_kv = G_new * tf([1, 0], [1]);
[num,den] = tfdata(G_for_kv, 'v');
Kv_for_independent_system = num(3)/den(4);

Error = 1 / Kv_for_independent_system
figure(2)
rlocus(G_new)

# Unity feedback closed loop transfer function calculation
C_closed_loop = feedback(G_new, 1);

# Plot step response
t = 0:0.0001:0.2;
figure(3)
step(C_closed_loop, t);