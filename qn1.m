
# Part (a)
num = [1];
den = [1,1,0];
Mp = 17;
ts = 3;
kv = 20;
G = tf(num,den);
z = sqrt((log(Mp/100.0)*log(Mp/100.0))/(pi*pi + (log(Mp/100.0)*log(Mp/100.0))));
wn = 4.0/(z*ts);
H_expected = tf([wn*wn], [1, 2*z*wn, wn*wn]);
[num_desired,den_desired] = tfdata(H_expected, 'v');
poles_expected = roots(den_desired)
curr_poles = roots(den);

# Part (b)
figure(1)
rlocus(G)

# Part (c)
phi = -pi + (arg(poles_expected(1) - curr_poles(1)) + arg(poles_expected(1) - curr_poles(2)))
theta = pi - arg(poles_expected(1));
gamma = 0.5*(pi - theta - phi);

zc = wn * (sin(gamma) / sin(theta + phi) );
pc = wn * (sin(gamma + phi) / (sin(theta + gamma + phi)));

alpha = zc / pc;
beta = 1 / alpha;

Gc_lead = tf([1, zc], [1, pc]);
Gc_lag = tf([1, 0.1], [1, 0.1/beta]); % Take tau_2 as 10 %
G_new = Gc_lead * Gc_lag * G;

# The bode plot is used to find magnitude of transfer function at a point as 
# octave doesn't allow direct substitution to transfer function.
[mag,phase] = bode( G_new, imag(poles_expected(1)) );
K = 1/abs(mag*( cos(phase)+1i*sin(phase) ));

# New G is plotted with K
G_new = K * G_new
figure(2)
rlocus(G_new)

# G_for_kv is s*Gc*G, after which s = 0 is substituted 
G_for_kv = K * G * Gc_lead * tf([1, 0], [1]);
[num,den] = tfdata(G_for_kv, 'v');

# From G_for_kv, when s = 0 is substituted, only the following coefficients remain. 
Kv_for_system = num(2)/den(3) * (1/alpha)


# Part (d)
beta_ind = (kv) / (Kv_for_system * alpha);
Gc_lag_ind = tf([1, 0.1], [1, 0.1/beta_ind]);
G_new = K * Gc_lead * Gc_lag_ind * G
figure(3)
rlocus(G_new)

G_for_kv = G_new * tf([1, 0], [1]);
[num,den] = tfdata(G_for_kv, 'v');
Kv_for_independent_system = num(3)/den(4)

# Plot step response
feedback(G*Gc, 1);
t = 0:0.01:10;
figure(4)
step(feedback(G*Gc, 1), t);
