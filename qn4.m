# Part (a)
G = 100 * tf([1, 1], [1, 2, 2]) * tf([1, 0.01], [1, 0.02, 0.0101]) * tf([1], [1, 10]);
figure(1)
margin(G)

% gain margin, phase margin, phase crossover, gain crossover % 
[gamma, phi, w_gamma, w_m] = margin(G) 

Mp = 10;
ts = 2;
z = sqrt((log(Mp/100.0)*log(Mp/100.0))/(pi*pi + (log(Mp/100.0)*log(Mp/100.0))));
wn = 4.0/(z*ts);
H_expected = tf([wn*wn], [1, 2*z*wn, wn*wn]);
[num_desired,den_desired] = tfdata(H_expected, 'v');
poles_expected = roots(den_desired)
G_expected = tf([wn*wn], [1, 2*z*wn, 0]);
[gamma_expected, phi_expected, w_gamma, w_phi] = margin(G_expected)

phi_m = (phi_expected - phi) * pi/180;
alpha = (1 - sin(phi_m))/(1 + sin(phi_m));
T = 1 / (sqrt(alpha) * w_m);
Gc = tf([T, 1], [alpha*T, 1])
figure(2)
margin(G*Gc)

# 58 degree is added here as the phase margin gotten till now is 20 degree while
# required is 58 degree.
phi_m = (phi_expected - phi + 58) * pi/180;
alpha = (1 - sin(phi_m))/(1 + sin(phi_m))
gain = 20*log10(1/sqrt(alpha))

# The magnitude of Gc at the peak is 20log(1/sqrt(alpha)), so we find the 
# frequency at which the original open loop transfer function gives the negative 
# of this magnitude, so after multiplying them, we would get the gain crossover
# frequency (Wm) at this frequency.
w_m = 6.515;
T = 1 / (sqrt(alpha) * w_m);

Gc = tf([T, 1], [alpha*T, 1])
figure(3)
margin(G*Gc)

feedback(G*Gc, 1);
t = 0:0.01:10;
figure(4)
step(10*feedback(G*Gc, 1), t);
