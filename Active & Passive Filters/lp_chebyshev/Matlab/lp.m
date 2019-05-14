%% Input-Output Signal Spectrum
Fs = 1e7;
t = 0:1/Fs:(10*5e-4);

pulsewidth = 2e-4;
pulseperiods = [0:10]*5e-4;

x = pulstran(t,pulseperiods,@rectpuls,pulsewidth);

plot(t,x)
axis([0 4e-3 -0.5 1.5])

Fs = 20000;          % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1000;            % Length of signal
time = (0:L-1)*T;     % Time vector           

Y = fft(x);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% Initialisation of the filter
fp = 4000;
fs = 1.75*fp;
amin = 24;
amax = 0.5;

wp = 2*pi*fp;
ws = 2*pi*fs;

Ws = ws/wp;

epsilon = (10^(amax/10)-1)^(1/2);
n = acosh(((10^(amin/10) -1)/(10^(amax/10) -1))^(1/2))/acosh(Ws);
n = ceil(n);
a = (1/n)*(asinh(1/epsilon));

WHP = (cosh((1/n) * acosh(1/epsilon)));
whp = WHP * wp;
fhp = whp / (2*pi);
ps12 = 22.5;
ps34 = 67.5;

real1 = -sinh(a)*cosd(ps12);
imag1 = cosh(a)*sind(ps12);

real2 = -(sinh(a)*cosd(ps34));
imag2 = cosh(a)*sind(ps34);

W012 = sqrt(real1^2 + imag1^2);
W034 = sqrt(real2^2 + imag2^2);
Q1 = W012/(-2*real1);
Q2 = W034/(-2*real2);

W012_real_norm = W012 * wp;
W034_real_norm = W034 * wp;


%% Unit1 
unit1_R1 = 1;
unit1_R2 = 1;
unit1_C1 = 1;
unit1_C2 = 1;
unit1_k = 3 - 1/Q1;
unit1_r1 = 1;
unit1_r2 = 2 - 1/Q1;

unit1_kf = W012_real_norm;
unit1_km = (10^6)*unit1_C1/unit1_kf;

unit1_scale_C1 = unit1_C1/(unit1_km*unit1_kf);
unit1_scale_C2 = unit1_C2/(unit1_km*unit1_kf);
unit1_scale_R1 = unit1_R1 * unit1_km;
unit1_scale_R2 = unit1_R2 * unit1_km;
unit1_scale_r1 = unit1_r1 * unit1_km;
unit1_scale_r2 = unit1_r2 * unit1_km;


%% Unit 2
unit2_R1 = 1;
unit2_R2 = 1;
unit2_C1 = 1;
unit2_C2 = 1;
unit2_k = 3 - 1/Q2;
unit2_r1 = 1;
unit2_r2 = 2 - 1/Q2;

unit2_kf = W034_real_norm;
unit2_km = (10^6)*unit2_C1/unit2_kf;

unit2_scale_C1 = unit2_C1/(unit2_km*unit2_kf);
unit2_scale_C2 = unit2_C2/(unit2_km*unit2_kf);
unit2_scale_R1 = unit2_R1 * unit2_km;
unit2_scale_R2 = unit2_R2 * unit2_km;
unit2_scale_r1 = unit2_r1 * unit2_km;
unit2_scale_r2 = unit2_r2 * unit2_km;


%% Bode Diagrams
T1 = tf((unit1_k * (W012_real_norm^2)),[1 W012_real_norm/Q1 W012_real_norm^2]);
T2 = tf((unit2_k * (W034_real_norm^2)),[1 W034_real_norm/Q2 W034_real_norm^2]);


plot_transfer_function( T1, [10 fp fs] );
plot_transfer_function( T2, [10 fp fs] );

k_total = unit1_k * unit2_k;
a_gain = (10^0.5)/k_total;
Z1_gain = unit1_scale_R1/a_gain;
Z2_gain = unit1_scale_R1/(1-a_gain);
T_total = series(T1,T2);
plot_transfer_function( T_total, [10 fp fs] );

ltiview({'bodemag'}, T1, T2, T_total);

invTLP = inv(T_total);
ltiview({'bodemag'}, invTLP);

T_total = a_gain * T_total;
invTLP = inv(T_total);
ltiview({'bodemag'}, invTLP);

%% Input and Output Signal Spectrum
sim = lsim(T_total,x,t);
Y = fft(sim);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
