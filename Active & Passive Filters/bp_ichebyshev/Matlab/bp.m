
f0 = 1700;
f1 = 1400+ 25*1;
f2 = f0^2/f1;
D = (3.5)*((f0^2 - f1^2)/f1);
f3 = (1/2)*(-D + sqrt(D^2 + (f0^2)*4));
f4 = f0^2/f3;
amin = 33 + 25/9;
amax = 0.4 + 1/36;

w0 = 2*pi*f0;
w1 = 2*pi*f1;
w2 = 2*pi*f2;
w3 = 2*pi*f3;
w4 = 2*pi*f4;

%% Multisim frequencies of the input signal
multisim_w0 = w0 - ((w0 - w1)/2);
multisim_w1 = w0 + ((w0 + w1)/3);
multisim_w2 = 0.4 * w3;
multisim_w3 = 2.5 * w4;
multisim_w4 = 3 * w4;

multisim_f0 = multisim_w0 / (2*pi);
multisim_f1 = multisim_w1 / (2*pi);
multisim_f2 = multisim_w2 / (2*pi);
multisim_f3 = multisim_w3 / (2*pi);
multisim_f4 = multisim_w4 / (2*pi);

syms t;
input_signal(t) = cos(multisim_w0*t) + 0.6*cos(multisim_w1*t) + cos(multisim_w2*t) ...
               + 0.8*cos(multisim_w3*t) + 0.4*cos(multisim_w4*t);
figure;
time_sym = 0:0.001:0.15;
ezplot(input_signal(t),time_sym);
Fs = 100000;          % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;            % Length of signal
time = (0:L-1)*T;     % Time vector           


input_signal_fft = cos(multisim_w0*time) + 0.6*cos(multisim_w1*time) + 1.5*cos(multisim_w2*time) ...
                   + 0.7*cos(multisim_w3*time) + 0.4*cos(multisim_w4*time);

Y = fft(input_signal_fft);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% Initialisation 
w0_check = sqrt(w2*w1);
Ws = (w4-w3)/(w2-w1);
bw = w2-w1;

n = acosh(((10^(amin/10) -1)/(10^(amax/10) -1))^(1/2))/acosh(Ws);
n = ceil(n);

epsilon = 1 / ((10^(amin/10) -1)^(1/2));

whp = 1 / (cosh((1/n) * acosh(1/epsilon)));
Whp = whp * bw;
a = (1/n)*(asinh(1/epsilon));

%Gia n = 4 exoume tis parakatw gwnies Butterworth
ps12 = 22.5;
ps34 = 67.5;


%Opote oi poloi tou Chebyshev prokiptoun ws exis
real1 = -sinh(a)*cosd(ps12);
imag1 = cosh(a)*sind(ps12);

real2 = -(sinh(a)*cosd(ps34));
imag2 = cosh(a)*sind(ps34);

W012 = sqrt(real1^2 + imag1^2);
W034 = sqrt(real2^2 + imag2^2);
Q12 = W012/(-2*real1);
Q34 = W034/(-2*real2);
W012_dilda = 1/W012;
W034_dilda = 1/W034;

W012_scaled = W012_dilda * Ws;
W034_scaled = W034_dilda * Ws;

WZ1 = sec(pi/8);
WZ2 = sec((3*pi)/8);

WZ1_dilda = WZ1 *Ws;
WZ2_dilda = WZ2 *Ws;

WZ1_triangle = WZ1_dilda;
WZ2_triangle = WZ2_dilda;

%Poles anodiavatis function
Sigma12_triangle = - W012_scaled / (2*Q12);
Sigma34_triangle = - W034_scaled / (2*Q34);
W12_triangle = sqrt(W012_scaled^2 - Sigma12_triangle^2);
W34_triangle = sqrt(W034_scaled^2 - Sigma34_triangle^2);

%Transformation of complex poles p12
qc12 = w0 / bw;
C1 = Sigma12_triangle^2 + W12_triangle^2;
D1 = -(2*Sigma12_triangle)/qc12;
E1 = 4 + C1/(qc12^2);
G1 = sqrt((E1^2) - 4*(D1^2));
Qcomplex1 = (1/D1)*sqrt((E1+G1)/2);
K1 = -(Sigma12_triangle * Qcomplex1)/qc12;
W1 = K1 + sqrt((K1^2) - 1);
Wcomplex1 = w0 * W1;
Wcomplex2 = w0 / W1;


%Transformation of complex poles p34
qc34 = w0 / bw;
C2 = Sigma34_triangle^2 + W34_triangle^2;
D2 = -(2*Sigma34_triangle)/qc34;
E2 = 4 + C2/(qc34^2);
G2 = sqrt((E2^2) - 4*(D2^2));
Qcomplex2 = (1/D2)*sqrt((E2+G2)/2);
K2 = -(Sigma34_triangle * Qcomplex2)/qc34;
W2 = K2 + sqrt((K2^2) - 1);
Wcomplex3 = w0 * W2;
Wcomplex4 = w0 / W2;

%Transformation of zero wz1
qc = qc12;
K_zero1 = 2 + (WZ1_triangle^2)/(qc^2);
x1 = (K_zero1 + sqrt(K_zero1^2 -4))/2;
wz1_zero1 = w0 * sqrt(x1);
wz2_zero2 = w0 / sqrt(x1);

%Transformation of zero wz2
qc = qc34;
K_zero2 = 2 + (WZ2_triangle^2)/(qc^2);
x2 = (K_zero2 + sqrt(K_zero2^2 -4))/2;
wz3_zero3 = w0 * sqrt(x2);
wz4_zero4 = w0 / sqrt(x2);


%Unit 1 -> LPN Boctor
unit1_w0 = 1;
unit1_wz = wz1_zero1 / Wcomplex1;
unit1_Q = Qcomplex1;
unit1_choose_k1 = 1 / (unit1_wz^2);
unit1_k1 = 0.8;
unit1_R1 = 2 / (unit1_k1*(unit1_wz^2) - 1);
unit1_R2 = 1 / (1 - unit1_k1);
unit1_R3 = (1/2)*((unit1_k1/(unit1_Q^2)) + (unit1_k1*(unit1_wz^2)) - 1);
unit1_R4 = 1 / unit1_k1;
unit1_R5 = 1;
unit1_R6 = 1;
unit1_C1 = unit1_k1 /(2*unit1_Q);
unit1_C2 = 2*unit1_Q;
unit1_k = 1 / (1 + unit1_R3);
unit1_kf = Wcomplex1;
unit1_km = unit1_C1 * (10^8) / unit1_kf;

unit1_scale_R1 = unit1_R1 * unit1_km;
unit1_scale_R2 = unit1_R2 * unit1_km;
unit1_scale_R3 = unit1_R3 * unit1_km;
unit1_scale_R4 = unit1_R4 * unit1_km;
unit1_scale_R5 = unit1_R5 * unit1_km;
unit1_scale_R6 = unit1_R6 * unit1_km;
unit1_scale_C2 = unit1_C2 / (unit1_km*unit1_kf);



%Unit 2 -> HPN
unit2_w0 = 1;
unit2_wz = wz2_zero2 / Wcomplex2;
unit2_Q = Qcomplex1;
unit2_decision_for_Boctor = 1 / (1 - unit2_wz^2);
if (1/unit2_Q <= 1 - unit2_wz^2)
    error('You should use a High Pass Notch.');
end
BoctorHighPass(wz2_zero2, Wcomplex2, unit2_Q);
unit2_k = circuit.H;

%Unit3 -> LPN Boctor
unit3_w0 = 1;
unit3_wz = wz3_zero3 / Wcomplex3;
unit3_Q = Qcomplex2;
unit3_choose_k1 = 1 / (unit3_wz^2);
unit3_k1 = 0.8;
unit3_R1 = 2 / (unit3_k1*(unit3_wz^2) - 1);
unit3_R2 = 1 / (1 - unit3_k1);
unit3_R3 = (1/2)*((unit3_k1/(unit3_Q^2)) + (unit3_k1*(unit3_wz^2)) - 1);
unit3_R4 = 1 / unit3_k1;
unit3_R5 = 1;
unit3_R6 = 1;
unit3_C1 = unit3_k1 /(2*unit3_Q);
unit3_C2 = 2*unit3_Q;
unit3_k = 1 / (1 + unit3_R3);
unit3_kf = Wcomplex3;
unit3_km = unit3_C1 * (10^8) / unit3_kf;

unit3_scale_R1 = unit3_R1 * unit3_km;
unit3_scale_R2 = unit3_R2 * unit3_km;
unit3_scale_R3 = unit3_R3 * unit3_km;
unit3_scale_R4 = unit3_R4 * unit3_km;
unit3_scale_R5 = unit3_R5 * unit3_km;
unit3_scale_R6 = unit3_R6 * unit3_km;
unit3_scale_C2 = unit3_C2 / (unit3_km*unit3_kf);


%Unit 4 -> HPN
unit4_w0 = 1;
unit4_wz = wz4_zero4 / Wcomplex4;
unit4_Q = Qcomplex2;
unit4_decision_for_Boctor = 1 / (1 - unit4_wz^2);
if (1/unit4_Q >= 1 - unit4_wz^2)
    error ('Use Boctor filter instead');
end

unit4_k1 = (Wcomplex4/wz4_zero4)^2 - 1;
unit4_R1 = 1;
unit4_R2 = (unit4_Q^2) * ((unit4_k1+2)^2);
unit4_R3 = 1;
unit4_R4 = (unit4_Q^2)*(unit4_k1+2);
unit4_C = 1 /(unit4_Q*(unit4_k1+2));
unit4_C1 = unit4_k1 * unit4_C;
unit4_k2 = unit4_R4 / (unit4_R4 + 1);
unit4_k = unit4_k2 * ((Wcomplex4/wz4_zero4)^2);

unit4_kf = Wcomplex4;
unit4_km = (unit4_C * (10^8)) / unit4_kf;

unit4_scale_R1 = unit4_R1 * unit4_km;
unit4_scale_R2 = unit4_R2 * unit4_km;
unit4_scale_R3 = unit4_R3 * unit4_km;
unit4_scale_R4 = unit4_R4 * unit4_km;
unit4_scale_C1 = unit4_C1 / (unit4_km*unit4_kf);

unit4_k_low_freq = unit4_k2;


k_total = unit1_k*unit2_k*unit3_k*unit4_k;
% a_gain = (10^0.5) / k_total;
%% Transfer functions
T1 = tf([1 0 wz1_zero1^2],[1 Wcomplex1/Qcomplex1 Wcomplex1^2]);
T1 = unit1_k * T1;

T2 = tf([1 0 wz2_zero2^2],[1 Wcomplex2/Qcomplex1 Wcomplex2^2]);
T2 = unit2_k * T2;

T3 = tf([1 0 wz3_zero3^2],[1 Wcomplex3/Qcomplex2 Wcomplex3^2]);
T3 = unit3_k * T3;

T4 = tf([1 0 wz4_zero4^2],[1 Wcomplex4/Qcomplex2 Wcomplex4^2]);
T4 = unit4_k * T4;

T1_norm_w0 = (unit1_k * (wz1_zero1^2 - w0^2)) / sqrt((Wcomplex1^2 - w0^2)^2 + (Wcomplex1*w0/Qcomplex1)^2);
T2_norm_w0 = (unit2_k * (wz2_zero2^2 - w0^2)) / sqrt((Wcomplex2^2 - w0^2)^2 + (Wcomplex2*w0/Qcomplex1)^2);
T3_norm_w0 = (unit3_k * (wz3_zero3^2 - w0^2)) / sqrt((Wcomplex3^2 - w0^2)^2 + (Wcomplex3*w0/Qcomplex2)^2);
T4_norm_w0 = (unit4_k * (wz4_zero4^2 - w0^2)) / sqrt((Wcomplex4^2 - w0^2)^2 + (Wcomplex4*w0/Qcomplex2)^2);

T_norm_total = abs(T1_norm_w0 * T2_norm_w0 * T3_norm_w0 * T4_norm_w0);
a_gain = (10^0.5) / T_norm_total;
% Bode Diagrams
plot_transfer_function( T1, [10 1000 10000] );
plot_transfer_function( T2, [10 1000 10000] );
plot_transfer_function( T3, [10 1000 10000] );
plot_transfer_function( T4, [10 1000 10000] );
 
T_1 = series(T1,T2);
T_2 = series(T3,T4);

T_total = series(T_1,T_2);
plot_transfer_function( T_total, [10 1000 10000] ); 
ltiview({'bodemag'}, T1, T2, T3, T4, T_total); 

invTLP = inv(T_total);
ltiview({'bodemag'}, invTLP);

T_total = a_gain*T_total;
invTLP = inv(T_total);
ltiview({'bodemag'}, invTLP);

sim = lsim(T_total,input_signal_fft,time);
Y = fft(sim);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

