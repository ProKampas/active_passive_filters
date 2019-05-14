%% Initialisation
fp = 5000;
fs = fp / 2.6;
amin = 24 + 6/9;
amax = 0.5 + 5/36;

wp = 2*pi*fp;
ws = 2*pi*fs;
Ws = wp/ws;
Wp = 1;
n = log((10^(amin/10) -1)/(10^(amax/10) -1))/(2*log(Ws));
n = ceil(n);

W0 = Wp / (10^(amax/10) -1)^(1/(2*n));
w0 = wp / W0;
f0 = w0 / (2*pi);

Q1 = 0.54;
Q2 = 1.31;

%% Multisim frequencies of the input signal
multisim_w0 = 0.2 * ws;
multisim_w1 = 0.7 * ws;
multisim_w2 = 1.6 * wp;
multisim_w3 = 2.4 * wp;
multisim_w4 = 3.5 * wp;

multisim_f0 = multisim_w0 / (2*pi);
multisim_f1 = multisim_w1 / (2*pi);
multisim_f2 = multisim_w2 / (2*pi);
multisim_f3 = multisim_w3 / (2*pi);
multisim_f4 = multisim_w4 / (2*pi);

syms t;
input_signal(t) = cos(multisim_w0*t) + 0.6*cos(multisim_w1*t) + 1.5*cos(multisim_w2*t) ...
               + 0.7*cos(multisim_w3*t) + 0.4*cos(multisim_w4*t);
figure;
time_sym = 0:0.0001:0.04;
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
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
% plot(time, Y);
%Design Strategy 1 for High Pass Sallen-Key
%% Unit1 
unit1_R1 = 1;
unit1_R2 = 1;
unit1_C1 = 1;
unit1_C2 = 1;
unit1_k = 3 - 1/Q1;
unit1_r1 = 1;
unit1_r2 = 2 - 1/Q1;

unit1_kf = w0;
unit1_km = (10^7)*unit1_C1/unit1_kf;

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

unit2_kf = w0;
unit2_km = (10^7)*unit2_C1/unit2_kf;

unit2_scale_C1 = unit2_C1/(unit2_km*unit2_kf);
unit2_scale_C2 = unit2_C2/(unit2_km*unit2_kf);
unit2_scale_R1 = unit2_R1 * unit2_km;
unit2_scale_R2 = unit2_R2 * unit2_km;
unit2_scale_r1 = unit2_r1 * unit2_km;
unit2_scale_r2 = unit2_r2 * unit2_km;

%% Implementation of the traansfer functions
k_total = unit1_k * unit2_k;
a_gain = (10^0.5)/k_total;

T1 = tf([1 0 0],[1 w0/Q1 w0^2]);
T1 = unit1_k * T1;

T2 = tf([1 0 0],[1 w0/Q2 w0^2]);
T2 = unit2_k * T2;

T_total = series(T1, T2);

plot_transfer_function( T1, [10 fs fp] );
plot_transfer_function( T2, [10 fs fp] );
plot_transfer_function( T_total, [10 fs fp] );
ltiview({'bodemag'}, T1, T2, T_total); 

%Before the gain configuration
invT = inv(T_total);
% ltiview({'bodemag'}, invT);

% After the gain configuration
T_total = a_gain*T_total;
invT = inv(T_total);
% ltiview({'bodemag'}, invT);

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