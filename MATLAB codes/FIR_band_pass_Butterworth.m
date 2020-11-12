close all;
del = 0.15; %tolerances
A = -20*log10(del);
extra_samples = 8;
f_s = 330;
f_p1 = 40.6;
f_p2 = 60.6; %all frequencies are in kHz
delta_w = 4*(2*pi/f_s);
w_c1 = (f_p1)*(2*pi/f_s) - delta_w/2;
w_c2 = (f_p2)*(2*pi/f_s) + delta_w/2;

N_min = ceil((A-8)/(4.57*delta_w));
N_min = N_min + extra_samples;
h_fir = ideal_lp_impulse(N_min,w_c2)-ideal_lp_impulse(N_min,w_c1);
stem(-N_min:N_min,h_fir,'filled','LineStyle','-.'); title('Time domain impulse resonse of FIR band-pass filter');
xlabel('Samples'); ylabel('Amplitude');

fvtool(h_fir); %frequency response
%magnitude response
[H,f] = freqz(h_fir,1,100000, f_s);
plot(f,abs(H)); title('Magnitude response of FIR band-pass filter');
ylabel('Magnitude'); xlabel('Frequency(in kHz)');
grid on;

