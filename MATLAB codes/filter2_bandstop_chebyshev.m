clear all;
%Filter specifications for the first filter
f_s = 260; %sampling frequency in kHz
m = 22; q_m = floor(0.1*m); r_m = m - 10*q_m;
BL = 25 + 1.9*q_m + 4.1*r_m; BH = BL + 20; 
BP_low = BL - 4; BP_high = BH + 4;
w_small = (2*pi/f_s)*[BP_low, BL, BH, BP_high];
Omega = tan(w_small/2);
del1 = 0.15; del2 = 0.15;
%choosing the parameters for the transformed low pass version
Omega_o = sqrt(Omega(1)*Omega(4));
B = Omega(4)-Omega(1);
Omega_lpf = (B*Omega)./(Omega_o.^2-Omega.^2); % aplying the analog frequency transform
%omega_s = min(abs(Omega_lpf(1)),abs(Omega_lpf(4))); omega_p = 1;
D1 = (1/(1-del1)^2) - 1; D2 = (1/del2)^2 - 1;
epsilon = sqrt(D1); omega_s = min(abs(Omega_lpf(2)),abs(Omega_lpf(3))); omega_p = abs(Omega_lpf(4));
N = ceil(acosh(sqrt(D2/D1))/acosh(omega_s/omega_p));

%Now we will find the poles of the magnitude squared function
B_k = asinh(1/epsilon)/N;
Hsq_poles = zeros(2*N,1);
for k = 1:(2*N)
    A_k = (2*k-1)*(pi/(2*N));
    Hsq_poles(k) = sin(A_k)*sinh(B_k) + 1i*cos(A_k)*cosh(B_k);
end
syms z_real z_imaginary; fimplicit((z_real/cosh(B_k)).^2 + (z_imaginary/sinh(B_k)).^2 - 1); hold on;
plot(real(Hsq_poles(5:8)),imag(Hsq_poles(5:8)),'g*',real(Hsq_poles(1:4)),imag(Hsq_poles(1:4)),'bo'); grid on; 
legend('Ellipse','Poles of analog LPF', 'Mirror images of poles of analog LPF');
title('Poles of magnitude squared function of analog Chebyshev LPF'); ylabel('Imaginary Axis'); xlabel('Real Axis');
xlim([-1.5,1.5]); ylim([-1.5,1.5]);
%plotting the frequecny response from the transfer function
my_poles = Hsq_poles(5:8);
s_cord = linspace(-4,4,100000); s = 1i*s_cord; 
den = s - my_poles(1);
for j = 2:size(my_poles,1)
    den = den.*(s - my_poles(j));
end
den = sqrt(1+D1)*den; H = prod(my_poles)./den;
H_mag = abs(H); H_phase = angle(H); %plot(s,H_phase); 
figure(); 
plot(s_cord,H_mag,'blue','LineWidth',2);  hold on; grid on;
plot(s_cord,H_phase,'red','LineWidth',1); legend('Magnitude response of the analog low pass','Phase response of the analog low pass')
xlabel('Analog Frequency Axis'); ylabel('Phase and Magnitude response');
title('Plot of the frequency response for Chebyshev analog low pass filter');

%Finding poles of the band stop transfer function
BS_poles = calc_poles_Chebyshev(Omega_o,B,my_poles); 
figure();
plot(BS_poles,'g*'); legend('Poles of analog BSF'); grid on;
title('Poles of analog band stop filter'); ylabel('Imaginary Axis'); xlabel('Real Axis');
xlim([-1.5,1.5]); ylim([-1.5,1.5]);
%Plotting the frequency response of the band-stop filter
den = s-BS_poles(1); 
for j = 2:size(BS_poles,1)
    den = den.*(s-BS_poles(j));
end
H = ((s.^2+Omega_o^2).^4)./den; H = H/sqrt(1+D1); 
H_mag = abs(H); H_phase = angle(H); %plot(s,H_phase);  
figure(); 
plot(s_cord,H_mag,'green','LineWidth',2);  hold on; grid on;
plot(s_cord,H_phase,'red','LineWidth',1); legend('Magnitude response of the analog band stop','Phase response of the analog band stop')
xlabel('Analog Frequency Axis'); ylabel('Phase and Magnitude response');
title('Plot of the frequency response for analog band stop filter');

%Finding poles of the discrete time system
zeros(1) = (1-1i*Omega_o)/(1+1i*Omega_o); zeros(2) = (1+1i*Omega_o)/(1-1i*Omega_o); 
z_poles = (1+BS_poles)./(1-BS_poles);   
figure();
plot(z_poles,'r*'); hold on; plot(zeros,'bo'); grid on;
syms z_real z_imaginary; fimplicit(z_real.^2 + z_imaginary.^2 - 1);
legend('Poles of Discrete Time BSF Filter','Zeros of Discrete Time BSF(Order = 4 each)','Unit Circle')
title('Poles and Zeros of Discrete BSF filter'); ylabel('Imaginary Axis'); xlabel('Real Axis');
xlim([-2 2]);
ylim([-2 2]);

%plotting the frequency response of the discrete time filter
w_cord = linspace(-pi,pi,10000);
E = exp(-1i*w_cord);
denom = 1-BS_poles(1) - (1+BS_poles(1))*E; 
for j = 2:size(BS_poles,1)
    denom = denom.*(1-BS_poles(j)-(1+BS_poles(j))*E);
end
numer = (1+1i*Omega_o + (1i*Omega_o-1)*E).^4; numer = numer.*((1-1i*Omega_o - (1i*Omega_o+1)*E).^4);
H_discrete = numer./denom;
H_mag_dis = abs(H_discrete); H_phase_dis = angle(H_discrete);  
figure(); 
plot(w_cord,H_mag_dis/max(H_mag_dis),'green','LineWidth',2);  hold on; grid on;
plot(w_cord,H_phase_dis,'red','LineWidth',1); legend('Magnitude response of the discrete band stop','Phase response of the discrete band stop')
xlabel('Discrete Domain Frequency Axis'); ylabel('Phase and Magnitude Response');
title('Plot of the frequency response for discrete band stop filter');



