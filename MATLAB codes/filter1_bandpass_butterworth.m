clear all;
%Filter specifications for the first filter
f_s = 330; %sampling frequency in kHz
m = 22; q_m = floor(0.1*m); r_m = m - 10*q_m;
BL1 = 25 + 1.7*q_m + 6.1*r_m; BH1 = BL1 + 20; 
BS_low = BL1 - 4; BS_high = BH1 + 4;
w_small = (2*pi/f_s)*[BS_low, BL1, BH1, BS_high];
Omega = tan(w_small/2);
del1 = 0.15; del2 = 0.15;
%choosing the parameters for the transformed low pass version
Omega_o = sqrt(Omega(2)*Omega(3));
B = Omega(3)-Omega(2);
Omega_lpf = (Omega.^2-Omega_o.^2)./(B*Omega); % aplying the analog frequency transform
omega_s = min(abs(Omega_lpf(1)),abs(Omega_lpf(4))); omega_p = 1; 
%Now we perform the calculations for the butterworth design
D1 = (1/(1-del1)^2) - 1; D2 = (1/del2)^2 - 1;
N = ceil((log(D2/D1))/(2*log(omega_s/omega_p)));
omega_c_upper = omega_s/(D2^(1/(2*N))); omega_c_lower = omega_p/(D1^(1/(2*N))); 
omega_c = (omega_c_upper+omega_c_lower)/2;
%roots of the denominator ..........
k = 0:(2*N-1); power = (2*k+1)*(pi/(2*N));
s_k = (1i)*omega_c*exp((1i)*power);
x_cord = real(s_k); y_cord = imag(s_k);
x_cord1 = zeros(0,1); y_cord1 = zeros(0,1); x_cord2 = zeros(0,1); y_cord2 = zeros(0,1); 
for j = 1:size(x_cord,2)
    if x_cord(j) < 0
        x_cord1(end+1) = x_cord(j); y_cord1(end+1) = y_cord(j);
    else
        x_cord2(end+1) = x_cord(j); y_cord2(end+1) = y_cord(j);
    end
end
plot(x_cord1,y_cord1,'g*',x_cord2,y_cord2,'bo'); legend('Poles of analog LPF filter', 'Mirror images of poles of analog LPF filter');
title('Poles of magnitude plot of analog LPF filter'); ylabel('Imaginary Axis'); xlabel('Real Axis');

s = linspace(-pi,pi,1000);
s = 1i*s; 
my_poles = x_cord1 + 1i*y_cord1;
s1 = s - my_poles(1); s2 = s - my_poles(2); s3 = s - my_poles(3); s4 = s - my_poles(4);
s5 = s - my_poles(5); s6 = s - my_poles(6); s7 = s - my_poles(7); s8 = s - my_poles(8);
denominator = s1.*s2.*s3.*s4.*s5.*s6.*s7.*s8;
s = -1i*s;
H = (omega_c^8)./denominator; 
H_mag = abs(H); H_phase = angle(H); %plot(s,H_phase); 
figure(); 
plot(s,H_mag,'blue','LineWidth',3);  hold on; grid on;
plot(s,H_phase,'red','LineWidth',1.5); legend('Magnitude response of the analog low pass','Phase response of the analog low pass')
xlabel('Analog Frequency Axis'); ylabel('Phase and Magnitude response');
title('Plot of the frequency response for analog low pass filter');

%calculating the frequency response for the band-pss filter
% s_cord = linspace(-pi,pi,1000);
% s = 1i*s_cord; s = (s.^2 + Omega_o^2)./(B*s);
% my_poles = x_cord1 + 1i*y_cord1;
% s1 = s - my_poles(1); s2 = s - my_poles(2); s3 = s - my_poles(3); s4 = s - my_poles(4);
% s5 = s - my_poles(5); s6 = s - my_poles(6); s7 = s - my_poles(7); s8 = s - my_poles(8);
% denominator = s1.*s2.*s3.*s4.*s5.*s6.*s7.*s8;
% H = (omega_c^8)./denominator; 
% H_mag = abs(H); H_phase = angle(H); %plot(s,H_phase); 
% figure();
% plot(s_cord,H_mag,'blue','LineWidth',3);  hold on;
% plot(s_cord,H_phase,'red','LineWidth',1.5); legend('Magnitude response of the analog low pass','Phase response of the analog low-pass')

BP_poles = calc_poles(Omega_o,B,my_poles); 
figure();
plot(BP_poles,'g*'); legend('Poles of analog BPF filter')
title('Poles of analog BPF filter'); ylabel('Imaginary Axis'); xlabel('Real Axis');

s_cord = linspace(-pi,pi,1000);
s = 1i*s_cord; 
den = s-BP_poles(1); 
for j = 2:size(BP_poles,2)
    den = den.*(s-BP_poles(j));
end
H = ((omega_c*B*s).^8)./den;
H_mag = abs(H); H_phase = angle(H); %plot(s,H_phase); 
figure(); 
plot(s_cord,H_mag,'green','LineWidth',3);  hold on; grid on;
plot(s_cord,H_phase,'red','LineWidth',1.5); legend('Magnitude response of the analog band pass','Phase response of the analog band pass')
xlabel('Analog Frequency Axis'); ylabel('Phase and Magnitude response');
title('Plot of the frequency response for analog band pass filter');

%Converting into discete filter by using bilateral transformation
z_poles = (1+BP_poles)./(1-BP_poles);   
figure();
plot(z_poles,'r*'); hold on; plot(1,0,'bo'); plot(-1,0,'bo'); grid on;
legend('Poles of Discrete Time BPF Filter','Zeros of Discrete Time BPF(Order = 8)')
title('Poles and Zeros of Discrete BPF filter'); ylabel('Imaginary Axis'); xlabel('Real Axis');
syms z_real z_imaginary; fimplicit(z_real.^2 + z_imaginary.^2 - 1);
xlim([-2 2]);
ylim([-2 2]);

%Plotting the frequency response of the discrete filter that we have designed
w = linspace(-pi,pi,10000);
E = exp(-1i*w);
denom = 1-BP_poles(1) - (1+BP_poles(1))*E; 
for j = 2:size(BP_poles,2)
    denom = denom.*(1-BP_poles(j)-(1+BP_poles(j))*E);
end
H_discrete = ((omega_c*B*(1-E.^2)).^8)./denom;
H_mag_dis = abs(H_discrete); H_phase_dis = angle(H_discrete);  
figure(); 
plot(w,H_mag_dis/max(H_mag_dis),'green','LineWidth',3);  hold on; grid on;
plot(w,H_phase_dis,'red','LineWidth',1.5); legend('Magnitude response of the discrete band pass','Phase response of the discrete band pass')
xlabel('Discrete Domain Frequency Axis'); ylabel('Phase and Magnitude Response');
title('Plot of the frequency response for discrete band pass filter');




