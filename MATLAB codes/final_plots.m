f_s = 260;
w = linspace(0,pi,10000); size_w = size(w,2);
E = exp(-1i*w);
denom = 1-BP_poles(1) - (1+BP_poles(1))*E; 
for j = 2:size(BP_poles,2)
    denom = denom.*(1-BP_poles(j)-(1+BP_poles(j))*E);
end
H_discrete = ((omega_c*B*(1-E.^2)).^8)./denom;
H_mag_dis = abs(H_discrete); H_phase_dis = angle(H_discrete); 
figure(); 
plot(w*(f_s/(2*pi)),H_mag_dis/max(H_mag_dis),'LineWidth',1.2);  hold on; grid on;
plot(w*(f_s/(2*pi)),0.15*ones(1,size_w),'r'); plot(w*(f_s/(2*pi)),0.85*ones(1,size_w),'r');
xlabel('Frequency Axis(in kHz)'); ylabel('Magnitude Response');
title('Plot of the magnitude response for band pass filter');
%do the same for w_small
E = exp(-1i*w_small);
denom2 = 1-BP_poles(1) - (1+BP_poles(1))*E; 
for j = 2:size(BP_poles,2)
    denom2 = denom2.*(1-BP_poles(j)-(1+BP_poles(j))*E);
end
H_discrete = ((omega_c*B*(1-E.^2)).^8)./denom2;
H_mag_dis = abs(H_discrete); H_phase_dis = angle(H_discrete); 
points = [(f_s/(2*pi))*w_small',H_mag_dis'];
for j = 1:4
    plot(points(j,1),points(j,2),'r*','MarkerSize',10)
end


