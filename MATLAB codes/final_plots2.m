f_s = 260;
w_cord = linspace(0,pi,10000); size_w = size(w_cord,2);
E = exp(-1i*w_cord);
denom = 1-BS_poles(1) - (1+BS_poles(1))*E; 
for j = 2:size(BS_poles,1)
    denom = denom.*(1-BS_poles(j)-(1+BS_poles(j))*E);
end
numer = (1+1i*Omega_o + (1i*Omega_o-1)*E).^4; numer = numer.*((1-1i*Omega_o - (1i*Omega_o+1)*E).^4);
numer = numer/sqrt(1+D1);
H_discrete = numer./denom;
H_mag_dis = abs(H_discrete); 
figure(); 
plot((f_s/(2*pi))*w_cord,H_mag_dis,'LineWidth',1.2);  hold on; grid on;
plot(w_cord*(f_s/(2*pi)),0.15*ones(1,size_w),'r'); plot(w_cord*(f_s/(2*pi)),0.85*ones(1,size_w),'r');
xlabel('Frequency Axis(in kHz)'); ylabel('Magnitude Response');
title('Plot of the magnitude response for discrete band-stop filter');
%do the same for w_small
E = exp(-1i*w_small);
denom = 1-BS_poles(1) - (1+BS_poles(1))*E; 
for j = 2:size(BS_poles,1)
    denom = denom.*(1-BS_poles(j)-(1+BS_poles(j))*E);
end
numer = (1+1i*Omega_o + (1i*Omega_o-1)*E).^4; numer = numer.*((1-1i*Omega_o - (1i*Omega_o+1)*E).^4);
numer = numer/sqrt(1+D1);
H_discrete = numer./denom;
H_mag_dis = abs(H_discrete);
points = [(f_s/(2*pi))*w_small',H_mag_dis'];
for j = 1:4
    plot(points(j,1),points(j,2),'r*','MarkerSize',15)
end


