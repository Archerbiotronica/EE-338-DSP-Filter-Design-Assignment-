function poles = calc_poles_Chebyshev(omega,B,p)
    p_size = size(p,1);
    poles = zeros(2*p_size,1);
    for j = 1:p_size
        root = sqrt((B/p(j))^2 - 4*(omega)^2);
        poles(2*j-1) = (B/(p(j)) + root)/2; poles(2*j) = (B/(p(j)) - root)/2; 
    end
end