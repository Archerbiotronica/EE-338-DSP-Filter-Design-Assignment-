function poles = calc_poles(omega,B,p)
    p_size = size(p,2);
    poles = zeros(1,2*p_size);
    for j = 1:p_size
        root = sqrt((p(j)*B)^2 - 4*(omega)^2);
        poles(2*j-1) = ((p(j)*B) + root)/2; poles(2*j) = ((p(j)*B) - root)/2; 
    end
end