function dzdt = EX3_NLode_2025(t,z,M,m,k,c,cd1,cd2,cd3,A,T)
 % computes the derivative of the states z
    F = 0;
    P = A*sin(2*pi*t/T); % COMPLETE
    if t>T/2
        P = 0; % COMPLETE
    end

    dzdt(1,1) = z(3); % COMPLETE
    dzdt(2,1) = z(4); % COMPLETE
    dzdt(3,1) = 1/m*((z(4)-z(3))*c+(z(2)-z(1))*k+F); % COMPLETE
    dzdt(4,1) = 1/M*((z(3)-z(4))*c+(z(1)-z(2))*k-cd1*z(4)-cd2*z(4)^2-cd3*z(4)^3+P); % COMPLETE

end