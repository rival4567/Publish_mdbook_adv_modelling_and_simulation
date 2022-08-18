
// **********************************************************************
// System order determination for parameter estimation

clear;
clc;
clf;
figure(0)
title("PTn System")
T = 0.05;
t = [0:T:4];
input_step = ones(1, length(t));

function [order, num_r, den_r] = system_order_determination(num, den)
    num_r = roots(num)
    den_r = roots(den)
    n = length(num_r)
    d = length(den_r)
    //format(25)
    if n >= d then order = n
    else order = d 
    end
    for i = 1:n
        for j = 1:d
            real_err = real(num_r(i)) - real(den_r(j))
            imag_err = imag(num_r(i)) - imag(den_r(j))
            
            // Set precision as 0.005
            if (abs(real_err) < 0.005 && abs(imag_err) < 0.005)
                order = order - 1
            end
        end
    end
endfunction

z = poly(0, 'z');
num = -3 -5*z +2*z^2;
den = -10 -6*z +7*z^2;
[order, num_r, den_r] = system_order_determination(num, den)
disp("System model order: ", order, num_r, den_r)

f_z = syslin('d', num, den);
disp(f_z)
yd = flts(input_step, tf2ss(f_z));
plot(t, yd, "r-")
