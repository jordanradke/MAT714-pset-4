% flux function for godunov method, based on sol'ns of Reimann problem
% at the flux boundary. Specialized to f(u)=u^2 & f'(u) = 0 @ u=0

function [fx] = flux(Ql,Qr)
% parameters: 
% Ql is cell average to left of boundary
% Qr is cell average to right of boundary

% shock speed
s = (Ql + Qr)/2;

u = NaN;

% info propagating to the right
if Ql >= 0 && Qr >= 0
    u = Ql;
% info propagating to the left
elseif Ql < 0 && Qr <= 0
    u = Qr;
% shock: use rankine-hugoniot to determine speed
elseif Ql >= 0 && Qr <= 0
    % shock propatagating w/ positive speed
    if s > 0
        u = Ql;
    else
        u = Qr;
    end
% transonic rarefaction: info to left propagating leftward,
% info to right propagating rightward
elseif Ql < 0 && Qr > 0
    u = 0;
end
    
fx = u.^2/2;

end





