
function [x,b]=prad_inv(X,I,I0, KB, x0, b0)
%% PRAD_INV: solve prad inverse problem
% 
% Given the proton mapping x' = x + b/KB,
% where b is the line-integrated magnetic field, and KB a constant,
% the detected fluence I is related to the input fluence
% I0(x) by  the Jacobion transformation ,
%
% I(x') = I0(x) / |dx'/dx|
%
% The inverse problem solves for b(x) given I and I0 as data
%
% [x,b] =  PRAD_INV (X, I, I0, KB, <optional> x0, <optional> b0)
%
% inputs:
% X: spatial coordinate of data
% I: measured proton fluence at points X
% I0: undisturbed proton fluence at points X
%     I0 can be a singleton, then I0 is assumed constant over domain
% KB = deflection parameter (constant)
% optional (x0, b0) = starting point from which to start integration;
%   defaults to (x0=x(1),b0=0) if not specified.
%
% outputs:
% x = points in plasma plane where b is determined
% b = line-integrated field from solving inverse problem
%
% The routine throws a warning if x' or x are generated which
% are out of the data domain of X.   Simple extrapolation is used 
% to allow moderately out of domain data in this case, but you are warned.
%
% See also prad_inv_bc, prad_fwd


% Helper function for ODE solver
    function [dbdx] = dbdx_prad(x, b, X0, I, I0, KB)
    
        xp = x + b/KB;
        
        if (xp > max(X0))
            disp(sprintf('db/dx warning: out of bounds at x = %f, b=%f, x'' = %f',  x, b, xp))
        end 
        if (xp < min(X0))
            disp(sprintf('db/dx warning: out of bounds at x = %f, b=%f, x'' = %f', x, b, xp))
        end

        
        I0 = interp1(X0, I0, x, 'linear', 'extrap');
        Ip = interp1(X0, I, xp, 'linear', 'extrap');
        
        dbdx = KB * (I0./Ip-1);
    
    end

% Input processing and cleanups

if nargin < 5
  x0 = X(1);      % integrate from first point if not specified
  b0 = 0;
end

if nargin < 6
    b0 = 0;
end
    

% transpose to required orientation
if size(X,1) == 1
    X = X';
end

if size(I,1) == 1
    I = I';
end

if any( size(I) ~= size(X) )
    error('prad_inv: Mismatch in size(X0) and size(I)')
end


% if I0 is a scalar, assume user means constant fluence of value I0
if (length(I0) == 1)
   I0 = I0 + 0*I;
end

% Requested final domains -- make sure to include starting point x0
Xdomain = unique( sort ([X; x0]) );


% integrate the right domain (x0 <= x) 


x1 = [];
b1 = [];

if (x0 < Xdomain(end))

    X1 = Xdomain( Xdomain >= x0) ;

    [x1, b1] = ode23(@(x,B) dbdx_prad(x,B, X, I, I0,KB), X1, b0, odeset('AbsTol', 1e-6));

    % clean up end of integration region
    r = isnan(b1);
    if any(r)
        disp('Nan''s')
        last = find(r);
        disp('nan at x0 = %f' \ {x0(last)} ),
    end
    x1 = x1(~r);
    b1 = b1(~r);
end


% interate the left domain X(1) <= x <= x0)
x2 = [];
b2= [];
if (x0 > Xdomain(1))

    X2 = Xdomain( Xdomain <= x0) ;

    [x2, b2] = ode23(@(x,B) dbdx_prad(x,B, X, I, I0,KB), X2(end:-1:1), b0, odeset('AbsTol', 1e-6));

    % clean up end of integration region
    r = isnan(b2);
    if any(r)
        disp('Nan''s')
        last = find(r);
        disp('nan at x0 = %f' \ {x0(last)} ),
    end
    x2 = x2(~r);
    b2 = b2(~r);
end

%glue domains together    
x = [x2(end:-1:2); x1];
b = [b2(end:-1:2); b1];  
    
end
