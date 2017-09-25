function [x,term]=Broyden(f,x0,varargin);

%
% function x=Broyden(f,x0,P1,P2,...);
%
% f : function for which we want to find a zero
% x0 : initial condition for x
% P1,... : parameters of the function
%
% x : solution
% Term : Termination status (1->OK, 0-> Failure)
%

eps1 = 1e-8;
eps2 = 1e-8;
x0 = x0(:);
y0 = feval(f,x0,varargin{:});
S = eye(size(x0,1));
err = 1;
while err>0;
    d = -S\y0;
    x = x0+d;
    y = feval(f,x,varargin{:});
    S = S+((y-y0)-S*d)*d'/(d'*d);
    tmp = sqrt((x-x0)'*(x-x0));
    err = tmp-eps1*(1+abs(x));
    ferr = sqrt(y'*y);
    x0 = x;
    y0 = y;
end

if ferr<eps2;
    term = 1;
else
    term = 0;
end
