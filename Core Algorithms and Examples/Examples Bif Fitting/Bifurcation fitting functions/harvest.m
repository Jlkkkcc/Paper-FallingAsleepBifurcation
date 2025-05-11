function dxdt = harvest(t,x,params)

r = params(1);
K = params(2);
h= params(3);
m = params(4);
dxdt = zeros(2,1);
dxdt(1) = r*x(1)*(1-x(1)/K) - x(2)*(x(1)^2/(x(1)^2+h^2));
dxdt(2) = m;

end
