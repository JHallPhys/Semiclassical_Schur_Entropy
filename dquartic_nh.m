function pdr=dquartic_nh(t,y,F,omega,hbar_eff,gamma) 

q=y(1);
p=y(2);

% pdr=[-p;hbar_eff^2*q^3 - F/hbar_eff*cos(omega*t)-gamma*p;-gamma*p.^2*y(3)];
pdr=[-p;hbar_eff^2*q^3 - F/hbar_eff*cos(omega*t)-gamma*p;-gamma*p.^2];

end