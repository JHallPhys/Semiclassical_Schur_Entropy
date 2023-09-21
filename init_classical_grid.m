function [q,p,dq,dp,qmesh,pmesh]=init_classical_grid(q_R,p_R,N,hbar_eff)

q=linspace(-q_R,q_R,N); % The initial condition for your q
p=linspace(-p_R,p_R,N); % The range of your p - set to 0 for simplicity can be a range
q=q./hbar_eff;
p=p./hbar_eff;
dq=abs(q(2)-q(1));
dp=abs(p(2)-p(1));
[qmesh,pmesh]=meshgrid(q,p);
end