function dxdt = Equation_Solve(t,x,omega,A_0,m,psi,t_mech)

%%%%second order equation has been split into two first order equations,
%%%%To be solved using ODE45.  

if t <= t_mech
	F_t = A_0*sin((pi/t_mech)*t);
else
	F_t = 0;
end

dxdt = [x(2); (-F_t - omega^2*x(1)  - 2*psi*omega*x(2))];

