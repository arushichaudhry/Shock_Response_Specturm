%
% @Author: pradeep
% @Date:   2019-10-03 18:26:43
% @Last Modified by:   pradeep
% @Last Modified time: 2019-10-03 18:26:43
%

%%%%%%This file is ued to evaluate the shock response specturm of the marman clamp system. 
%%%%%% But, this would work in general to other spring based systems. 

%%%%The system is simplified as a spring-mass system with stiffness K and damping C
%%%%The governing equation is m(d2x/dt2) + c(dx/dt) + kx = cdy/dt + ky
%%%%where, x is the displacement of the mass m, and y is the displacement of the 
%%%%platform the mass m has been connected to. 
%%%%This is a second order differential equation which can be solved using Runge-Kutta method

clear all
close all

%%%ISRO IBL 298 separation system
% v_injection = 1.2;     
% t_mech = 10e-3;
% A_0 = v_injection/t_mech;
% m = 100;
% c = 0.05;

% %%%%ISRO IBL 230 separation system
% v_injection = 1.0;     
% t_mech = 10e-3;
% A_0 = v_injection/t_mech;
% m = 100;
% c = 0.05;

% %%%%ISIS M3S separation system
% v_injection = 1.0;     
% t_mech = 5e-3;
% A_0 = v_injection/t_mech;
% m = 100;
% c = 0.05;

% %%%%Dassault 5SSASAP5S separation system
% v_injection = 3.0;     
% t_mech = 1.5e-3;
% A_0 = v_injection/t_mech;
% m = 100;
% c = 0.05;

%%%%Planetary corp Light band Mk2
v_injection = 3.0;     
t_mech = 0.8e-3;
A_0 = v_injection/t_mech;
m = 100;
c = 0.05;

nat_Freq_vec(1) = 1;
% nat_Freq_upp_vec(1) = nat_Freq_vec(1)*sqrt(2^(1/3));
% nat_Freq_low_vec(1) = nat_Freq_vec(1)/sqrt(2^(1/3));

nat_Freq_highest = 1;
count = 1;

while nat_Freq_highest <= 4000
	count = count+1;
	nat_Freq_vec(count) = nat_Freq_vec(count-1)*2^(1/3);
	% nat_Freq_upp_vec(count) = nat_Freq_vec(count)*sqrt(2^(1/3));
	% nat_Freq_low_vec(count) = nat_Freq_vec(count)/sqrt(2^(1/3));
	nat_Freq_highest = nat_Freq_vec(count);
end

for ii = 1 :length(nat_Freq_vec)
	nat_Freq = nat_Freq_vec(ii);
	% nat_Freq_upper = nat_Freq_upp_vec(ii);
	% nat_Freq_lower = nat_Freq_low_vec(ii);
	
	omega = nat_Freq*2*pi;
	% omega_upp = nat_Freq_upper*2*pi;
	% omega_low = nat_Freq_lower*2*pi;

	K = omega^2*m;
	psi = c/(2*sqrt(K*m));

	% K_upp = omega_upp^2*m;
	% psi_upp = c/(2*sqrt(K_upp*m));

	% K_low = omega_low^2*m;
	% psi_low = c/(2*sqrt(K_low*m));

	options = odeset('RelTol',1e-8,'AbsTol',1e-8,'MaxStep',1e-6);
	[t,x] = ode45(@(t,x)Equation_Solve(t,x,omega,A_0,m,psi,t_mech),[0,0.06],[0;0],options);
	% [t_upp,x_upp] = ode45(@(t,x)Equation_Solve(t,x,omega_upp,A_0,m,psi_upp),[0,0.06],[0;0],options);
	% [t_low,x_low] = ode45(@(t,x)Equation_Solve(t,x,omega_low,A_0,m,psi_low),[0,0.06],[0;0],options);

	for jj = 1:length(t)
		t_eval = t(jj);
		% t_eval_upp = t_upp(jj);
		% t_eval_low = t_low(jj);
		if t_eval <= t_mech
			A_t = A_0*sin((pi/t_mech)*t_eval);
		else
			A_t = 0;
		end
		acc_part(jj,1) = A_t;
	end

	% acc = -x(:,1)*omega^2;
	% max_acc(ii) = max(abs(acc));

	for jj = 1:length(t)
		t_eval = t(jj);
		x_eval = x(jj,:);
		acc1(:,jj) = Equation_Solve(t_eval,x_eval,omega,A_0,m,psi,t_mech); 
	end

	acc1_total = acc1(2,:)' + acc_part(:,1);

	max_acc1(ii) = max(abs(acc1_total));
	
	% acc_upp = -x_upp(:,1)*A_0*omega_upp^2 ;
	% max_acc_upp(ii) = max(abs(acc_upp));

	% acc_low = -x_low(:,1)*A_0*omega_low^2 ;
	% max_acc_low(ii) = max(abs(acc_low));
	
end

figure(1)
loglog(nat_Freq_vec,max_acc1);
% loglog(nat_Freq_vec,max_acc);
xlabel('natural frequency (Hz)')
ylabel('acceleration (g)')
% ylim([10 10^3])
grid on
hold on
% legend('Disp','Equation')

hold on
plot(nat_Freq_vec(10),max_acc1(10),'r*');
str1 = strcat('(',num2str(nat_Freq_vec(10),'%.0f'),',',num2str(max_acc1(10),'%.0f'),')');
text(nat_Freq_vec(10)+0.5,max_acc1(10),str1);

str2 = strcat('(',num2str(nat_Freq_vec(20),'%.0f'),',',num2str(max_acc1(20),'%.0f'),')');
text(nat_Freq_vec(20),max_acc1(20)+20,str2);
plot(nat_Freq_vec(20),max_acc1(20),'r*');

str3 = strcat('(',num2str(nat_Freq_vec(length(nat_Freq_vec)),'%.0f'),',',num2str(max_acc1(length(max_acc1)),'%.0f'),')');
text(nat_Freq_vec(length(nat_Freq_vec)),max_acc1(length(max_acc1))+20,str3);
plot(nat_Freq_vec(length(nat_Freq_vec)),max_acc1(length(max_acc1)),'r*');

app_line_data_x = [ nat_Freq_vec(5) nat_Freq_vec(10) nat_Freq_vec(20) nat_Freq_vec(length(nat_Freq_vec)) ];
app_line_data_y = [ max_acc1(5) max_acc1(10) max_acc1(20) max_acc1(length(max_acc1)) ];
str3 = strcat('(',num2str(nat_Freq_vec(5),'%.0f'),',',num2str(max_acc1(5),'%.0f'),')');
text(nat_Freq_vec(5)+0.2,max_acc1(5),str3);
plot(nat_Freq_vec(5),max_acc1(5),'r*');

plot(app_line_data_x,app_line_data_y,'linewidth',2);

dlmwrite('PlanetaryCorp_Mk2', [app_line_data_x' app_line_data_y'], 'delimiter', '\t');

