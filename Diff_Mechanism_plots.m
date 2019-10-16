% @Author: pradeep
% @Date:   2019-10-11 18:34:01
% @Last Modified by:   pradeep
% @Last Modified time: 2019-10-11 18:34:22


Dassault_data = load('Dassault');
ISIS_M3S_data = load('ISIS_M3S');     
ISRO_230_data = load('ISRO_230');
ISRO_298_data = load('ISRO_298');
Planetary_Corp_data = load('PlanetaryCorp_Mk2');


loglog(Dassault_data(:,1),Dassault_data(:,2),'-*')
hold all                                     
loglog(ISIS_M3S_data(:,1),ISIS_M3S_data(:,2),'-*')          
loglog(ISRO_230_data(:,1),ISRO_230_data(:,2),'-*')
loglog(ISRO_298_data(:,1),ISRO_298_data(:,2),'-*')
loglog(Planetary_Corp_data(:,1),Planetary_Corp_data(:,2),'-*')
legend('Dassault','ISIS M3S','ISRO 230','ISRO 298','Planetary Corp');
grid on
xlabel('Frequency(Hz)');
ylabel('Acceleration(m/s^2)');  
