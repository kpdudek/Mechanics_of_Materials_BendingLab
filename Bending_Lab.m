function Bending_Lab
close all
clear
clc

%%% - - - Data Definitions - - - %%%
%%% Cantilever beam
cant_mass = [20,40,60,80]; %g
cant_lengths = [9.0,9.0,10.5]; %in. al,br,ss

cant_al = [.131,.252,.396,.500]; %in.

cant_br = [.100,.190,.265,.342]; %in.

cant_ss = [.046,.091,.128,.164]; %in.

cant = [cant_mass;cant_al;cant_br;cant_ss];

%%% Three point bending 6"
tpb_6in_mass = [50,100,150,200]; %g

tpb_6in_al = [.034,.066,.101,.131]; %in.
tpb_6in_br = [.034,.063,.087,.107]; %in.
tpb_6in_ss = [.014,.035,.044,.047]; %in.

tpb_6in = [tpb_6in_mass;tpb_6in_al;tpb_6in_br;tpb_6in_ss];

%%% Three point bending 7.5"
tpb_75in_mass = [50,100,150,200]; %g

tpb_75in_al = [.066,.130,.195,.252];
tpb_75in_br = [.05,.112,.114,.186];
tpb_75in_ss = [.023,.044,.065,.087];

tpb_7in = [tpb_75in_mass;tpb_75in_al;tpb_75in_br;tpb_75in_ss];

%%% Three point bending 9"
tpb_9in_mass = [50,100,150,200]; %g

tpb_9in_al = [.102,.211,.323,.440]; %in.
tpb_9in_br = [.084,.151,.209,.287]; %in.
tpb_9in_ss = [.037,.075,.113,.148]; %in.

tpb_9in = [tpb_9in_mass;tpb_9in_al;tpb_9in_br;tpb_9in_ss];

%%% Three point bending 9" deflection along length :: 100g
tpb_9in_deflection100_mass = 100; %g
tpb_9in_deflection100_lengths = [1,2,3,4,5,6,7,8,9]; %in.

tpb_9in_deflection100_al = [.025,.067,.102,.145,.174,.198,.215,.223,.230]; %in.
tpb_9in_deflection100_br = [.025,.038,.051,.061,.084,.104,.116,.124,.140]; %in.
tpb_9in_deflection100_ss = [.020,.028,.039,.047,.055,.063,.069,.071,.077]; %in.

tpb_deflection_100 = [tpb_9in_deflection100_al;tpb_9in_deflection100_br;tpb_9in_deflection100_ss];

%%% Three point bending 9" deflection along length :: 200g
tpb_9in_deflection200_mass = 200; %g
tpb_9in_deflection_lengths = [1,2,3,4,5,6,7,8,9]; %in.

tpb_9in_deflection200_al = [.077,.145,.212,.281,.330,.384,.424,.446,.466]; %in.
tpb_9in_deflection200_br = [.051,.090,.119,.159,.194,.240,.260,.289,.301]; %in.
tpb_9in_deflection200_ss = [.035,.055,.074,.093,.110,.128,.139,.146,.151]; %in.

tpb_deflection_200 = [tpb_9in_deflection200_al;tpb_9in_deflection200_br;tpb_9in_deflection200_ss];

%%% - - -  Analysis - - - %%%
% Plot the force versus deflection for all materials and all tests
plot_ALL_force_vs_deflection(cant,tpb_6in,tpb_7in,tpb_9in)

% Plot and apply linear fit to each test
fprintf('\nLinear Equations and 0 Deflection...')
k_cant = plot_force_vs_deflection(cant,'Cantilever');
k_tpb_6 = plot_force_vs_deflection(tpb_6in,'Simple Support 6"');
k_tpb_7 = plot_force_vs_deflection(tpb_7in,'Simple Support 7.5"');
k_tpb_9 = plot_force_vs_deflection(tpb_9in,'Simple Support 9"');

% Compute the Euler-Bernoulli value
fprintf('\nK_Bending Values...\n')
euler_cant = bending_stiffness(k_cant,cant_lengths,'Cantilever');
euler_tpb_6 = bending_stiffness(k_tpb_6,6.0,'Simple Support 6"');
euler_tpb_7 = bending_stiffness(k_tpb_7,7.5,'Simple Support 7.5"');
euler_tpb_9 = bending_stiffness(k_tpb_9,9.0,'Simple Support 9"');

% Calculate the average values
average_euler_values(euler_cant,euler_tpb_6,euler_tpb_7,euler_tpb_9)

% Calculate beam deflection
beam_deflection(tpb_deflection_200,tpb_deflection_100)
    

function plot_ALL_force_vs_deflection(cant,tpb_6,tpb_7,tpb_9)
force_cant = ((cant(1,:) .* .001).*9.8) .* 0.224808943; % Convert g to N, then N to lbf
force_tpb6 = ((tpb_6(1,:) .* .001).*9.8) .* 0.224808943; % Convert g to N, then N to lbf
force_tpb7 = ((tpb_7(1,:) .* .001).*9.8) .* 0.224808943; % Convert g to N, then N to lbf
force_tpb9 = ((tpb_9(1,:) .* .001).*9.8) .* 0.224808943; % Convert g to N, then N to lbf

%Plot raw data
figure('Name','All Data')
plot(cant(2,:),force_cant,'ok',cant(3,:),force_cant,'ob',cant(4,:),force_cant,'or')
legend('Cantilever Al','Cantilever Brass','Cantilever Steel')
hold on
plot(tpb_6(2,:),force_tpb6,'+k',tpb_6(3,:),force_tpb6,'+b',tpb_6(4,:),force_tpb6,'+r')
legend('Simple Support 6" Al','Simple Support 6" Brass','Simple Support 6" Steel')
hold on
plot(tpb_7(2,:),force_tpb7,'*k',tpb_7(3,:),force_tpb7,'*b',tpb_7(4,:),force_tpb7,'*r')
legend('Simple Support 7.5" Al','Simple Support 7.5" Brass','Simple Support 7.5" Steel')
hold on
plot(tpb_9(2,:),force_tpb9,'sk',tpb_9(3,:),force_tpb9,'sb',tpb_9(4,:),force_tpb9,'sr')

legend('Cantilever Al','Cantilever Brass','Cantilever Steel',...
    'Simple Support 6" Al','Simple Support 6" Brass','Simple Support 6" Steel',...
    'Simple Support 7.5" Al','Simple Support 7.5" Brass','Simple Support 7.5" Steel',...
    'Simple Support 9" Al','Simple Support 9" Brass','Simple Support 9" Steel','Location','southeast')

% Add axes labels
label = sprintf('Force vs Deflection for all Tests');
title(label)
xlabel('Deflection (in.)')
ylabel('Force (lbf)')

function k = plot_force_vs_deflection(data,NAME)
force = ((data(1,:) .* .001).*9.8) .* 0.224808943; % Convert g to N, then N to lbf
al = data(2,:);
br = data(3,:);
ss = data(4,:);
label = sprintf('Force vs Deflection for %s',NAME);

% Plot the raw data
figure('Name',label)
plot(al,force,'ok',br,force,'*b',ss,force,'+r')
hold on

fprintf('\nLinear Fits for %s\n',NAME)
% linear fit for aluminum
lin_al = fitlm(al,force,'linear');
b_al = table2array(lin_al.Coefficients(1,'Estimate'));
m_al = table2array(lin_al.Coefficients(2,'Estimate'));
al_fit = al .* m_al + b_al;

al_eq = sprintf('Aluminum: P = %.3f*delta + %.3f',m_al,b_al);
fprintf('%s\n',al_eq)

delta_offset_al = -b_al/m_al;
fprintf('Aluminum: Delta offset = %.3f\n',delta_offset_al)

% linear fit for brass
lin_br = fitlm(br,force,'linear');
b_br = table2array(lin_br.Coefficients(1,'Estimate'));
m_br = table2array(lin_br.Coefficients(2,'Estimate'));
br_fit = br .* m_br + b_br;

br_eq = sprintf('Brass: P = %.3f*delta + %.3f',m_br,b_br);
fprintf('%s\n',br_eq)

delta_offset_br = -b_br/m_br;
fprintf('Brass: Delta offset = %.3f\n',delta_offset_br)

% linear fit for steel
lin_ss = fitlm(ss,force,'linear');
b_ss = table2array(lin_ss.Coefficients(1,'Estimate'));
m_ss = table2array(lin_ss.Coefficients(2,'Estimate'));
ss_fit = ss .* m_ss + b_ss;

ss_eq = sprintf('Steel: P = %.3f*delta + %.3f',m_ss,b_ss);
fprintf('%s\n',ss_eq)

delta_offset_ss = -b_ss/m_ss;
fprintf('Steel: Delta offset = %.3f\n',delta_offset_ss)

% Plot the linear fits
plot(al,al_fit,'k',br,br_fit,'b',ss,ss_fit,'r')

% Add axes labels to plots
title(label)
if NAME == "Cantilever"
    xlabel('Tip Deflection (in.)')
else
    xlabel('Midpoint Deflection (in.)')
end
ylabel('Force (lbf)')
legend(al_eq,br_eq,ss_eq,'Location','southeast')
xlim([0 max([al,br,ss])+.05])
ylim([0 max(force)+.05])

k = [m_al,m_br,m_ss];

function euler = bending_stiffness(k,l,NAME)
k_al = k(1);
k_br = k(2);
k_ss = k(3);

if NAME == "Cantilever"
    k_bending_al = (1/3) * k_al;
    k_bending_br = (1/3) * k_br;
    k_bending_ss = (1/3) * k_ss;
    l_al = l(1);
    l_br = l(2);
    l_ss = l(3);
else
    k_bending_al = (1/48) * k_al;
    k_bending_br = (1/48) * k_br;
    k_bending_ss = (1/48) * k_ss;
    l_al = l;
    l_br = l;
    l_ss = l;
end

fprintf('The K_Bending values for %s are:\n')
fprintf('Aluminum: %.3e\n',k_bending_al)
fprintf('Brass: %.3e\n',k_bending_br)
fprintf('Steel: %.3e\n',k_bending_ss)

euler_al = (1/l_al) * (1/k_bending_al)^(1/3);
euler_br = (1/l_br) * (1/k_bending_br)^(1/3);
euler_ss = (1/l_ss) * (1/k_bending_ss)^(1/3);

euler = [euler_al,euler_br,euler_ss];

function average_euler_values(e1,e2,e3,e4)
e_al = [e1(1),e2(1),e3(1),e4(1)];
e_br = [e1(2),e2(2),e3(2),e4(2)];
e_ss = [e1(3),e2(3),e3(3),e4(3)];

average_al = mean(e_al);
average_br = mean(e_br);
average_ss = mean(e_ss);

fprintf('\nPercent Differences...\n')
percent_differences(e_al,average_al,'ALuminum')
fprintf('\n')
percent_differences(e_br,average_br,'Brass')
fprintf('\n')
percent_differences(e_ss,average_ss,'Steel')

function percent_differences(vals,average,NAME)
names = {'Cantilever','Simple 6"','Simple 7.5"','Simple 9"'};
for i = 1:length(vals)
    diff = (abs(average-vals(i)) / average) * 100;
    fprintf('The percent difference for %s, %s is: %.2f\n',NAME,names{i},diff)
end


function beam_deflection(d_200,d_100)
delta_al = d_200(1,:) - d_100(1,:);
delta_br = d_200(2,:) - d_100(2,:);
delta_ss = d_200(3,:) - d_100(3,:);

x = 1:9;
L = 18;

v = v_x(x,L);

delta_max_al = sum(v'*delta_al)/sum(v'*v);
delta_max_br = sum(v'*delta_br)/sum(v'*v);
delta_max_ss = sum(v'*delta_ss)/sum(v'*v);

theta_max_al = (3*delta_max_al)/L;
theta_max_br = (3*delta_max_br)/L;
theta_max_ss = (3*delta_max_ss)/L;

g = 5/32; %in.

fit_error_al = theta_max_al * g;
fit_error_br = theta_max_br * g;
fit_error_ss = theta_max_ss * g;

% Display fit error of the micrometer
fprintf('\nFit errors...')
errors = [fit_error_al,fit_error_br,fit_error_ss];
disp(errors)

%%% Plot the beam deflections with error bars and theoretical
figure('Name','Deflection AL')
errorbar(x,delta_al,fit_error_al.*ones(1,length(x))); hold on;
plot(x,delta_max_al*v)
xlim([0 10])
xlabel('X position (in.)')
ylabel('\delta (in.)')
title('Beam Deflection for Aluminum')
legend('Measured Beam Deflection','Predicted Beam Deflection','Location','northwest')

figure('Name','Deflection Br')
errorbar(x,delta_br,fit_error_br.*ones(1,length(x))); hold on;
plot(x,delta_max_br*v)
xlim([0 10])
xlabel('X position (in.)')
ylabel('\delta (in.)')
title('Beam Deflection for Brass')
legend('Measured Beam Deflection','Predicted Beam Deflection','Location','northwest')

figure('Name','Deflection SS')
errorbar(x,delta_ss,fit_error_ss.*ones(1,length(x))); hold on;
plot(x,delta_max_ss*v)
xlim([0 10])
xlabel('X position (in.)')
ylabel('\delta (in.)')
title('Beam Deflection for Steel')
legend('Measured Beam Deflection','Predicted Beam Deflection','Location','northwest')

function v = v_x(x,L)
v = (1/L^3) * ((3*L^2 - 4.*x.^2).*x + 4*(abs(x-(L/2)).^3 + (x - (L/2)).^3));






