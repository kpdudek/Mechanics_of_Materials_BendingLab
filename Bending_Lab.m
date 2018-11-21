function Bending_Lab
%%% - - - Data Definitions - - - %%%
%%% Cantilever beam
cant_mass = [20,40,60,80]; %g

cant_al = [.131,.252,.396,.500]; %in.
cant_al_len = 9.0; %in.

cant_br = [.100,.190,.265,.342]; %in.
cant_br_len = 9.0; %in.

cant_ss = [.046,.091,.128,.164]; %in.
cant_ss_len = 10.5; %in.

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

%%% Three point bending 9" deflection along length :: 200g
tpb_9in_deflection200_mass = 200; %g
tpb_9in_deflection_lengths = [1,2,3,4,5,6,7,8,9]; %in.

tpb_9in_deflection200_al = [.077,.145,.212,.281,.330,.384,.424,.446,.466]; %in.
tpb_9in_deflection200_br = [.051,.090,.119,.159,.194,.240,.260,.289,.301]; %in.
tpb_9in_deflection200_ss = [.035,.055,.074,.093,.110,.128,.139,.146,.151]; %in.

%%% - - -  Analysis - - - %%%
% Plot the force versus deflection for all materials and all tests
plot_ALL_force_vs_deflection(cant,tpb_6in,tpb_7in,tpb_9in)

% Plot and apply linear fit to each test
plot_force_vs_deflection(cant,'Cantilever')
plot_force_vs_deflection(tpb_6in,'Simple Support 6"')
plot_force_vs_deflection(tpb_7in,'Simple Support 7.5"')
plot_force_vs_deflection(tpb_9in,'Simple Support 9"')


function plot_ALL_force_vs_deflection(cant,tpb_6,tpb_7,tpb_9)
force_cant = ((cant(1,:) .* .001).*9.8) .* 0.224808943; % Convert g to N, then N to lbf
force_tpb6 = ((tpb_6(1,:) .* .001).*9.8) .* 0.224808943; % Convert g to N, then N to lbf
force_tpb7 = ((tpb_7(1,:) .* .001).*9.8) .* 0.224808943; % Convert g to N, then N to lbf
force_tpb9 = ((tpb_9(1,:) .* .001).*9.8) .* 0.224808943; % Convert g to N, then N to lbf

figure
plot(cant(2,:),force_cant,'ok',cant(3,:),force_cant,'ob',cant(4,:),force_cant,'or')
hold on
plot(tpb_6(2,:),force_tpb6,'+k',tpb_6(3,:),force_tpb6,'+b',tpb_6(4,:),force_tpb6,'+r')
hold on
plot(tpb_7(2,:),force_tpb7,'*k',tpb_7(3,:),force_tpb7,'*b',tpb_7(4,:),force_tpb7,'*r')
hold on
plot(tpb_9(2,:),force_tpb9,'sk',tpb_9(3,:),force_tpb9,'sb',tpb_9(4,:),force_tpb9,'sr')


label = sprintf('Force vs Deflection for all Tests');
title(label)

xlabel('Deflection (in.')
ylabel('Force (lbf)')

function plot_force_vs_deflection(data,NAME)
force_data = ((data(1,:) .* .001).*9.8) .* 0.224808943; % Convert g to N, then N to lbf
label = sprintf('Force vs Deflection for %s',NAME);

figure('Name',label)
plot(data(2,:),force_data,'ok',data(3,:),force_data,'*b',data(4,:),force_data,'+r')

title(label)

if NAME == "Cantilever"
    xlabel('Tip Deflection (in.)')
else
    xlabel('Midpoint Deflection (in.)')
end

ylabel('Force (lbf)')




