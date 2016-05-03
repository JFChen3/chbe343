function heat_exchanger

%Information needed/solved for:
%Inlet, outlet temperatures
%Flow rates
%HX dimensions
%Thermal data for CO2, IL: viscosity, thermal conductivity, heat capacity
%Conductivity of tube material

%%%DATA, listed in [IL, CO2] (Will be different for each stream integration stage)
mass_flow = [1649.52, 596000/3600]; %kg/s
density = [1115, 29.18]; %kg/m^3
mu = [5.5e-2,3.44e-5]; %Dynamic viscosity, Pa*s
nu = mu./density; %Kinematic viscosity
k = [0.16, 0.056]; %Thermal conductivity, W/(m*K)
Cp = [1800, 1170]; %Heat capacity, J/(kg*K)

L_range = linspace(1,10)';

q_total = 1649.52*1200*91;
q_stage = q_total/4;

dT = q_stage./(mass_flow.*Cp)
eff_NTU_calc(mass_flow, density, mu, nu, k, Cp, L_range, [70+273, 634], q_stage)

end

function [L_exchanger] = eff_NTU_calc(mass_flow, density, mu, nu, k, Cp, L, Ti, q_stage)
%Eff-NTU calculation to get q vs L

%%%Tube dimensions
%Shell
D_large = 4.0;
A_large = pi*(D_large/2)^2;
%Inner tube
Di = 0.03;
Do = 0.0315;
A_small = pi*(Di/2)^2;
A_tube = pi*(Do/2)^2;
N_tubes = floor(A_large/(A_tube*2)); %Number of tubes that fit
% N_tubes = 1000;
fprintf('Number of tubes: %.0f\n', N_tubes)
kt = 385; %Stainless steel tube

%Convert mass flow rates to volumetric flow rates and velocities
volume_flow = mass_flow ./density;
v = volume_flow ./ ([(A_small*N_tubes), (A_large-N_tubes*A_tube)]); %Assume IL on inside

Re = reynolds(nu, v, L); %Reynolds number

Pr = prandtl(Cp, mu, k); %Prandtl number

%%CO2 pressure drop, check to see if need to include
f = 0.184*(Re(2))^-0.2;
dP = (f.*density(2).*(v(2))^2.*L)./(2*D_large);


hi = (k(1)*nusselt(Re(:,1),Pr(1)))./L; %Convection coefficient
ho = (k(2)*nusselt(Re(:,2),Pr(2)))./L; 

qmax = min(mass_flow.*Cp)*(Ti(2)-Ti(1)); %Max possible heat transfer

UA = heat_transfer_coeff(Di, Do, L, kt, hi, ho); %Total heat transfer coeffs
NTU = UA./(min(Cp.*mass_flow)); %Number of transfer units
Cr = min(Cp.*mass_flow)/max(Cp.*mass_flow); %Heat capacity ratio

eps = eff_shell(NTU, Cr);
n=2;
x = (1-eps.*Cr)./(1-eps).^n;
eps_shell = (x - 1)./(x - Cr); %Multi shell pass
%Heat duty
q_shell = N_tubes*eps_shell.*qmax;

plot(L, q_shell)
xlabel('Exchanger Tube Length')
ylabel('Heat Duty')

L_sufficient = L(q_shell > q_stage);
L_exchanger = L_sufficient(1);

end

function [Re] = reynolds(nu, v, L)

Re = L*(v./nu);

end

function [Pr] = prandtl(Cp, mu, k)

Pr = (Cp.*mu)./k;

end

function [Nu] = nusselt(Re, Pr)

Nu = 0.023*Re.^(4/5).*Pr^(0.3); %Check to make sure this is right

end

function [UA] = heat_transfer_coeff(Di, Do, L, kt, hi, ho)

Ai = pi.*Di.*L;
Ao = pi.*Do.*L;

Ri = 1./(hi.*Ai);
Rt = log(Do/Di)./(2.*pi.*kt.*L);
Ro = 1./(ho.*Ao);

UA = (Ri + Rt + Ro).^-1;

end

function [epsilon] = eff_shell(NTU, Cr)

x = sqrt(1+Cr^2);
A = 1+exp(-NTU.*x);
B = 1-exp(-NTU.*x);
epsilon = 2*(1+Cr+x.*(A./B)).^-1;

end
