function [T,S,E,D,nD,U,R0,Age_structure, FinalState] = ODEs_without_vaccine(M_from_to_H, M_from_to_O, alpha, gamma, sigma, d, tau, HHQ, N0, MaxTime, InitialState) %#codegen
%% function to numerically solve system of ODEs. No vaccination included.

% Inputs:
% M_from_to_H, M_from_to_O - Mixing matrices for household & non-household (other) settings 
% alpha - 1/latent_period (rate of leaving latent class)
% gamma - Recovery rate
% sigma - Vector, age-dependent susceptibility 
% d - Vector, age-dependent symptomatic probability 
% tau - Vector, age-dependent relative infectiousness of an asymptomatic
% HHQ - Probability of quarantining if symptomatic
% N0 - Vector of population size per age group
% MaxTime - End timepoint the ODEs are to be evaluated over.
% InitialState - Vector of initial conditions, entry for each compartment in ODE system

% Outputs:
% T - Vector of timepoints
% S,E,D,nD,U - Arrays with row per timepoint & column per age group. 
%               Number of susceptible, latent, symptomatic (detectable),
%               new symtomatic/detectable, asymptomatic (undetectable)
% R0 - Basic reproduction number
% Age_structure - Eigenvector of next generation matrix (used to compute R0)
% FinalState - Values at last timepoint

% Declare number of latent stages
number_E_states=3;

% Number of age groups
L=length(N0);

% Shape sigma and symptomatic probability (d) based on if age-independent
% (length of 1) or age-dependent.
if length(sigma)==1
    sigma=sigma+0*N0;
else
    sigma=reshape(sigma,L,1);
end
if length(d)==1
    d=d+0*N0;
else
    d=reshape(d,L,1);
end
if length(tau)==1
    tau=tau+0*N0;
else
    tau=reshape(tau,L,1);
end

% Shape hhq
hhq = zeros(size(N0,1),size(N0,2));
if length(HHQ)==1
    hhq=HHQ*ones(size(N0,1),size(N0,2));
end
if length(HHQ)==4
   hhq=HHQ([1 2 2 2 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4]);
end

% Get overall mixing matrix. Sum of household and non-household settings
M_from_to=M_from_to_H+M_from_to_O;

% Get dominant eigenvalue of next generation matrix
M_from_to_HAT=zeros(L,L);
for f=1:L
    for t=1:L
        M_from_to_HAT(f,t)=M_from_to(f,t)*d(t)*sigma(t)*(1+tau(f)*(1-d(f))/d(f));
    end
end
[V,D]=eig(M_from_to_HAT);
[R0 i]=max(abs(diag(D))); R0=R0/gamma;

% Assign eigenvector associated with dominant eigenvalue to variable
Age_structure=abs(V(:,i));

% If required, set up initial conditions
if InitialState(1)<0
   E0=Age_structure; D0=Age_structure; U0=Age_structure; S0=N0;
   InitialState=zeros(1,441); L=21;
   InitialState(1:L)=S0'; m=L;
   for i=1:number_E_states
       InitialState(m+[1:3*L])=[E0'/number_E_states 0*E0' 0*E0']; m=m+3*L;
   end
   InitialState(m+[1:4*L])=[D0' zeros(1,length(E0)*2) U0'];
end

% The main iteration 
options = odeset('RelTol', 1e-8);
[t, pop]=ode45(@Diff_3_4,[0:1:MaxTime],[InitialState],options,[L number_E_states reshape(M_from_to_H,1,[]) reshape(M_from_to_O,1,[]) sigma' d' tau' alpha gamma N0' hhq']);

% Assign ODE solutions to output variables
T=t; S=pop(:,1:L); EF=pop(:,L+[1:L]); ES1=pop(:,2*L+[1:L]); ES2=pop(:,3*L+[1:L]); 
DF=pop(:,4*L+[1:L]); DS1=pop(:,5*L+[1:L]); DS2=pop(:,6*L+[1:L]);
UF=pop(:,7*L+[1:L]); US=pop(:,8*L+[1:L]); 
EQ=pop(:,9*L+[1:L]);  DQF=pop(:,10*L+[1:L]); DQS=pop(:,11*L+[1:L]); UQ=pop(:,12*L+[1:L]); 

m=number_E_states;

S=pop(:,1:L);  EF=0*S; ES1=0*S; ES2=0*S; EQ=0*S;
for i=1:m
    EF=EF+pop(:,L+[1:L]+3*(i-1)*L); ES1=ES1+pop(:,2*L+[1:L]+3*(i-1)*L); ES2=ES2+pop(:,3*L+[1:L]+3*(i-1)*L);
    EQ=EQ+pop(:,6*L+[1:L]+3*m*L+(i-1)*L);
end
DF=pop(:,1*L+[1:L]+3*m*L); DS1=pop(:,2*L+[1:L]+3*m*L); DS2=pop(:,3*L+[1:L]+3*m*L); 
UF=pop(:,4*L+[1:L]+3*m*L); US=pop(:,5*L+[1:L]+3*m*L); 
DQF=pop(:,6*L+[1:L]+4*m*L); DQS=pop(:,7*L+[1:L]+4*m*L); UQ=pop(:,8*L+[1:L]+4*m*L);

E=EF+ES1+ES2+EQ;  D=DF+DS1+DS2+DQF+DQS; U=UF+US+UQ;

% Find all new detectable / symptomatic infections
nD=zeros(size(pop,1),length(d));
for a=1:length(d)
    nD(:,a)=d(a)*alpha*m*(pop(:,L+a+3*(m-1)*L)+pop(:,2*L+a+3*(m-1)*L)+pop(:,3*L+a+3*(m-1)*L)+pop(:,6*L+a+3*m*L+(m-1)*L));
end

% Assign quantities at final timepoint to output variable
FinalState=pop(end,:);



% Calculates the differential rates used in the integration.
function dPop=Diff_3_4(t, pop, parameter)

% Deaggregate parameter into separate variables 
L=parameter(1);  % Number of age groups
m=parameter(2);  % Number of latent stages

% M_from_to_H, M_from_to_O - Mixing matrices for household & non-household (other) settings 
% sigma - Vector, age-dependent susceptibility 
% d - Vector, age-dependent symptomatic probability 
% tau - Vector, age-dependent relative infectiousness of an asymptomatic
% alpha - 1/latent_period (rate of leaving latent class)
% gamma - Recovery rate
% N - Vector of population size per age group
% hhq - Probability of quarantining if symptomatic
s=3; l=L*L;
l=L*L;  M_from_toH=reshape(parameter(s+[1:l]-1),L,L);    s=s+l;
l=L*L;  M_from_toO=reshape(parameter(s+[1:l]-1),L,L);    s=s+l;
l=L;    sigma = parameter(s+[1:l]-1);           s=s+l;
l=L;    d = parameter(s+[1:l]-1);               s=s+l;
l=L;    tau = parameter(s+[1:l]-1);             s=s+l;

l=1;    alpha=parameter(s+[1:l]-1);                 s=s+l;
l=1;    gamma=parameter(s+[1:l]-1);                 s=s+l; 
l=L;    N=parameter(s+[1:l]-1);                     s=s+l;
l=L;    hhq=parameter(s+[1:l]-1);                   s=s+l;

% S,E,D,U - Susceptible, latent, symptomatic (detectable), asymptomatic (undetectable)
% Q - Quarantined
% DS1 are infected by a detected and cannot cause lock-down.
% DS2 are infected by undected and can. 
% In the write-up these are D^SD and D^SU
S=pop(1:L)'; 
EF=zeros(m,L); ES1=zeros(m,L); ES2=zeros(m,L); EQ=zeros(m,L);
for i=1:m
    EF(i,:)=pop(L+[1:L]+3*(i-1)*L)';  ES1(i,:)=pop(2*L+[1:L]+3*(i-1)*L)'; ES2(i,:)=pop(3*L+[1:L]+3*(i-1)*L)';
end
DF=pop(1*L+[1:L]+3*m*L)'; DS1=pop(2*L+[1:L]+3*m*L)'; DS2=pop(3*L+[1:L]+3*m*L)'; 
UF=pop(4*L+[1:L]+3*m*L)'; US=pop(5*L+[1:L]+3*m*L)'; 
for i=1:m
    EQ(i,:)=pop(6*L+[1:L]+3*m*L+(i-1)*L)';
end
DQF=pop(6*L+[1:L]+4*m*L)'; DQS=pop(7*L+[1:L]+4*m*L)'; UQ=pop(8*L+[1:L]+4*m*L)'; 

% if hhq is in place - put individals straight into Q classes.

% Initialise output variable
dPop=zeros(length(pop),1);

% Effective number of first infecteds (IF) & second infecteds (IS)
IF=DF + tau.*UF;   IS=(DS1+DS2) + tau.*US; 

a=[1:L]; % Indexes with length equal to number of age groups

% Force of infection to cause first infecteds (InfF),
% secondary infecteds caused by a symptomatic (InfS1),
% secondary infecteds caused by an asymptomatic (InfS2),
% and quarantined infecteds (InfSQ)
InfF=(sigma(a).*((IF+ IS)*M_from_toO(:,a))) .* (S(a)./N(a));
InfS1=(sigma(a).*((DF)*M_from_toH(:,a))) .* (S(a)./N(a)) ;
InfS2=(sigma(a).*((tau.*UF)*M_from_toH(:,a))) .* (S(a)./N(a)) ;
InfSQ=(sigma(a).*((DQF)*M_from_toH(:,a))) .* (S(a)./N(a)) ;

% The ODEs giving the rate of change of each epidemiological state
s=0;
%dS
dPop(a+s)= - InfF - InfS1 -InfS2 - InfSQ;       s=s+L;

%dEF   dES1  and dES2 (first ones)
dPop(a+s)=InfF  - m*alpha*EF(1,a);                s=s+L;
dPop(a+s)=InfS1 - m*alpha*ES1(1,a);               s=s+L;
dPop(a+s)=InfS2 - m*alpha*ES2(1,a);               s=s+L;

%dEF  dES1  & dES2 (subsequent ones)
for i=2:m
    dPop(a+s)=m*alpha*EF(i-1,a) - m*alpha*EF(i,a);      s=s+L;
    dPop(a+s)=m*alpha*ES1(i-1,a) - m*alpha*ES1(i,a);    s=s+L;
    dPop(a+s)=m*alpha*ES2(i-1,a) - m*alpha*ES2(i,a);    s=s+L;
end

%dDF  dDS1  and dDS2
dPop(a+s)=d(a) .* (1-hhq(a)) .* alpha .* EF(m,a) * m - gamma*DF(a);       s=s+L;
dPop(a+s)=d(a) .* alpha .* ES1(m,a) * m - gamma*DS1(a);                   s=s+L;
dPop(a+s)=d(a) .* (1-hhq(a)) .* alpha .* ES2(m,a) * m - gamma*DS2(a);     s=s+L;

%dUF and dUS
dPop(a+s)=(1-d(a)) .* alpha .* EF(m,a) * m - gamma*UF(a);                   s=s+L;
dPop(a+s)=(1-d(a)) .* alpha .* (ES1(m,a) + ES2(m,a)) * m - gamma*US(a);     s=s+L;

%dEQ
dPop(a+s) = InfSQ - alpha .* EQ(1,a) * m;       s=s+L;
for i=2:m
    dPop(a+s) = alpha .* EQ(i-1,a) * m - alpha .* EQ(i,a) * m;       s=s+L;
end

%dDQF
dPop(a+s) = d(a) .* hhq(a) .* alpha .* EF(m,a) * m - gamma .* DQF(a);       s=s+L;

%dDQS
dPop(a+s) = d(a) .* hhq(a) .* alpha .* ES2(m,a) * m  +  d(a) .* alpha .* EQ(m,a) * m - gamma*DQS(a);       s=s+L;

%dUQ
dPop(a+s) = (1-d(a)) .* alpha .* EQ(m,a) * m - gamma*UQ(a);       s=s+L;







