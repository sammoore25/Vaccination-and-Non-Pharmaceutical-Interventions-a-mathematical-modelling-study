function [T,S,E,D,nD,U,R0,Age_structure, FinalState,Vacd,Vacd2] = ODEs_with_vaccine(Vacd,Vacd2,einc,effic,M_from_to_H, M_from_to_O, alpha, gamma, sigma, d, tau, HHQ, N0, MaxTime, InitialState)
%%


number_E_states=3;

L=length(N0);

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

hhq = zeros(size(N0,1),size(N0,2));
if length(HHQ)==1
    hhq=HHQ*ones(size(N0,1),size(N0,2));
end
if length(HHQ)==4
   hhq=HHQ([1 2 2 2 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4]);
end

M_from_to=M_from_to_H+M_from_to_O;

M_from_to_HAT=zeros(L,L);
for f=1:L
    for t=1:L
        M_from_to_HAT(f,t)=M_from_to(f,t)*d(t)*sigma(t)*(1+tau(f)*(1-d(f))/d(f));
    end
end


[V,D]=eig(M_from_to_HAT);
[R0 i]=max(abs(diag(D))); R0=R0/gamma;
Age_structure=abs(V(:,i));


if isempty(InitialState)
   E0=Age_structure; D0=Age_structure; U0=Age_structure; S0=N0;
   InitialState=[S0'];
   InitialState=[InitialState, S0'*0];
   for i=1:number_E_states
       InitialState=[InitialState E0'/number_E_states 0*E0' 0*E0'];
   end
   InitialState=[InitialState D0' zeros(1,length(E0)*2) U0' zeros(1,length(E0)*(4+number_E_states))];
else
    InitialState=[InitialState(1:L)-Vacd-Vacd2,Vacd,Vacd2,InitialState((L+1):end)];
end



options = odeset('RelTol', 1e-5);


[t, pop]=ode45(@Diff_3_4,[0:1:MaxTime],[InitialState],options,[L number_E_states reshape(M_from_to_H,1,[]) reshape(M_from_to_O,1,[]) sigma' d' tau' alpha gamma N0' hhq',effic,einc]);

Vacd=pop(end,(L+1):(2*L));Vacd2=pop(end,(2*L+1):(3*L));
pop=[pop(:,1:L)+pop(:,(L+1):(2*L))+pop(:,(2*L+1):(3*L)),pop(:,(3*L+1):end)];

T=t; 


m=number_E_states;

S=pop(:,1:L);  EF=0*S; ES1=0*S; ES2=0*S; EQ=0*S;
for i=1:m
    EF=EF+pop(:,L+[1:L]+m*(i-1)*L); ES1=ES1+pop(:,2*L+[1:L]+m*(i-1)*L); ES2=ES2+pop(:,3*L+[1:L]+m*(i-1)*L);
    EQ=EQ+pop(:,6*L+[1:L]+3*m*L+(i-1)*L);
end
DF=pop(:,1*L+[1:L]+3*m*L); DS1=pop(:,2*L+[1:L]+3*m*L); DS2=pop(:,3*L+[1:L]+3*m*L); 
UF=pop(:,4*L+[1:L]+3*m*L); US=pop(:,5*L+[1:L]+3*m*L); 
DQF=pop(:,6*L+[1:L]+4*m*L); DQS=pop(:,7*L+[1:L]+4*m*L); UQ=pop(:,8*L+[1:L]+4*m*L);

E=EF+ES1+ES2+EQ;  D=DF+DS1+DS2+DQF+DQS; U=UF+US+UQ;

% Find all new detectable / symptomatic infections
nD=zeros(size(pop,1),length(d));
for a=1:length(d)
    nD(:,a)=d(a)*m*alpha*(pop(:,L+a+3*(m-1)*L)+pop(:,2*L+a+3*(m-1)*L)+pop(:,3*L+a+3*(m-1)*L)+pop(:,6*L+a+3*m*L+(m-1)*L));
end

FinalState=pop(end,:);



% Calculates the differential rates used in the integration.
function dPop=Diff_3_4(t, pop, parameter)

L=parameter(1);
m=parameter(2);


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
l=L;    effic=parameter(s+[1:l]-1);                   s=s+l;
l=L;    einc=parameter(s);    

%DS1 are infected by a detected and cannot cause lock-down, DS2 are
%infected by undected and can. In the write-up these are D^SD and D^SU
S=pop(1:L)'; V=pop((L+1):(2*L))'; V2=pop((2*L+1):(3*L))'; EF=zeros(m,L); ES1=zeros(m,L); ES2=zeros(m,L); EQ=zeros(m,L);
for i=1:m
    EF(i,:)=pop(3*L+[1:L]+m*(i-1)*L)'; ES1(i,:)=pop(4*L+[1:L]+m*(i-1)*L)'; ES2(i,:)=pop(5*L+[1:L]+m*(i-1)*L)';
end
DF=pop(3*L+[1:L]+3*m*L)'; DS1=pop(4*L+[1:L]+3*m*L)'; DS2=pop(5*L+[1:L]+3*m*L)'; 
UF=pop(6*L+[1:L]+3*m*L)'; US=pop(7*L+[1:L]+3*m*L)'; 
for i=1:m
    EQ(i,:)=pop(8*L+[1:L]+3*m*L+(i-1)*L)';
end
DQF=pop(8*L+[1:L]+4*m*L)'; DQS=pop(9*L+[1:L]+4*m*L)'; UQ=pop(10*L+[1:L]+4*m*L)'; 

% if hhq is in place - put individals straight into Q classes.

dPop=zeros(11*L+4*m*L,1);

IF=DF + tau.*UF;   IS=(DS1+DS2) + tau.*US; 

a=[1:L];
%for a=1:L
    InfF=(sigma(a).*((IF+ IS)*M_from_toO(:,a))) .* (S(a)./N(a));
    InfS1=(sigma(a).*((DF)*M_from_toH(:,a))) .* (S(a)./N(a)) ; 
    InfS2=(sigma(a).*((tau.*UF)*M_from_toH(:,a))) .* (S(a)./N(a)) ; 
    InfSQ=(sigma(a).*((DQF)*M_from_toH(:,a))) .* (S(a)./N(a)) ; 
    
      InfFv=(sigma(a).*((IF+ IS)*M_from_toO(:,a))) .* (V(a)./N(a)).*effic;
    InfS1v=(sigma(a).*((DF)*M_from_toH(:,a))) .* (V(a)./N(a)) .*effic; 
    InfS2v=(sigma(a).*((tau.*UF)*M_from_toH(:,a))) .* (V(a)./N(a)) .*effic; 
    InfSQv=(sigma(a).*((DQF)*M_from_toH(:,a))) .* (V(a)./N(a)) .*effic; 
    
      InfFv2=(sigma(a).*((IF+ IS)*M_from_toO(:,a))) .* (V2(a)./N(a)).*(effic-einc);
    InfS1v2=(sigma(a).*((DF)*M_from_toH(:,a))) .* (V2(a)./N(a)) .*(effic-einc); 
    InfS2v2=(sigma(a).*((tau.*UF)*M_from_toH(:,a))) .* (V2(a)./N(a)) .*(effic-einc); 
    InfSQv2=(sigma(a).*((DQF)*M_from_toH(:,a))) .* (V2(a)./N(a)) .*(effic-einc); 
    
    s=0;
    %dS
    dPop(a+s)= - InfF - InfS1 -InfS2 - InfSQ;       s=s+L;
    
     %dV
    dPop(a+s)= - InfFv - InfS1v -InfS2v - InfSQv;       s=s+L;
     %dV2
    dPop(a+s)= - InfFv2 - InfS1v2 -InfS2v2 - InfSQv2;       s=s+L;

    %dEF   dES1  and dES2 (first ones)
    dPop(a+s)=InfF +InfFv+InfFv2 - m*alpha*EF(1,a);                s=s+L;
    dPop(a+s)=InfS1 +InfS1v+InfS1v2- m*alpha*ES1(1,a);               s=s+L;
    dPop(a+s)=InfS2+InfS2v+InfS2v2 - m*alpha*ES2(1,a);               s=s+L;
    
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
    dPop(a+s) = InfSQ+InfSQv +InfSQv2- alpha .* EQ(1,a) * m;       s=s+L;
    for i=2:m
        dPop(a+s) = alpha .* EQ(i-1,a) * m - alpha .* EQ(i,a) * m;       s=s+L;
    end
    
    %dDQF
    dPop(a+s) = d(a) .* hhq(a) .* alpha .* EF(m,a) * m - gamma .* DQF(a);       s=s+L;

    %dDQS
    dPop(a+s) = d(a) .* hhq(a) .* alpha .* ES2(m,a) * m  +  d(a) .* alpha .* EQ(m,a) * m - gamma*DQS(a);       s=s+L;
    
    %dUQ
    dPop(a+s) = (1-d(a)) .* alpha .* EQ(m,a) * m - gamma*UQ(a);       s=s+L;







