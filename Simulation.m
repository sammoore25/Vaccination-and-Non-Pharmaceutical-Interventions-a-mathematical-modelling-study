function [nDC,nHospital, inHospital, nICU, inICU, nDeaths, PD_Lockdown, INF]=Simulation(Impact_of_New_Var,S_New_Var,delay,EffT,Disease_efficacy,WW,vstep,VDEL,vpropm,whov,Region, tau, ALPHA, a, S_Scale, Factor, h_factor, i_factor, d_factor, h_stretch, ...
    i_stretch, Lag, Start_Date,  WALES_FLAG, ComplianceT, Run_stop, Detection, Susceptibility, gamma)
%% Main simulation script, calls in ODEs and determines disease outcomes

S=size(Run_stop);
if S(2)==1
   Run_stop=Run_stop';
   ComplianceT=ComplianceT';
end

load Starting_Data
load Distributions2.mat
load Probabilities.mat
load Regional_PP.mat
Region_PP(Region_PP==0)=10;

%distribitions of outcomes by age
UK_PP=UK_PP'*sum(Region_PP(Region,:))/sum(UK_PP);
UK_from_toH = UK_from_toH .* ( ones(21,1) * (Region_PP(Region,:)./UK_PP));
UK_from_toW = UK_from_toW .* ( ones(21,1) * (Region_PP(Region,:)./UK_PP));
UK_from_toS = UK_from_toS .* ( ones(21,1) * (Region_PP(Region,:)./UK_PP));
UK_from_toO = UK_from_toO .* ( ones(21,1) * (Region_PP(Region,:)./UK_PP));

%starting values
Early_R0=3.8; 
Early_gamma=gamma/1.5;
Z=1.3;


%generate starting state from ODES and initial R
[T2,S,E,D,nD2,U , R0, Da, FinalState]=ODEs_without_vaccine(UK_from_toH*0 , gamma*(UK_from_toH + UK_from_toW + UK_from_toS + UK_from_toO), a, Early_gamma, Susceptibility, Detection, tau, 0, Region_PP(Region,:)' , 2, -ones(1,441));

Susceptibility2=S_Scale*Susceptibility*Early_R0/R0;
Susceptibility3=S_Scale*Susceptibility;

%% EARLY SET-UP
%Run Hot to 12th March (day 71)
FinalState(22:end)=FinalState(22:end)*Factor; 
 %Run ODES
[T,S,E,D,nD,U , R0, Da, FinalState]=ODEs_without_vaccine(gamma*UK_from_toH*Z , gamma*(UK_from_toW + UK_from_toS + UK_from_toO), a, Early_gamma, Susceptibility2, Detection, tau, 0, Region_PP(Region,:)', 71-Start_Date, [FinalState]);
T=T+Start_Date;

%Self Isolation for 4 days
[t,s,e,d,nd,u , R0, Da, FinalState]=ODEs_without_vaccine(gamma*UK_from_toH*Z, gamma*(UK_from_toW + UK_from_toS + UK_from_toO), a, gamma, Susceptibility3, Detection, tau, 0, Region_PP(Region,:)' , 4, [FinalState]);
%record outputs
T=[T; T(end)+t(2:end)]; S=[S; s(2:end,:)]; E=[E; e(2:end,:)]; D=[D; d(2:end,:)]; nD=[nD; nd(2:end,:)]; U=[U; u(2:end,:)];

%Work from Home + Self Isolation for 4 days
[new_UK_from_toH, new_UK_from_toW, new_UK_from_toS, new_UK_from_toO] = Return_New_Matrices(0, 0.3, 0, ComplianceT(:,1), UK_from_toH, UK_from_toW, UK_from_toS, UK_from_toO);
[t,s,e,d,nd,u , R0, Da, FinalState]=ODEs_without_vaccine(gamma*new_UK_from_toH*Z, gamma*(new_UK_from_toW + new_UK_from_toS + new_UK_from_toO), a, gamma, Susceptibility3, Detection, tau, 0, Region_PP(Region,:)', 4, [FinalState]);
T=[T; T(end)+t(2:end)]; S=[S; s(2:end,:)]; E=[E; e(2:end,:)]; D=[D; d(2:end,:)]; nD=[nD; nd(2:end,:)]; U=[U; u(2:end,:)];

%Some x (Social distancing + HHQ + School Closures + Work from Home + Self Isolatio)n for 3 days
[new_UK_from_toH, new_UK_from_toW, new_UK_from_toS, new_UK_from_toO] = Return_New_Matrices(0.8, 0.5, 0.3, ComplianceT(:,1), UK_from_toH, UK_from_toW, UK_from_toS, UK_from_toO);
[t,s,e,d,nd,u , R0, Da, FinalState]=ODEs_without_vaccine(gamma*new_UK_from_toH*Z, gamma*(new_UK_from_toW + new_UK_from_toS + new_UK_from_toO), a, gamma, Susceptibility3, Detection, tau, 0.3*ComplianceT(1,1), Region_PP(Region,:)', 3, [FinalState]);
T=[T; T(end)+t(2:end)]; S=[S; s(2:end,:)]; E=[E; e(2:end,:)]; D=[D; d(2:end,:)]; nD=[nD; nd(2:end,:)]; U=[U; u(2:end,:)];


%% PRE VACCINATION

PD_Lockdown = 0;
for W=1:WW %W corresponds to NPI period
    
    %update NPI compliance
    [aC,new_UK_from_toH, new_UK_from_toW, new_UK_from_toS, new_UK_from_toO] = npi_change(S_New_Var,Impact_of_New_Var,W,T,Region,ComplianceT,Region_PP,PD_Lockdown,Run_stop, UK_from_toH, UK_from_toW, UK_from_toS, UK_from_toO);
    %Run ODES
    [t,s,e,d,nd,u , R01, Da, FinalState]=ODEs_without_vaccine(gamma*new_UK_from_toH*Z, gamma*(new_UK_from_toW + new_UK_from_toS + new_UK_from_toO), a, gamma, Susceptibility3, Detection, tau, 0.8*max(aC), Region_PP(Region,:)' , Run_stop(W)-T(end), [FinalState]);
    if Run_stop(W)-T(end)==1 % if the interval is only 1 day, it outputs more than just two points, so just take the end timepoint
        T=[T; T(end)+t(end)]; S=[S; s(end,:)]; E=[E; e(end,:)]; D=[D; d(end,:)]; nD=[nD; nd(end,:)]; U=[U; u(end,:)];
    else
        T=[T; T(end)+t(2:end)]; S=[S; s(2:end,:)]; E=[E; e(2:end,:)]; D=[D; d(2:end,:)]; nD=[nD; nd(2:end,:)]; U=[U; u(2:end,:)];
    end
end

%% Vaccination 

%VACCINATION setup

%vaccinated susceptibles
Vacd=Region_PP(Region,:).*0;
Vacd2=Region_PP(Region,:).*0;


nV1=vpropm.*Region_PP(Region,:);
V1=vpropm.*s(end,:);

%who is in each vaccination group
vgroups=zeros(length(whov),21);
for i=1:length(whov)
vgroups(i,whov{i}/5+1)=1;
end

%numbers to vaccinate dose 1
gV1=zeros(length(whov),1);
for i=1:length(whov)
gV1(i)=sum(nV1.*vgroups(i,:));
end

%numbers to vaccinate dose 2 
V2=zeros(1000,21);
gV2=zeros(1000,1);


%calculate [disease efficacies 1 dose, disease efficacies 2 dose increase,transmission efficacies 1 dose, transmission efficacies 2 dose increase]
[efficD,eincD,efficT,eincT] = efficacies(gV1,vgroups,Region_PP,Region,EffT,Disease_efficacy(1),Disease_efficacy(2),Disease_efficacy(3),Disease_efficacy(4));


nvdel=1;
%adjust dose numbers proportional to region size
VDEL=VDEL*sum(Region_PP(Region,:))/sum(sum(Region_PP(2:11,:)));
vdel=VDEL(nvdel);

W=WW+1; 

 
 %RUN VACCINATION PERIOD
for i=1:length(VDEL)-1
    
      %find who to vaccinate next
    [gV1,gV2,V1,V2,Vacd,Vacd2,vdel] = Allocate_dose(delay,T(end),i,gV1,gV2,V1,V2,Vacd,Vacd2,vdel,vstep,vgroups,whov);
    
    %update NPI compliance
    [aC,new_UK_from_toH, new_UK_from_toW, new_UK_from_toS, new_UK_from_toO] = npi_change(S_New_Var,Impact_of_New_Var,W,T,Region,ComplianceT,Region_PP,PD_Lockdown,Run_stop,UK_from_toH, UK_from_toW, UK_from_toS, UK_from_toO);
    
    Vacdb4=Vacd;
    Vacd2b4=Vacd2;
    sb4=s(end,:);
    %run ODES
    [t,s,e,d,nd,u , R01, Da, FinalState,Vacd,Vacd2]=ODEs_with_vaccine(Vacd,Vacd2,eincT,efficT,gamma*new_UK_from_toH*Z, gamma*(new_UK_from_toW + new_UK_from_toS + new_UK_from_toO), a, gamma, Susceptibility3, Detection, tau, 0.8*abs(ComplianceT(:,W)), Region_PP(Region,:)' , vstep, [FinalState]);
    %find reduction of symptomatics due to vaccinated disease efficacy
    nd=nd.*min(((Vacdb4-Vacd)./(sb4-s(end,:))).*efficD+((Vacd2b4-Vacd2)./(sb4-s(end,:))).*eincD+(1-(Vacdb4-Vacd+Vacd2b4-Vacd2)./(sb4-s(end,:))),1);
    %record new state
    if vstep==1 % if the interval is only 1 day, it outputs more than just two points, so just take the end timepoint
        T=[T; T(end)+t(end)]; S=[S; s(end,:)]; E=[E; e(end,:)]; D=[D; d(end,:)]; nD=[nD; nd(end,:)]; U=[U; u(end,:)];
    else
        T=[T; T(end)+t(2:end)]; S=[S; s(2:end,:)]; E=[E; e(2:end,:)]; D=[D; d(2:end,:)]; nD=[nD; nd(2:end,:)]; U=[U; u(2:end,:)];
    end
    %reduction in Vaccinated class due to infection
    loss=max(vpropm.*(sb4-s(end,:))-(Vacdb4-Vacd+Vacd2b4-Vacd2),0);
    V1=max(V1-loss,0);
    
    if T(end)>=Run_stop(W) %new NPI period
        W=find(Run_stop>T(end),1,'first');
        if isempty(W)
            break;
        end
        
    end
      vdel=VDEL(i+1)+vdel; 
    
    if sum(gV1)+sum(gV2(i+1:end))<100 %if no more to vaccinate, end period
        break
    end

end


%% POST VACCINATION

%remain npi period
if W<=length(Run_stop)  
     [aC,new_UK_from_toH, new_UK_from_toW, new_UK_from_toS, new_UK_from_toO] = npi_change(S_New_Var,Impact_of_New_Var,W,T,Region,ComplianceT,Region_PP,PD_Lockdown,Run_stop,UK_from_toH, UK_from_toW, UK_from_toS, UK_from_toO);
     Vacdb4=Vacd;
    Vacd2b4=Vacd2;
    sb4=s(end,:);
  %run ODES
    [t,s,e,d,nd,u , R01, Da, FinalState,Vacd,Vacd2]=ODEs_with_vaccine(Vacd,Vacd2,eincT,efficT,gamma*new_UK_from_toH*Z, gamma*(new_UK_from_toW + new_UK_from_toS + new_UK_from_toO), a, gamma, Susceptibility3, Detection, tau, 0.8*abs(ComplianceT(:,W)), Region_PP(Region,:)' ,Run_stop(W)-T(end), [FinalState]);
    %find reduction of symptomatics due to vaccinated disease efficacy
    nd=nd.*min(((Vacdb4-Vacd)./(sb4-s(end,:))).*efficD+((Vacd2b4-Vacd2)./(sb4-s(end,:))).*eincD+(1-(Vacdb4-Vacd+Vacd2b4-Vacd2)./(sb4-s(end,:))),1);
    %record new state
    if Run_stop(W)-T(end)==1 % if the interval is only 1 day, it outputs more than just two points, so just take the end timepoint
        T=[T; T(end)+t(end)]; S=[S; s(end,:)]; E=[E; e(end,:)]; D=[D; d(end,:)]; nD=[nD; nd(end,:)]; U=[U; u(end,:)];
    else
        T=[T; T(end)+t(2:end)]; S=[S; s(2:end,:)]; E=[E; e(2:end,:)]; D=[D; d(2:end,:)]; nD=[nD; nd(2:end,:)]; U=[U; u(2:end,:)];
    end
end

%remaining simulation time
for i= W:length(Run_stop)-1
    W=W+1;
  [aC,new_UK_from_toH, new_UK_from_toW, new_UK_from_toS, new_UK_from_toO] = npi_change(S_New_Var,Impact_of_New_Var,W,T,Region,ComplianceT,Region_PP,PD_Lockdown,Run_stop,UK_from_toH, UK_from_toW, UK_from_toS, UK_from_toO);
   Vacdb4=Vacd;
    Vacd2b4=Vacd2;
    sb4=s(end,:);
  %run ODES
    [t,s,e,d,nd,u , R01, Da, FinalState,Vacd,Vacd2]=ODEs_with_vaccine(Vacd,Vacd2,eincT,efficT,gamma*new_UK_from_toH*Z, gamma*(new_UK_from_toW + new_UK_from_toS + new_UK_from_toO), a, gamma, Susceptibility3, Detection, tau, 0.8*abs(ComplianceT(:,W)), Region_PP(Region,:)' , Run_stop(W)-T(end), [FinalState]);
    %find reduction of symptomatics due to vaccinated disease efficacy
    nd=nd.*min(((Vacdb4-Vacd)./(sb4-s(end,:))).*efficD+((Vacd2b4-Vacd2)./(sb4-s(end,:))).*eincD+(1-(Vacdb4-Vacd+Vacd2b4-Vacd2)./(sb4-s(end,:))),1);
    %record new state 
    if Run_stop(W)-T(end)==1 %  if the interval is only 1 day, it outputs more than just two points, so just take the end timepoint
        T=[T; T(end)+t(end)]; S=[S; s(end,:)]; E=[E; e(end,:)]; D=[D; d(end,:)]; nD=[nD; nd(end,:)]; U=[U; u(end,:)];
    else
        T=[T; T(end)+t(2:end)]; S=[S; s(2:end,:)]; E=[E; e(2:end,:)]; D=[D; d(2:end,:)]; nD=[nD; nd(2:end,:)]; U=[U; u(2:end,:)];
    end
end
  
% Detectable Cases
clear nDC;
nDC(T+Lag,:)=nD;

Assumed_Delay_Reporting_Deaths=zeros(11,1);
Distribution_Hopital_to_Death(29:end)=0;

%calculate outcomes
[nDC, nHospital, inHospital, nICU, inICU, nDeaths, INF]=ODE_to_Observables(T, nDC, E, Region, h_factor,i_factor,d_factor,h_stretch,i_stretch,a,...
    Assumed_Delay_Reporting_Deaths, Distribution_Hosp_Time', Distribution_HospICU_Time', Distribution_ICU_Time', Distribution_Symptoms_to_Hospital, Distribution_Symptoms_to_ICU, Distribution_Hopital_to_Death,...
    Hosp_2_Death, Sympt_2_critcal, Sympt_2_hosp,  WALES_FLAG); 

clear nDC;
nDC(T,:)=nD;

end

%%


