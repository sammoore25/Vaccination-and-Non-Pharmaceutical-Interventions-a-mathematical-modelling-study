function [efficD,eincD,efficT,eincT] = efficacies(gV1,vgroups,Region_PP,Region,EffT,Phzr_efficacy,AZ_efficacy,PhzrDose1_eff,AZDose1_eff)
%% function to setup efficacies dependent on proportions of each vaccine delivered to each group and interaction between infection protection and disease

% Inputs:
% gV1 - Numbers per age group vaccinated during only Pfizer vacc phase
% vgroups - Binary indicator of the age groups desiganted to receive vacc
%           during Pfizer only phase.
% Region_PP - Population by region and age.
% Region - ID for region currently being analysed
% EffT - Assumed prevention of infection efficacy following two doses
% Phzr_efficacy,AZ_efficacy - Efficacy against symptomatic disease after two doses (Pfizer, AZ) 
% PhzrDose1_eff,AZDose1_eff - Efficacy against symptomatic disease after one dose (Pfizer, AZ)

% Outputs:
% efficD - Scaling to symptomatic disease outcomes after one dose (1 - disease efficacy)
% eincD - Increase in symptomatic disease efficacy with second dose
% efficT - Scaling to force of infection due to prevention of infection efficacy (infection blocking) after one dose (1 - infection_prevention_efficacy)
% eincT - Increase in prevention of infection efficacy with second dose

%Set disease efficacies:

% first 1M doses PHZR alone
phzrdose=1e6*sum(Region_PP(Region,:))/sum(sum(Region_PP(2:11,:))); %proportional PFZR delivery to region size
[~,effthresh]=find(cumsum(gV1)>phzrdose,1,'first'); %which groups recieve these first doses

% then 90%/10% AZR/PHZR mix
mixeff1=(0.1*PhzrDose1_eff+0.9*AZDose1_eff);
mixeff2=(0.1*(Phzr_efficacy-PhzrDose1_eff)+0.9*(AZ_efficacy-AZDose1_eff));

% set age dependent efficacy dependent on proportions PHZR alone and mixed
% doses in each group
efficD=zeros(1,21);eincD=zeros(1,21);
if isempty(effthresh) %if all groups receive PFZR alone
    efficD=efficD+PhzrDose1_eff;
    eincD=eincD+Phzr_efficacy-PhzrDose1_eff;
else
    if effthresh>1 %if some groups recieve just PFZR
        eincD=eincD+min(sum(vgroups(1:effthresh-1,:)),1)*(Phzr_efficacy-PhzrDose1_eff);
        efficD=efficD+min(sum(vgroups(1:effthresh-1,:)),1)*PhzrDose1_eff;
        rem=(phzrdose-gV1(effthresh-1))/gV1(effthresh);
    else%if only some of first group recieve just PFZR
        rem=phzrdose/gV1(effthresh);
    end
    %Group that partially was initially all PFZR, then received mixture
    eincD=eincD+vgroups(effthresh,:)*(rem*(Phzr_efficacy-PhzrDose1_eff)+(1-rem)*mixeff2);
    efficD=efficD+vgroups(effthresh,:)*(rem*PhzrDose1_eff+(1-rem)*mixeff1);

    %Remaining groups recieve mixed dose efficacies
    eincD=eincD+min(sum(vgroups(effthresh+1:end,:)),1)*mixeff2;
    efficD=efficD+min(sum(vgroups(effthresh+1:end,:)),1)*mixeff1;
end

%transmission efficacy 80% after dose 1, then 100% after dose 2
efficT=1-0.8*EffT*ones(1,21); % Scaling to force of infection from first dose (1 - 0.8*inf_prevention_efficacy)
eincT=0.2*EffT*ones(1,21);    % Increase in infection prevention efficacy from second dose

%adjust disease efficacy dependent on transmission efficacy
eincD=eincD+efficD;
efficD=1-min((1-efficD)./(efficT),1);
eincD=1-min((1-eincD)./(efficT-eincT),1);
efficD=1-efficD;
eincD=1-eincD;
end

