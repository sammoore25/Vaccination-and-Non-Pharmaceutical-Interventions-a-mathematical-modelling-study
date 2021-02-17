function [efficD,eincD,efficT,eincT] = efficacies(gV1,vgroups,Region_PP,Region,EffT,Phzr_efficacy,AZ_efficacy,PhzrDose1_eff,AZDose1_eff)
%% function to setup efficacies dependent on proportions of each vaccine delivered to each group and interaction between infection protection and disease
% gV1 = total number left to vaccinate with dose 1
% vgroups = vaccination groups by age
% Region_PP = size of region by age
% Region = current region
% EffT = 2dose transmission eficacy
% Phzr_efficacy = 2 dose Phzr efficacy
% AZ_efficacy = 2 dose AZ efficacy
% PhzrDose1_eff = 1 dose Phzr efficacy
% AZDose1_eff = 1 dose AZ efficacy
%
% efficD = disease efficacies by age dose 1
% eincD =  disease efficacies by age dose 2 increase
% efficT = transmission efficacies by age dose 1
% eincT = transmission efficacies by age dose 2 increase

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
    %Remaining roups recieve mixed dose efficacies
    eincD=eincD+vgroups(effthresh,:)*(rem*(Phzr_efficacy-PhzrDose1_eff)+(1-rem)*mixeff2);
    eincD=eincD+min(sum(vgroups(effthresh+1:end,:)),1)*mixeff2;
    efficD=efficD+vgroups(effthresh,:)*(rem*PhzrDose1_eff+(1-rem)*mixeff1);
    efficD=efficD+min(sum(vgroups(effthresh+1:end,:)),1)*mixeff1;
end

%transmission efficacy 80% after dose 1, then 100% after dose 2
efficT=1-0.8*EffT*ones(1,21);
eincT=0.2*EffT*ones(1,21);

%adjust disease efficacy dependent on transmission effcacy (disease
%efficacy
eincD=eincD+efficD;
efficD=1-min((1-efficD)./(efficT),1);
eincD=1-min((1-eincD)./(efficT-eincT),1);
efficD=1-efficD;
eincD=1-eincD;
end

