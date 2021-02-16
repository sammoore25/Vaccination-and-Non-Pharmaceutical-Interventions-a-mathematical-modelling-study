function [efficD,eincD,efficT,eincT] = efficacies(gV1,vgroups,Region_PP,Region,EffT,Phzr_efficacy,AZ_efficacy,PhzrDose1_eff,AZDose1_eff)


%Set disease efficacies: first 1M doses PHZR alone then 90%/10% AZR/PHZR mix
     mixeff1=(0.1*PhzrDose1_eff+0.9*AZDose1_eff);
    mixeff2=(0.1*(Phzr_efficacy-PhzrDose1_eff)+0.9*(AZ_efficacy-AZDose1_eff));
phzrdose=1e6*sum(Region_PP(Region,:))/sum(sum(Region_PP(2:11,:)));
efficT=1-0.8*EffT*ones(1,21);
eincT=0.2*EffT*ones(1,21);
efficD=zeros(1,21);eincD=zeros(1,21);
[~,effthresh]=find(cumsum(gV1)>phzrdose,1,'first');
if isempty(effthresh)
    efficD=efficD+PhzrDose1_eff;
    eincD=eincD+Phzr_efficacy-PhzrDose1_eff;
else
    if effthresh>1
        eincD=eincD+min(sum(vgroups(1:effthresh-1,:)),1)*(Phzr_efficacy-PhzrDose1_eff);
        efficD=efficD+min(sum(vgroups(1:effthresh-1,:)),1)*PhzrDose1_eff;
        rem=(phzrdose-gV1(effthresh-1))/gV1(effthresh);
    else
        rem=phzrdose/gV1(effthresh);
    end
    eincD=eincD+vgroups(effthresh,:)*(rem*(Phzr_efficacy-PhzrDose1_eff)+(1-rem)*mixeff2);
    eincD=eincD+min(sum(vgroups(effthresh+1:end,:)),1)*mixeff2;
     efficD=efficD+vgroups(effthresh,:)*(rem*PhzrDose1_eff+(1-rem)*mixeff1);
    efficD=efficD+min(sum(vgroups(effthresh+1:end,:)),1)*mixeff1;
end
eincD=eincD+efficD;
efficD=1-min((1-efficD)./(efficT),1);
eincD=1-min((1-eincD)./(efficT-eincT),1);
efficD=1-efficD;
eincD=1-eincD;
end

