function [gV1,gV2,V1,V2,Vacd,Vacd2,vdel] = Allocate_dose(vdelay,T,i,gV1,gV2,V1,V2,Vacd,Vacd2,vdel,vstep,vgroups,whov)
%% assign new susceptibles as vaccinated in each timestep
% vdelay = time between doses
% T = time
% i = dose number
% gV1 = total number left to vaccinate with dose 1
% gV2 = total number waiting for dose 2 in this time step
% V1 = number of susceptibles left to vaccinate with dose 1
% V2 = number of susceptibles waiting for dose 2 in this time step
% Vacd = number of susceptibles having received dose 1
% Vacd2 = number of susceptibles having received dose 2
% vdel = number of doses to assign
% vstep = time step
% vgroups = vaccination groups by age
% whov = vaccination group ordering

if gV2(i)>0 %find those ready for 2nd doses first
    propv=min(vdel/gV2(i),1);
    V2(i,:)=min(V2(i,:),Vacd);
    newV2=(round(propv*V2(i,:)));
    Vacd2=Vacd2+newV2;
    Vacd=Vacd-newV2;
    V2(i+1,:)=V2(i+1,:)+round((1-propv)*V2(i,:));
    gV2(i+1)=gV2(i+1)+round((1-propv)*gV2(i));
    vdel=vdel*(1-min(sum(gV2(i))/vdel,1));
end

for j=1:length(whov) %assign remaining vaccine as dose 1 in priority order
    if vdel==0
        break
    end
    if gV1(j)>0
        propv=min(vdel/gV1(j),1);
        if propv==1
            grouptime(j)=T;
        end
        gV2(i+round(vdelay/vstep))=gV2(i+round(vdelay/vstep))+propv*gV1(j);%*0.75/0.85;
        Vacd=Vacd+round(propv*V1.*vgroups(j,:));
        V2(i+round(vdelay/vstep),:)=V2(i+round(vdelay/vstep),:)+round(propv*V1.*vgroups(j,:));%*0.75/0.85;
        gV1(j)=gV1(j)*(1-propv);
        vdel=vdel*(1-min(gV1(j)/vdel,1));
        V1=V1-V1.*vgroups(j,:).*propv;
    end
end
end

