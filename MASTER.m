%%MASTER FILE TO PRODUCE FULL EPIDEMIC SIMULATION WITH VACCINATION AND
%%RELAXING RESTRICTIONS
clear all

%LOAD IN STARTING DATA
load parameter_set.mat
load Starting_Data
load Distributions2.mat
load Probabilities.mat
load Regional_PP.mat
Region_PP(Region_PP==0)=100;


REGIONS=2:11; %Select UK regions to simulate
%Regions are {'England','East of England','London','Midlands','North East and Yorkshire','North West','South East','South West','Wales','Scotland','Northern Ireland','United Kingdom'};


% generate run_stops by region
Num_Comp=size(COMPLIANCE,1);
RUN_STOPs = zeros(11,Num_Comp+1);
for Region = REGIONS
    RUN_STOPs(Region,1:Num_Comp) = RUN_STOP;
    Z=(RUN_STOP(end-1)+7):7:(ceil(datenum(2022,1,14)+1-datenum(2020,1,1)));  % Run every 7 days until NPI relaxation complete.
    RUN_STOPs(Region,Num_Comp+[1:length(Z)]-1)=Z;
end
RUN_STOPs=[RUN_STOPs,1600*ones(11,1)];%add simulation end point
maxtime = max(RUN_STOPs(end,:))+30; % needs to be big enough to take account of any lags
RUN_STARTs=[82*ones(11,1) RUN_STOPs(:,1:(end-1))];



%Outputs to record
daily_deaths=zeros(maxtime,1);
daily_hospital_admissions=zeros(maxtime,1);
daily_hospital_occupancy=zeros(maxtime,1);


%Vaccination parameters
target_group_order={80:5:100,75,70,65,60,55,50,10:5:45};
transmission_efficacy=0.6;
Disease_efficacy=[0.94,0.88,0.88,0.7]; %[PHZR 2dose,AZR 2dose, PHZR 1dose, AZR 1dose] 
uptake= [0,0,0,0.75*2/5,ones(1,6)*0.75,ones(1,6)*0.85,ones(1,5)*0.95];
vstep=2; %simulation stepsize
Delivery=[(0.786/3)*1e6,(0.786/3)*1e6,(0.786/3)*1e6,0.32*1e6,1.225e6,1.61e6,2.2515e6,2.75e6,2e6,1.75e6,repmat(2.5e6,1,50)]; %Weekly doses accross all regions
Delivery=repelem(Delivery,7)/7;
if Delivery>1
    Delivery=sum(reshape(Delivery(1:vstep*floor(length(Delivery)/vstep)),vstep,floor(length(Delivery)/vstep)));
end
Vacc_start_date=datenum(2020,12,6)+14-datenum(2020,1,1);
dose_delay=84; %time between doses


%choose relaxation parameters
Relaxation_start=datenum(2021,2,1)+1-datenum(2020,1,1);
Relaxation_end=datenum(2021,12,1)+1-datenum(2020,1,1);
Final_Compliance=0;

% generate compliance timeline by region
Comps=zeros(max(REGIONS),length(RUN_STOPs));
Num_Comp=size(COMPLIANCE,1);

for Region = REGIONS
    Comps(Region,1:Num_Comp)=COMPLIANCE(1:Num_Comp,Region);
    Comps(Region,Num_Comp:size(RUN_STOPs,2))=Comps(Region,Num_Comp);
    m=find(RUN_STARTs(Region,:)>=Relaxation_start & RUN_STARTs(Region,:)<=Relaxation_end);
    Comps(Region,m)=Comps(Region,m).*[1:-(1-Final_Compliance/Comps(Region,m(1)))/(length(m)-1):Final_Compliance/Comps(Region,m(1))];
    Comps(Region,(max(m)+1):end)=Final_Compliance;
end


%% Run simulation


for Region=REGIONS
    
    %Add vaccination start point
    WW=find(RUN_STOPs(Region,:)>Vacc_start_date,1,'first')-1;
    if RUN_STOPs(Region,WW)~=Vacc_start_date
        RUN_STOP2=[RUN_STOPs(Region,1:WW),Vacc_start_date,RUN_STOPs(Region,WW+1:end)];
        Comps2=[Comps(Region,1:WW),Comps(Region,WW),Comps(Region,WW+1:end)];
        WW=WW+1;
    else
        RUN_STOP2=RUN_STOPs(Region,:);
        Comps2=Comps(Region,:);
    end
    
    %simulation
    [~, nHospital, inHospital, ~, ~, nDeaths, ~, ~]=Simulation(Impact_of_New_Var(Region),S_New_Var(Region,:) ,dose_delay,transmission_efficacy,Disease_efficacy,WW,vstep,Delivery,uptake,target_group_order,Region, TAU, ALPHA, INC_P, SCALING(Region), FACTOR(Region), H_FACTOR(Region), I_FACTOR(Region), D_FACTOR(Region),  H_STRETCH(Region), I_STRETCH(Region), LAG(Region), START_DATE(Region)+1, 0, Comps2, RUN_STOP2, Detection, Susceptibility, gamma);
    
    
    %record output
    padding = maxtime-size(nDeaths,1);
    daily_deaths=daily_deaths+sum([(nDeaths);zeros(padding,21)],2);
    daily_hospital_admissions=daily_hospital_admissions+sum([(nHospital);zeros(padding,21)],2);
    daily_hospital_occupancy=daily_hospital_occupancy+sum([(inHospital);zeros(padding,21)],2);
end


%% Plot results
plot_data=daily_deaths;

maxp=max(ceil(plot_data/1000))*1000;
period=[datenum(2020,12,1):datenum(2022,7,1)]+1-datenum(2020,1,1); %plot period

figure;hold on;

X=datenum(2020,1,1)+period;

p1=patch([Relaxation_start,Relaxation_end,Relaxation_end,Relaxation_start]+datenum(2020,1,1),[maxp 0 0 0], [1,1,1]*0.87,'LineStyle','none','HandleVisibility','off');
p2=patch([period(1),Relaxation_start,Relaxation_start,period(1)]+datenum(2020,1,1),[maxp maxp 0 0], [1,1,1]*0.87,'LineStyle','none','HandleVisibility','off');
set(gca, 'Layer', 'top');



plot(X,plot_data(period),'LineWidth',1.5,'color',[0.4940 0.1840 0.5560]);
ylabel('Daily deaths');
datetick('x','mmm-yy','keeplimits');
set(gca, 'fontname', 'times','fontsize',14);
xlim([datenum(2020,1,1)+period(1),datenum(2020,1,1)+period(end)]);
ylim([0,maxp]);
box on
ax = gca;
ax.YGrid = 'on';
ax.GridLineStyle = '-';
