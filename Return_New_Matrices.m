function [new_UK_from_toH, new_UK_from_toW, new_UK_from_toS, new_UK_from_toO] = Return_New_Matrices(SCHOOL_CLOSURES, WORK_CLOSURES, SOCIAL_DISTANCING, Compliance, UK_from_toH, UK_from_toW, UK_from_toS, UK_from_toO)
%% function to amend the mixing matrices

% Inputs:
% SCHOOL_CLOSURES, WORK_CLOSURES - Maximum proportional reduction in
%                                       contacts in school & work settings.
% SOCIAL_DISTANCING - Maximum proportional reduction in contacts in other settings
% Compliance - Adherence to NPIs
% UK_from_toH, UK_from_toW, UK_from_toS, UK_from_toO - 
%   The original mixing matrices/contact arrays.
%   Order is Household, Work, School, Other

% Outputs:
% new_UK_from_toH, new_UK_from_toW, new_UK_from_toS, new_UK_from_toO -
%   The amended mixing matrices/contact arrays.
%   Order is Household, Work, School, Other

% Reshape Compliance input if needed
% Four age groupings: PreSchool, School, Adult, Elderly. 
if length(Compliance)==1
    Compliance=[1 1 1 1]'*Compliance;
end

% Map four age groupings to the five-year age bands
% Preschool: 0-4yrs
% School: 5-19yrs
% Adult: 20-69yrs
% Elderly: 70-100+yrs
if length(Compliance)==4
    Compliance=Compliance([1 2 2 2 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4]);
end

if size(Compliance,2)>1
    Compliance=Compliance';
end

% Raw Uncontrolled Scaling
% Indexes: [Household, Work, School, Other]
ReductionPreSchool=[1.0 1.0 1.0 1.0];
ReductionSchool=[1.0 1.0 1.0 1.0];
ReductionAdult=[1.0 1.0 1.0 1.0];   workers_interacting_with_other=0.3;
ReductionElderly=[1.0 1.0 1.0 1.0];
hhq=0;
txt=[];

% Maximum impact of school closures on each of the four age groupings
if SCHOOL_CLOSURES
    ReductionPreSchool(3)=1-SCHOOL_CLOSURES;  ReductionSchool(3)=1-SCHOOL_CLOSURES; ReductionAdult(3)=1-SCHOOL_CLOSURES; ReductionElderly(3)=1-SCHOOL_CLOSURES;
    ReductionSchool(1)=1+0.25*SCHOOL_CLOSURES; % Home contacts increase for school children
    txt='SC';
end

% Maximum impact of work closures on each of the four age groupings
if WORK_CLOSURES
    ReductionPreSchool(2)=1-WORK_CLOSURES;  ReductionSchool(2)=1-WORK_CLOSURES; ReductionAdult(2)=1-WORK_CLOSURES; ReductionElderly(2)=1-WORK_CLOSURES;
    ReductionAdult(1)=1+0.25*WORK_CLOSURES; % Home contacts increase for adults
    if isempty(txt)  txt='WC'; else txt=[txt '+WC']; end
end

% Maximum impact of social distancing/reduction of contacts in other settings 
% on each of the four age groupings
if SOCIAL_DISTANCING
    ReductionPreSchool(4)=1-SOCIAL_DISTANCING;  ReductionSchool(4)=1-SOCIAL_DISTANCING; ReductionAdult(4)=1-SOCIAL_DISTANCING; ReductionElderly(4)=1-SOCIAL_DISTANCING;
    
    % Home contacts increase in all age age groups
    ReductionPreSchool(1)=1+0.25*SOCIAL_DISTANCING;  ReductionSchool(1)=ReductionSchool(1)+0.25*SOCIAL_DISTANCING; ReductionAdult(1)=ReductionAdult(1)+0.25*SOCIAL_DISTANCING; ReductionElderly(1)=1.0+0.25*SOCIAL_DISTANCING;
    if isempty(txt)  txt='SD'; else txt=[txt '+SD']; end
end


if isempty(txt)  txt='None'; end

% Get maximum reduction scale factors, w.r.t adherence values.
ReductionPreSchool=ReductionPreSchool'*Compliance(1)' + [1 1 1 1]'*(1-Compliance(1));
ReductionSchool=ReductionSchool'*Compliance(2:4)' + [1 1 1 1]'*(1-Compliance(2:4)');
ReductionAdult=ReductionAdult'*Compliance(5:14)' + [1 1 1 1]'*(1-Compliance(5:14)');
ReductionElderly=ReductionElderly'*Compliance(15:21)' + [1 1 1 1]'*(1-Compliance(15:21)');
OtherVector=[ReductionPreSchool(4) ReductionSchool(4,:) ReductionAdult(4,:) ReductionElderly(4,:)];

% Iterate over each five-year age band and populate output mixing matrices

f=1; %pre-school
new_UK_from_toH(f,:)=UK_from_toH(f,:)*ReductionPreSchool(1);
new_UK_from_toW(f,:)=UK_from_toW(f,:)*ReductionPreSchool(2);
new_UK_from_toS(f,:)=UK_from_toS(f,:)*ReductionPreSchool(3);
new_UK_from_toO(f,:)=UK_from_toO(f,:)*ReductionPreSchool(4).*OtherVector;

f=2:4; %school
new_UK_from_toH(f,:)=UK_from_toH(f,:).*(ReductionSchool(1,f+1-f(1))'*ones(1,21));
new_UK_from_toW(f,:)=UK_from_toW(f,:).*(ReductionSchool(2,f+1-f(1))'*ones(1,21));
new_UK_from_toS(f,:)=UK_from_toS(f,:).*(ReductionSchool(3,f+1-f(1))'*ones(1,21));
new_UK_from_toO(f,:)=UK_from_toO(f,:).*(ReductionSchool(4,f+1-f(1))'*OtherVector);

f=5:14; %adults
new_UK_from_toH(f,:)=UK_from_toH(f,:).*(ReductionAdult(1,f+1-f(1))'*ones(1,21));
new_UK_from_toW(f,:)=(1-workers_interacting_with_other)*UK_from_toW(f,:).*(ReductionAdult(2,f+1-f(1))'*ones(1,21)) ...
    + workers_interacting_with_other*UK_from_toW(f,:).*(ReductionAdult(2,f+1-f(1))'*OtherVector);
new_UK_from_toS(f,:)=UK_from_toS(f,:).*(ReductionAdult(3,f+1-f(1))'*ones(1,21));
new_UK_from_toO(f,:)=UK_from_toO(f,:).*(ReductionAdult(4,f+1-f(1))'*OtherVector);

f=15:21; %elderly
new_UK_from_toH(f,:)=UK_from_toH(f,:).*(ReductionElderly(1,f+1-f(1))'*ones(1,21));
new_UK_from_toW(f,:)=(1-workers_interacting_with_other)*UK_from_toW(f,:).*(ReductionElderly(2,f+1-f(1))'*ones(1,21)) ...
    + workers_interacting_with_other*UK_from_toW(f,:).*(ReductionElderly(2,f+1-f(1))'*OtherVector);
new_UK_from_toS(f,:)=UK_from_toS(f,:).*(ReductionElderly(3,f+1-f(1))'*ones(1,21));
new_UK_from_toO(f,:)=UK_from_toO(f,:).*(ReductionElderly(4,f+1-f(1))'*OtherVector);

end
