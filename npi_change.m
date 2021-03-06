function [aC,new_UK_from_toH, new_UK_from_toW, new_UK_from_toS, new_UK_from_toO] = npi_change(S_New_Var,Impact_of_New_Var,W,T,Region,ComplianceT,Region_PP,PD_Lockdown,Run_stop, UK_from_toH, UK_from_toW, UK_from_toS, UK_from_toO)
%% Function to update mixing dependent on current compliance and effect of current level of evolved more contagious variants

% Inputs:
% S_New_Var - Fraction of virus in circulation that is more transmissible variant
% Impact_of_New_Var - Increase in transmissability due to new varaint
% W - Current index for accessing adherence array
% T - Time
% Region - ID for current region being analysed
% ComplianceT - Adherence to NPIs
% Region_PP - Population per five-year age group in the current region
% Run_stop - Time at which this phase of NPI strength ends.
% UK_from_toH, UK_from_toW, UK_from_toS, UK_from_toO -
%   The original mixing matrices/contact arrays.
%   Order is Household, Work, School, Other

% Outputs:
% aC - Amended adherence values
% new_UK_from_toH, new_UK_from_toW, new_UK_from_toS, new_UK_from_toO -
%   The amended mixing matrices/contact arrays.
%   Order is Household, Work, School, Other


% Initialise amended adherence vector (aC)
mmC=max(max(abs(ComplianceT)));
if mmC>0
    if length(ComplianceT(:,W))==1
        aC=ones(21,1)*abs(ComplianceT(:,W));
    end
    if length(ComplianceT(:,W))==4
        aC=abs(ComplianceT([1 2 2 2 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4],W));
    end
    if length(ComplianceT(:,W))==21
        aC=reshape(abs(ComplianceT(:,W)),21,1);
    end

end

% Based on time T, specify if schools are open by binary variable InSchoolFlag
InSchoolFlag=0;
if Region==10
    if T(end)>=230  % ie Mid-August put schools back in Scotland
        InSchoolFlag=1;
    end
    if T(end)>=292 && T(end)<=294
        InSchoolFlag=0;  % Unless its Half Term
    end
    if (T(end)>=355 && T(end)<=433)
        InSchoolFlag=0;   % ... EXTENDED CHRISTMAS BREAK
    end
else
    if T(end)>=244 % ie 1st Sept put schools back
        InSchoolFlag=1;
    end
    if T(end)>=299 && T(end)<=301
        InSchoolFlag=0;   % Unless its Half Term
    end

    if (T(end)>=355 && T(end)<=433)
        InSchoolFlag=0;   % ... EXTENDED CHRISTMAS BREAK
    end
end

if T(end)>=458 && T(end)<472
    InSchoolFlag=0; % Easter
end

% If schools are open, amend School age groups aC (age groups 2-4, so 5-19yrs)
if InSchoolFlag & ComplianceT(1,W)>=0
    aC(2:4)=0.2*aC(2:4);
end


PD_Lockdown=PD_Lockdown+(Run_stop(W)-T(end)).*(Region_PP(Region,:)*aC)./mmC;

% Get revised mixing matrices
[new_UK_from_toH, new_UK_from_toW, new_UK_from_toS, new_UK_from_toO] = Return_New_Matrices(0.95, 0.8, 0.95, aC, UK_from_toH, UK_from_toW, UK_from_toS, UK_from_toO);

% Compute scaling factor due to new variant
Multi_Factor = 1 + S_New_Var(T(end)) * Impact_of_New_Var;  % where T(end) should be the current time.

% Scale contacts according to transmissability scaling factor of variant
new_UK_from_toH = new_UK_from_toH * Multi_Factor;
new_UK_from_toW = new_UK_from_toW * Multi_Factor;
new_UK_from_toS = new_UK_from_toS * Multi_Factor;
new_UK_from_toO = new_UK_from_toO * Multi_Factor;

end
