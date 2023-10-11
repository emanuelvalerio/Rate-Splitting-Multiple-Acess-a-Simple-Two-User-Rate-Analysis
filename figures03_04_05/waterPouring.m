%% 
function [WF] = waterPouring(Ppk,g,Nt)
       initPowerAllocation = ((Ppk + sum(1./g))/Nt)   -  1./(g);
       while (length(find(initPowerAllocation < 0)) > 0)  % while there's negative power alocated ...
            indxNeg = find(initPowerAllocation < 0);       % store the indice of negative power allocated
            indxPos = find(initPowerAllocation >= 0);      % store the indice of positive or null power allocated
            nPositiveSubChannel = length(indxPos);         % stores the number of positive power allocated
            initPowerAllocation(indxNeg) = 0;              % change the value of negative power to zero.
            newSNR  = g(indxPos);                        % takes only the SNRs in which the positive power was allocated
            auxPowerAllocation = ((Ppk + sum(1./newSNR))/nPositiveSubChannel)- 1./(newSNR);
            initPowerAllocation(indxPos) = auxPowerAllocation; % adds to the vector of powers allocated from the current positive index
       end
       WF = initPowerAllocation;
      
end