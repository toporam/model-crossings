function indxL1 = mkL1(Z4Rvect, Z4Lvect, nUnitsL1, density)

indxL1 = NaN(nUnitsL1*2,max(density));

for i = 1:nUnitsL1*2
    for d = 1:max(density)
        if i <= nUnitsL1 %if unit is in left hemisphere L1
            indx1 = mnrnd(1,Z4Rvect); %generate single location according to Z4Rvect probability distribution
        elseif i > nUnitsL1 %if unit is in right hemisphere L1
            indx1 = mnrnd(1,Z4Lvect);
        end
        indxL1(i,d) = find(indx1==1);
        clear indx1
    end 
end 