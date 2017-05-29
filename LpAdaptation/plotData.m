function plotData(cntVec, currentCount, currentSaveInd, muVec, rVec, VerboseModuloGen, N, eigVals, r, P_empVecAll, P_empVecWindow, P_empAll, P_empWindow)
persistent lastInd lastSave lastEigVals;

if any(mod(currentCount,VerboseModuloGen)==0)
    
    % Initialize data and figure (will be done only once)
    if isempty(lastInd) || currentCount == VerboseModuloGen
        lastInd = 1;
        lastSave = 1;
        lastEigVals = ones(N,1);
    end
    
    currSaves = lastSave:currentSaveInd-1;
    currIndices = cntVec(currSaves);
    currInd = currIndices(end);
    
    figure(42);
    hold on;
    
    subplot(2,2,1);
    plot(currIndices,muVec(currSaves,:));
    title('mean');
    grid on;
    hold on
    
    subplot(2,2,2);
    semilogy(currIndices,rVec(currSaves,:)');
    title(['Current step size: ',num2str(r)]);
    grid on;
    
    hold on
    
    subplot(2,2,3)
    
    plot(currIndices,P_empVecWindow(currSaves),'r');
    plot(currIndices,P_empVecAll(currSaves));
    title(['Current acceptance prob all: ',num2str(P_empAll),' window (red): ',num2str(P_empWindow)]);
    grid on;
    
    hold on
    
    subplot(2,2,4)
    eigValsD=diag(eigVals);
    if N==2
        semilogy([lastInd,currInd],[1./sort(lastEigVals),1./sort(eigValsD)]','-');
    else
        semilogy([lastInd,currInd],[1./sort(lastEigVals),1./sort(eigValsD)],'-');
    end
    title(['Condition of C: ',num2str(max(eigValsD)/min(eigValsD))]);
    grid on;
    
    lastEigVals = eigValsD;
    lastInd = currIndices(end);
    lastSave = currentSaveInd-1;
    
    drawnow
end
end
