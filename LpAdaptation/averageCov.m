function [out_averageCov,warningCell]= averageCov(out,numLast)

p_test = out.opts.pn;
dim = size(out.xAcc,2);

% save estimated volume from MC at each iteration (of the last numLast iter)
Vol_test_all = nan(numLast,1);
% save how close eigenvectors are to first set of eigenvectors
Vec_similar_first = nan(numLast-1,dim);
Vec_similar_previous = nan(numLast-1,dim);

% save average eigenvalues and eigenvectors
eigVal_m_previous = nan(dim,1);
eigVec_m_previous = nan(dim);

eigVal_m_first = nan(dim,1);
eigVec_m_first = nan(dim);

eigValAll_first = nan(numLast,dim);
eigValAll_previous = nan(numLast,dim);
warningCell=cell(numLast,1);
warningCnt = 1;
for k=1:numLast
    %numLast
    %size(out.QCell)
    
    %get MC Covariance in each step
    Q = out.QCell{end-numLast+k};
    r = out.rVec(end-numLast+k);
    C = r^2*(Q*Q');
    
    %estimated Volume from MC
    Vol_test_all(k) = Vol_lp(dim,r,p_test)*out.P_empVecAll(end-numLast+k);
    
    %eigen decomposition
    [eigVec,eigVals] = eig(C); 
    
    if k==1
        %to save average eigenvalues and eigenvectors
        eigVal_m_previous=diag(eigVals);
        eigVec_m_previous=eigVec;
        
        eigVal_m_first = diag(eigVals);
        eigVec_m_first = eigVec;
        
        %eigenvalues and eigenvectors from previous iteration
        eigVal_previous = diag(eigVals);
        eigVec_previous = eigVec;
        %eigenvalues and eigenvectors from first iteration
        eigVal_first = diag(eigVals);
        eigVec_first = eigVec;
        
        eigValAll_previous(k,:) = eigVal_previous;
        eigValAll_first(k,:) = eigVal_first;
    else
        check_allUsed_previous = zeros(dim,1);
        check_allUsed_first = zeros(dim,1);
        
        order_previous = zeros(dim,1);
        order_first = zeros(dim,1);
        
        %signs
        si = ones(dim,1);
        
        for v1=1:dim
            
            %for each eigenvector, find corresponding eigenvector from
            %previous iteration and from . 
            eigVec1 = eigVec(:,v1);
            
            %get vectors that are most alike from previous iteration
            [maxi_previous,idx_previous]=max(abs(eigVec1'*eigVec_previous)); 
            check_allUsed_previous(idx_previous) = check_allUsed_previous(idx_previous)+1;
            order_previous(v1)=idx_previous;
            
            %get vectors that are most alike from first iteration
            [maxi_first,idx_first] = max(abs(eigVec1'*eigVec_first));
            check_allUsed_first(idx_first) = check_allUsed_first(idx_first)+1;
            order_first(v1) = idx_first; 
            
            %add eigenvalues 
            diag_eigVals = diag(eigVals);
            eigVal_m_first(idx_first) = eigVal_m_first(idx_first)+diag_eigVals(v1);
            eigVal_m_previous(idx_previous) = eigVal_m_previous(idx_previous)+diag_eigVals(v1);
            
            
%             % check that vectors are pointing in same direction
%             % check angle between them
%             % multiply eigVec1 * -1 if angle too big
%             angle_v = acos(eigVec1'*eigVec_old(:,idx));
%             if angle_v >= pi/4
%                 eigVec1 = eigVec1*(-1);
%                 disp('here');
%             end
            % check that vectors are pointing in same direction
            projection_previous = eigVec1'*eigVec_previous(:,idx_previous); 
            if projection_previous < 0
                eigVec1_previous = eigVec1*(-1);
                si(v1)=-1;
            else
                eigVec1_previous = eigVec1;
            end
            
             projection_first = eigVec1'*eigVec_first(:,idx_first); 
            if projection_first < 0
                eigVec1_first = eigVec1*(-1);
            else
                eigVec1_first = eigVec1;
            end
            
            % add eigenvectors
            eigVec_m_previous(:,idx_previous)=eigVec_m_previous(:,idx_previous)+eigVec1_previous;
            eigVec_m_first(:,idx_first)=eigVec_m_first(:,idx_first)+eigVec1_first;
            
            %vec similar
            Vec_similar_previous(k-1,idx_previous) = maxi_previous;
            Vec_similar_first(k-1,idx_first) = maxi_first;
            
            
        end
        if ~isempty(find(check_allUsed_previous==0))% any(iszero(check_allUsed_previous))
            k
            check_allUsed_previous
            error('not all eigenvectors matched to');
        end
        if ~isempty(find(check_allUsed_first==0))%any(iszero(check_allUsed_first))
            k
            warn.k=k;
            warn.check_allUsed_first=check_allUsed_first;
            warn.order_previous = order_previous; 
            warn.order_first=order_first;
            warn.warning='when matching to the first iteration, at least one eigenvector is not matched to';
%             k
%             check_allUsed_first
%             order_previous
%             order_first
%             warning('when matching to the first iteration, at least one eigenvector is not matched to');
            warningCell{warningCnt} = warn;
            warningCnt = warningCnt+1; 
        end
        if ~all(order_previous == order_first)
            k
            warn.k=k;
            warn.order_previous = order_previous';
            warn.order_first = order_first';
            warn.warning = 'the matchings are not the same';
            
            warningCell{warningCnt} = warn;
            warningCnt = warningCnt+1; 
%             k
%             order_previous'
%             order_first'
%             warning('the matchings are not the same');
        end
        
        eigVal_previous = eigVals(order_previous);
        eigValAll_previous(k,:) = eigVal_previous;
        eigVec_previous = eigVec(:,order_previous); 
        eigValAll_first(k,:) = eigVal_first;
        
        for v1=1:dim
            if si(v1)==-1
                eigVec_previous(:,v1)= -1*eigVec_previous(:,v1);
            end
        end
%         order_previous'
%         order_first'
        
    end
    
       
end

%divide by numLast --> get mean
eigVec_previous = eigVec_m_previous./numLast;
eigVals_previous = eigVal_m_previous./numLast;
% eigVal_m_previous
% eigVal_m_first
eigVec_first = eigVec_m_first./numLast;
eigVals_first = eigVal_m_first./numLast;

%normalize eigvectors v/norm(v)
eigVec_previous = eigVec_previous./norm(eigVec_previous); % for k=1:dim eigVec2(:,k) = eigVec(:,k)./norm(eigVec(:,k));end
eigVec_first = eigVec_first./norm(eigVec_first);

%get covariance
% previous
C_previous = eigVec_previous *diag(eigVals_previous)* inv(eigVec_previous)
%symmetrize C
C_previous_sym = 0.5*(C_previous + C_previous');
% asymmetric part
C_previous_asym = 0.5*(C_previous - C_previous');

% first
C_first = eigVec_first * diag(eigVals_first) * inv(eigVec_first);
%symmetrize C
C_first_sym = 0.5*(C_first + C_first');
% asymmetric part
C_first_asym = 0.5*(C_first - C_first');

%
% r_norm_previous = nthroot(det(C_previous),dim*2);
% r_norm_first = nthroot(det(C_first),dim*2);
% 
% Vol_previous = Vol_lp(dim,r_norm_previous,p_test)*out.P_empAll;
% Vol_first = Vol_lp(dim,r_norm_first,p_test)*out.P_empAll;


out_averageCov.Vec_similar_previous = Vec_similar_previous;
out_averageCov.Vec_similar_first = Vec_similar_first; 
out_averageCov.Vol_test_all = Vol_test_all;
out_averageCov.C_previous = C_previous_sym; 
out_averageCov.C_first = C_first_sym; 
out_averageCov.C_previous_asym = C_previous_asym;
out_averageCov.C_first_asym = C_first_asym;
out_averageCov.p_test = p_test; 
out_averageCov.dim = dim; 
out_averageCov.eigValAll_first = eigValAll_first; 
out_averageCov.eigValAll_previous = eigValAll_previous; 

% if warningCnt > 1
%     out_warningCell = warningCell{1:(warningCnt-1),:};
% end
warningCell = warningCell(1:(warningCnt-1));
out_averageCov.warningCell = warningCell;