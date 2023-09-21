function CD=get_cdensity(qmesh,pmesh,q,p,dq,dp,sigma_qp,Norm_hm_0,nfq,delta_upper_old,delta_lower_old,itt_max,eps_stop,efn_counter,efn_max)
dz=dq*dp;

display(' ') % Formating in terminal
msg=sprintf('State %d/%d', efn_counter, efn_max); % Tell you how many are done
reverseStr = ''; % String for counter
fprintf([reverseStr,msg]);
display(' ') % Formating in terminal
for itt=1:itt_max

  msg = sprintf('Constructing Density Itteration %d/%d', itt, itt_max); % Tell you how many are done
  fprintf([reverseStr, msg]);
  reverseStr = repmat(sprintf('\b'), 1, length(msg));
   
    % Reset arrays
    dnew=abs(delta_upper_old-delta_lower_old)/2; 
    delta_upper_new=delta_upper_old-dnew; % Move the right hand side 
    % Partition norm landscape
    eps_fwd=delta_upper_new;
    eps_bwd=eps_fwd;
    [Normp] = partition_loss(Norm_hm_0,eps_fwd,eps_bwd,'G');
    for ittq=1:length(q) % Integral dq
        ittq;
        for ittp=1:length(p)% Integrtal dp
%             ittp=26
    
            mu = [q(ittq) p(ittp)]; % mu=[mux mup]
          
            G = mvnpdf([qmesh(:) pmesh(:)],mu,sigma_qp);
            G = reshape(G,length(q),length(p)); % What the fuck is this?
%               G=G*sqrt(det(sigma_qp))*2*pi;
            CD(ittp,ittq)=sum(sum(Normp.*G*dz));
    

        end
    end

    nfc=sum(sum(CD*dq*dp));
  

        if abs(nfq-nfc)<eps_stop % Passes
                display('FINISHED ')
                break
        else % Fails and update interval
            if nfc>nfq
                delta_lower_old=delta_upper_new;
            elseif nfc<nfq
                delta_upper_old=delta_upper_new;
            end
        end

end


end