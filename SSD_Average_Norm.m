clear all 
close all

%==========================================================================
%   System Parameters
%==========================================================================

N=50;
t_final=10; 
gamma=0.001;
hbar_eff=0.12;
F =0.3; % The strength of the kicking term
omega=1;
q_R=1.5;
p_R=1.5; 
sigma_qq=hbar_eff;
sigma_pp=sigma_qq;
sigma_qp = [sigma_qq 0; 0 sigma_pp];
num_efn=30
%==========================================================================
%   Initialise grid and get norm landscape
%==========================================================================

[q,p,dq,dp,qmesh,pmesh]=init_classical_grid(q_R,p_R,N,hbar_eff); % Initalise rescaled grid
[zmesh,Norm_hm]=get_grid_and_norm(N,t_final,q,p,F,omega,hbar_eff,gamma); % Calculate dynamcis
[qmesh,pmesh]=meshgrid(q,p);
Norm_hm_av=zeros(N,N);
Norm_hm=exp(Norm_hm);
for j = 1:t_final
    Norm_hm_av=Norm_hm_av+Norm_hm(:,:,j);
end
Norm_hm_av=Norm_hm_av./t_final;

% figure
% imagesc(q,p,Norm_hm_av)
% colorbar
% title('classical Entropy')
% colormap(viridis)
% set(gca,'YDir','normal')
% xlabel('q')
% ylabel('p')
% c = colorbar('eastoutside');
% 
% return

%==========================================================================
%   Sort the norm landscape 
%==========================================================================
Normp=ones(N); % The full point set is the uniform distribution on the phase space
[i1,i2]=find(Normp==1); % The inidices of each point in the phase space
Norm_unsort=Norm_hm_av(:); % Take the unsorted normscape
Norm_hm_sort=sort(Norm_unsort,'descend'); % Sort the Normscape to get index set
% !!!!!!!!!

%The number of states need not be equal to N

tic
SS=zeros(N,N);
dz=dq*dp;
dcell=floor(2*pi/dz); % Maybe should be ceil
num_efn=floor(N^2/dcell)
for itt_state=1:num_efn
    
SS(:,:)=0;
SS=Norm_hm_av;
% Norm_hm_sort holds the partition parameters at each m*N element
SS(SS<=Norm_hm_sort(itt_state*dcell))=NaN; % contributions for all states < itt_state
if itt_state>1
    SS(SS>=Norm_hm_sort((itt_state-1)*dcell))=NaN; % Remove contributions from states not equal to itt_state
end
 SS(~isnan(SS))=1;
 SS(isnan(SS))=0;
 PS(:,:,itt_state)=SS;

end
toc


%==========================================================================
%   Now do the whole smoothing at once
%==========================================================================
tic
for ittq=1:length(q) % Integral dq
    ittq;
    for ittp=1:length(p)% Integrtal dp

%     phi2(1,:,:)=phin; % Store this way to 3D multiplication    
    mu = [q(ittq) p(ittp)]; % mu=[mux mup]
  
    G = mvnpdf([qmesh(:) pmesh(:)],mu,sigma_qp);
    G = reshape(G,length(q),length(p)); % What the fuck is this?

    CD(ittp,ittq,:)=sum(sum(PS.*G*dz));


    end
end
toc
clear PS SS
%==========================================================================
%   Now Calcualte entropy
%==========================================================================

Entropy=zeros(N,N);
CD_state=zeros(N,N);
for itt_state=1:num_efn
    CD_state=CD(:,:,itt_state);
    CD_state(CD_state==0)=1;
    Entropy=Entropy-CD_state.*log(CD_state);
end


figure
imagesc(q,p,Entropy)
colorbar
title('classical Entropy')
colormap(viridis)
set(gca,'YDir','normal')
xlabel('q')
ylabel('p')
c = colorbar('eastoutside');
