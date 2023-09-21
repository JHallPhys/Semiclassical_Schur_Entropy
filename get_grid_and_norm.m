function [z_ps,Norm_hm]=get_grid_and_norm(N,t_final,q,p,F,omega,hbar_eff,gamma) 

z_ps=zeros(N,N,t_final);
Norm_hm=zeros(N,N,t_final);
z=zeros(3,1);
reverseStr = ''; % String for counter
display(' ') % Formating in terminal

for itt_q=1:N % Loop over initial conditions
  msg = sprintf('Constructing Normscape %d/%d', itt_q, N); % Tell you how many are done
  fprintf([reverseStr, msg]);
  reverseStr = repmat(sprintf('\b'), 1, length(msg));
   
    tic
    for itt_p = 1:N
      
      
        q_0=q(itt_q);
        p_0=p(itt_p);
        w_0=0;

        for itt_time=1:t_final

        
        t_i=2*pi*(t_final-1); 
        t_f=2*pi*t_final;
    
   
        z(1)=q_0;
        z(2)=p_0;
        z(3)=w_0;
       
        % Integrate
        
        [t,y] = ode89(@(t,y) dquartic_nh(t,y,F,omega,hbar_eff,gamma),[t_i t_f],z); % Integrate
       
        % Update

        q_0=y(end,1); 
        p_0=y(end,2);
        w_0=y(end,3);
   
        % Store

        z_ps(itt_p,itt_q,itt_time)=q_0+1i*p_0;
        Norm_hm(itt_p,itt_q,itt_time)=w_0;
       
        
        end
        
       
    
    end    
    
end
display(' ') % Formating in terminal


end