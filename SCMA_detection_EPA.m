function [symbol_indices_hat]=SCMA_detection_EPA(J,K,M,y,Nr,N0,H_matrix,F,max_iter,C)
pr=0.25*ones(M,J);
%q_i=zeros(M,J);%APP
%n=zeros(J);
%N_p=N0(e);
%initialization of FN-mean and var
m_f_ini=zeros(K*Nr,J);
v_f_ini=zeros(K*Nr,J);
% m_f_ini=complex(zeros(K*Nr,J));
% v_f_ini=complex(zeros(K*Nr,J));
for j=1:J    
     NE_j=find(F(:,j));%neigh REs
   for nr=1:Nr
       for k1=1:length(NE_j)
    v_f_ini((nr-1)*K+NE_j(k1),j)=1000; 
       end
   end
end
m_f=m_f_ini;% mean FN
v_f=v_f_ini;%mean VN
q=zeros(M,J);%N_APP
m_q=zeros(K,J);%mean APP
v_q=zeros(K,J);%Var APP
m_v=zeros(K*Nr,J);%mean VN
v_v=zeros(K*Nr,J);%Var VN
U=zeros(Nr*K,M*J);%FN update
%symbol_indices_hat=zeros(J,1);

% q=complex(zeros(M,J));%N_APP
% m_q=complex(zeros(K,J));%mean APP
% v_q=complex(zeros(K,J));%Var APP
% m_v=complex(zeros(K*Nr,J));%mean VN
% v_v=complex(zeros(K*Nr,J));%Var VN
% U=complex(zeros(K*Nr,M*J));%FN update
symbol_indices_hat=zeros(J,1);
for t=1:max_iter
   for nr=1:Nr
     m_f_nr=m_f((nr-1)*K+1:K*nr,:);
     v_f_nr=v_f((nr-1)*K+1:K*nr,:);
    for j=1:J
        pr_codebook=C((j-1)*K+1:j*K,:);
        del_j=find(F(:,j));%all REs connected to j-UE
        for k=1:length(del_j)
            pr_k=del_j(k);
            for m=1:M
                pr_codeword=pr_codebook(:,m);
                U((nr-1)*K+pr_k,(j-1)*M+m)=exp((-(abs(pr_codeword(pr_k)-m_f_nr(pr_k,j)))^2)/v_f_nr(pr_k,j));
            end
        end
    end
   end
    for j=1:J
        del_j=find(F(:,j));
        for m=1:M
            prod_term2=1;          
            for nr=1:Nr
                  prod_term1=1;                                 %Modified
                U_nr_j=U((nr-1)*K+1:nr*K,(j-1)*M+1:j*M);
                for i=1:length(del_j)
                     prod_term1=prod_term1* U_nr_j(del_j(i),m);
                end
              prod_term2=prod_term2*prod_term1;  
            end
            q(m,j)=pr(m,j)*prod_term2;
        end
    end
    % Normalization
    for j=1:J
        nf=0;
        for m=1:M
            nf=nf+q(m,j);
        end
        for m=1:M
            q(m,j)=q(m,j)/nf;
        end
    end
    %APP_Mean
    for j=1:J
        pr_codebook=C((j-1)*K+1:j*K,:);
        del_j=find(F(:,j));
        for k=1:length(del_j)
            pr_k=del_j(k);
            nf1=0;
            for m=1:M
                pr_codeword=pr_codebook(:,m);
                nf1=nf1+q(m,j)*pr_codeword(pr_k);
            end
            m_q(pr_k,j)=nf1;
        end
    end
    %APP_Variance
    for j=1:J
        pr_codebook=C((j-1)*K+1:j*K,:);

        del_j=find(F(:,j));
        for k=1:length(del_j)
            pr_k=del_j(k);
            nf1=0;
            for m=1:M
                pr_codeword=pr_codebook(:,m);
               % v_q(pr_k,j)=nf1+q(pr_k,j)*(abs(pr_codeword(pr_k)-m_q(pr_k,j)))^2;
               nf1=nf1+q(m,j)*(abs(pr_codeword(pr_k)-m_q(pr_k,j)))^2;
            end
            v_q(pr_k,j)=nf1;
        end
    end
   %VN mean and varianve update 
%    v_v_nr=zeros(K,J);
%    m_v_nr=zeros(K,J);
    for nr=1:Nr
        m_f_nr=m_f((nr-1)*K+1:K*nr,:);
        v_f_nr=v_f((nr-1)*K+1:K*nr,:);
        %v_v_nr=v_v((nr-1)*K+1:K*nr,:); 
            v_v_nr=zeros(K,J);
            m_v_nr=zeros(K,J);
    for j=1:J
        del_j=find(F(:,j));
        for k=1:length(del_j)
            pr_k=del_j(k);
            v_v_nr(pr_k,j)=((1/v_q(pr_k,j))-(1/v_f_nr(pr_k,j)))^-1;
            m_v_nr(pr_k,j)=((m_q(pr_k,j)/v_q(pr_k,j))-(m_f_nr(pr_k,j)/v_f_nr(pr_k,j)))*v_v_nr(pr_k,j);
        end
    end 
            v_v((nr-1)*K+1:nr*K,:)=v_v_nr;
            m_v((nr-1)*K+1:nr*K,:)=m_v_nr;
    end
   
       
    %FN mean and variance update
   for nr=1:Nr
          y_nr=y((nr-1)*K+1:nr*K);
          h_matrix_nr=H_matrix((nr-1)*K+1:K*nr,:); 
          m_v_nr=m_v((nr-1)*K+1:K*nr,:);
          v_v_nr=v_v((nr-1)*K+1:K*nr,:);
       for k=1:K
        del_k=find(F(k,:));%all UEs connect to k-RE      
        for j=1:length(del_k)
            pr_j=del_k(j);%present UE
            del_k_min_j=setdiff(del_k,pr_j);%extrinsic
            arg1=0;
            arg2=0;
            for i=1:length(del_k_min_j)
                del_kj= del_k_min_j(i);
                %for nr=1:Nr                    
                arg1=arg1+h_matrix_nr(k,del_kj)*m_v_nr(k,del_kj);
                arg2=arg2+((abs(h_matrix_nr(k,del_kj)))^2)*v_v_nr(k,del_kj);
                %end
            end
            m_f((nr-1)*K+k,pr_j)=(y_nr(k)-arg1)/h_matrix_nr(k,pr_j);
            v_f((nr-1)*K+k,pr_j)=(N0+arg2)/(abs(h_matrix_nr(k,pr_j)))^2;
        end
      end
    end
    
    
    
    
end
V_posterior=q;
max_values=max(V_posterior);
for j=1:J
    for m=1:M
        if V_posterior(m,j)==max_values(j)
            symbol_indices_hat(j)=m;
            break;
        end
    end
end    

       






            


