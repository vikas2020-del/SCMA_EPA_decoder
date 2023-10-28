 %clc;close all; 
load('ALC_DL_DO_150.mat','C')    % Load the codebook

M=size(C,2);  % size of the codebook for one user
m=log2(M); % no of bits in a symbol/block
K=4;   % No of orthogonal PREs/ OFDM subcarriers
J=6;   % No of users/layers
Nr=16;
F=get_indicator_matrix_factor_graph(C, J, K);
F1=zeros(Nr*K,J);
for nr=1:Nr
    
F1((nr-1)*K+1:nr*K,:)=F;
    
end
%% %%%%%%% Power Allocation %%%%%%%%%
power_users=ones(1,J);  
sqrt_power_matrix=ones(K,J);
%% Monte Carlo %%%%%%%%%%%%%%%%%%%%%%%%
%Eb_N0_dB=5:3:30;
%Eb_N0_dB=-5:2:5;
Eb_N0_dB=-9:2:1;
Eb_N0=10.^(Eb_N0_dB/10);
Es=sum(sum((abs(C)).^2))/length(C);
Eb=Es/m;   
N0=Eb./Eb_N0;
sigma=sqrt(N0/2);
max_block_errors=[  100 70 50 20 10 5];
%max_block_errors= 1  ;
%max_block_errors=100;
max_iter=10;  % maximum number of iterations for SCMA detection
SER=zeros(1,length(Eb_N0_dB));
%cf=500; % LLR clipping factor (used only in the LLR-domain MPA)
%% 
 for e=1:length(Eb_N0_dB)

    block_errors=0;
    symbol_errors=0;
    block=0;
    while block_errors<max_block_errors(e)
        %%   SCMA Encoding %%%%%%
        bits=randi([0 1],J,m);% blocks of bits for all users
        %mapping
        symbol_indices=bi2de(bits,'left-msb')+1; % symbols for all users
        SCMA_codewords=zeros(K,J);  % collection of the codewords for all users
        for j=1:J         % for each user
            present_codebook=C((j-1)*K+1:j*K,:);   % codebook for the jth user
            SCMA_codewords(:,j)=present_codebook(:,symbol_indices(j));
        end
        SCMA_codewords_MA=zeros(Nr*K,J);
        H_matrix=zeros(Nr*K,J);
        for nr=1:Nr
        SCMA_codewords_MA((nr-1)*K+1:nr*K,:)=SCMA_codewords;
        H_matrix((nr-1)*K+1:nr*K,:)=1/sqrt(2)*(randn(K,J)+1j*randn(K,J));
        end
        %% Transmission through Rayleigh fading channel %%
        AWGN=sigma(e)*(randn(K*Nr,1)+1j*randn(K*Nr,1));    % complex Gaussian noise
        %AWGN=0;
        
       % h_matrix=1/sqrt(2)*(randn(K,J)+1j*randn(K,J));  % complex Rayleigh fading coefficient vector for all users
        %h_matrix=ones(K,J);   % Simple AWGN channel
        %h=1/sqrt(2)*(randn(K,1)+1j*randn(K,1));
       %h_matrix=repmat(h,1,J);% complex Rayleigh fading coefficient vector for Down-link        
        y=sum(H_matrix.*SCMA_codewords_MA,2)+AWGN;
        %y=reshape(y,[K,Nr]);
        %% received SCMA codeword UP_link         
       
        %[Rx_ant_prod,T1,T2,V,U,V_posterior,symbol_indices_hat]=SCMA_detection_MPA_prob_domain_MU_MIMO(power_users,y,N0(e),H_matrix,F,max_iter,M,C,Nr);
        %[U,V,arg1,T1,T2,V_posterior,symbol_indices_hat]=SCMA_detection_MPA_LLR_domain(power_users,y,N0(e),H_matrix,F1,max_iter,M,C,cf,Nr);
        %[ U,q,m_q,v_q,m_v,v_v,m_f,v_f,V_posterior,symbol_indices_hat]=SCMA_detection_EPA2(J,K,M,y,Nr,N0(e),H_matrix,F,max_iter,C);
        [symbol_indices_hat]=SCMA_detection_EPA(J,K,M,y,Nr,N0(e),H_matrix,F,max_iter,C);
        error_locations=find(symbol_indices~=symbol_indices_hat);
        %demapping
        %symbol_indices_hat=symbol_indices_hat-1;
%         bits_hat=de2bi(symbol_indices_hat,'left-msb');
%         bits_error=sum(sum(bits~=bits_hat));
        if ~isempty(error_locations)
            block_errors=block_errors+1;
            symbol_errors=symbol_errors+length(error_locations);
            %fprintf('\n  %d error collected',block_errors); 
        end     
        block=block+1;
    end
    SER(e)=symbol_errors/block/J;% Eack block contains J-symbols
    %BER(e)=SER(e)/m;
    fprintf('\n Simulation done for %d dB',Eb_N0_dB(e)); 
end
semilogy(Eb_N0_dB,SER,'b-*','LineWidth',2) ;
xlabel('Eb/No');
ylabel('SER');
legend('Nr=4,M=4,200 ');





