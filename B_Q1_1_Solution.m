clear all;                                            
close all;
%���贫��10000������
N_symbol=10000;

for snrdb=1:1:25
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%qpskԭʼ�ź�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %qpskԭʼ������λ��
    si_qpsk=(1/(2^(1/2)))*2*(round(rand(1,N_symbol))-0.5);                      
    sq_qpsk=(1/(2^(1/2)))*2*(round(rand(1,N_symbol))-0.5);                                    
    s_qpsk=si_qpsk + 1i*sq_qpsk;  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%8qamԭʼ�ź�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dot_8qam = [1+1j -1+1j -1-1j 1-1j 1i*3 -1i*(1+3^(1/2)) (1+3^(1/2)) -(1+3^(1/2))];  %8qam����ӳ���
    %������������% ӳ��16qam
    rand0_7 = randi(8,1,N_symbol);
    s_8qam = zeros(1,N_symbol);   
    for a=1:N_symbol
        s_8qam(1,a)=dot_8qam(1,rand0_7(1,a));
    end
    %8qamԭʼ������λ��
    s_8qam_oct=qamdemod_8qam(s_8qam,8); %��ʱ�����������8�����ź�
    s_8qam_bit_stream=de2bi(s_8qam_oct,3,'left-msb');              %ת��Ϊ��Ӧ�Ķ����Ʊ�����
    s_8qam_bit_stream=reshape(s_8qam_bit_stream.',numel(s_8qam_bit_stream),1');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%16qamԭʼ�ź�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dot_16qam = [1+1j -1+1j -1-1j 1-1j 3+1j 3+3j 1+3j -1+3j -3+3j -3+1j -3-1j -3-3j -1-3j 1-3j 3-3j 3-1j];  %16qam����ӳ���
    %������������% ӳ��16qam
    rand0_15 = randi(16,1,N_symbol);
    s_16qam = zeros(1,N_symbol);   
    for a=1:N_symbol
        s_16qam(1,a)=dot_16qam(1,rand0_15(1,a));
    end
    %16qamԭʼ������λ��
    s_16qam_hex=demodulate(modem.qamdemod(16),s_16qam); %��ʱ�����������16�����ź�
    s_16qam_bit_stream=de2bi(s_16qam_hex,4,'left-msb');              %ת��Ϊ��Ӧ�Ķ����Ʊ�����
    s_16qam_bit_stream=reshape(s_16qam_bit_stream.',numel(s_16qam_bit_stream),1');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%�Ӹ�˹������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %qpsk�Ӹ�˹������
    w_qpsk = awgn(s_qpsk,snrdb,'measured');
    r_qpsk = w_qpsk;       
    %8qam�Ӹ�˹������
    w_8qam = awgn(s_8qam,snrdb,'measured');
    r_8qam = w_8qam;
    %16qam�Ӹ�˹������
    w_16qam = awgn(s_16qam,snrdb,'measured');
    r_16qam = w_16qam;
    %�Լ�����֮����źż���������
    %qpsk BER
    si__qpsk=sign(real(r_qpsk));                                
    sq__qpsk=sign(imag(r_qpsk));                               
    ber1_qpsk=(N_symbol-sum((si_qpsk/(1/(2^(1/2))))==si__qpsk))/N_symbol;                          
    ber2_qpsk=(N_symbol-sum((sq_qpsk/(1/(2^(1/2))))==sq__qpsk))/N_symbol;   
    ber_qpsk(snrdb)=mean([ber1_qpsk ber2_qpsk]);            
    %8qam BER
    r_8qam_oct=qamdemod_8qam(r_8qam,8);                          %��ʱ�����������8�����ź�
    r_8qam_bit_stream=de2bi(r_8qam_oct,3,'left-msb');              %ת��Ϊ��Ӧ�Ķ����Ʊ�����
    r_8qam_bit_stream=reshape(r_8qam_bit_stream.',numel(r_8qam_bit_stream),1');
    ber_8qam(snrdb)=biterr(s_8qam_bit_stream,r_8qam_bit_stream)/(N_symbol*3);
    
    %16qam BER
    r_16qam_hex=demodulate(modem.qamdemod(16),r_16qam); %��ʱ�����������16�����ź�
    r_16qam_bit_stream=de2bi(r_16qam_hex,4,'left-msb');              %ת��Ϊ��Ӧ�Ķ����Ʊ�����
    r_16qam_bit_stream=reshape(r_16qam_bit_stream.',numel(r_16qam_bit_stream),1');
    ber_16qam(snrdb)=biterr(s_16qam_bit_stream,r_16qam_bit_stream)/(N_symbol*4); 
end
%QPSK����ͼ
figure(1);clf;hold on;
subplot(1,2,1);plot(s_qpsk,'b.');title('QPSK��������ͼ');axis equal;axis([-2 2 -2 2])
subplot(1,2,2);plot(r_qpsk,'b.');title('QPSK��������ͼ(25db)');axis equal;axis([-2 2 -2 2])

%8qam����ͼ
hold on;
figure(2);clf;hold on;
subplot(1,2,1);plot(s_8qam,'b.');title('8QAM��������ͼ');axis equal;axis([-5 5 -5 5])
subplot(1,2,2);plot(r_8qam,'b.');title('8QAM��������ͼ(25db)');axis equal;axis([-5 5 -5 5])
%16qam����ͼ

figure(3);clf;hold on;
subplot(1,2,1);plot(s_16qam,'b.');title('16QAM��������ͼ');axis equal;axis([-5 5 -5 5])
subplot(1,2,2);plot(r_16qam,'b.');title('16QAM��������ͼ(25db)');axis equal;axis([-5 5 -5 5])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%qpsk 8qam 16qam BERvsSNR%%%%%%%%%%%%%%%%%%%%%%%
figure(4);
snrdb=1:1:25;
semilogy(snrdb,ber_qpsk,'-bo',snrdb,ber_8qam,'-mh',snrdb,ber_16qam,'-gd',snrdb,0.02*ones(size(snrdb)),'r--','linewidth',2)
title('QPSK 8QAM 16QAM  with awgn');
xlabel('Signal to noise ratio');
ylabel('Bit error rate'); 
legend('QPSK','8QAM','16QAM');
axis([-2 25 -5 2]);
grid on;

function [out] = qamdemod_8qam(In,mod)
switch mod
    case 8
        Symmbol_8qam = [1+1j -1+1j -1-1j 1-1j 1i*3 -1i*(1+3^(1/2)) (1+3^(1/2)) -(1+3^(1/2))];
        for a=1:length(In)
            distance = abs(Symmbol_8qam-In(a));
            [~,Ind] = min(distance);
            out(a) = Ind-1;
        end
end
end



