clear all;                                            
close all;
clc;
%���贫��10000������
N_symbol=10000;
%�Ŵ����Է��������� ����
h = 6.62606896*(1e-34);
f = 193e12;
B = 50e9;
NF = 4;

%���ŵ�λ����
Transmit_Power  = 1e-3;
Transmit_Amp = sqrt(Transmit_Power);
Span_num = (10:10:250);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80kmΪһ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
One_span_distance = 80;   %��λkm 
for span_i=1:length(Span_num)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%qpskԭʼ�ź�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %qpskԭʼ������λ��
    si_qpsk=(1/(2^(1/2)))*2*(round(rand(1,N_symbol))-0.5);                      
    sq_qpsk=(1/(2^(1/2)))*2*(round(rand(1,N_symbol))-0.5);                                    
    s_qpsk=si_qpsk + 1i*sq_qpsk;  
    s_qpsk = s_qpsk*Transmit_Amp;
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
    s_8qam = s_8qam*Transmit_Amp;
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
    s_16qam = s_16qam*Transmit_Amp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���紫��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r_qpsk = span(s_qpsk,One_span_distance,Span_num(span_i),h,f,B,NF)/Transmit_Amp; 
    r_8qam = span(s_8qam,One_span_distance,Span_num(span_i),h,f,B,NF)/Transmit_Amp;
    r_16qam = span(s_16qam,One_span_distance,Span_num(span_i),h,f,B,NF)/Transmit_Amp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %����δ���֮����źż���������
    
    %qpsk BER
    si__qpsk=sign(real(r_qpsk));                                
    sq__qpsk=sign(imag(r_qpsk));                               
    ber1_qpsk=(N_symbol-sum((si_qpsk/(1/(2^(1/2))))==si__qpsk))/N_symbol;                          
    ber2_qpsk=(N_symbol-sum((sq_qpsk/(1/(2^(1/2))))==sq__qpsk))/N_symbol;   
    ber_qpsk(span_i)=mean([ber1_qpsk ber2_qpsk]);        
    %8qam BER
    r_8qam_oct = qamdemod_8qam(r_8qam,8);             %��ʱ�����������8�����ź�
    r_8qam_bit_stream=de2bi(r_8qam_oct,3,'left-msb');              %ת��Ϊ��Ӧ�Ķ����Ʊ�����
    r_8qam_bit_stream=reshape(r_8qam_bit_stream.',numel(r_8qam_bit_stream),1');
    ber_8qam(span_i)=biterr(s_8qam_bit_stream,r_8qam_bit_stream)/(N_symbol*3);
    
    %16qam BER
    r_16qam_hex=demodulate(modem.qamdemod(16),r_16qam); %��ʱ�����������16�����ź�
    r_16qam_bit_stream=de2bi(r_16qam_hex,4,'left-msb');              %ת��Ϊ��Ӧ�Ķ����Ʊ�����
    r_16qam_bit_stream=reshape(r_16qam_bit_stream.',numel(r_16qam_bit_stream),1');
    ber_16qam(span_i)=biterr(s_16qam_bit_stream,r_16qam_bit_stream)/(N_symbol*4); 
end


figure(1);clf;hold on;
plot(Span_num,ber_qpsk,'r-*');
plot(Span_num,ber_8qam,'b-o');
plot(Span_num,ber_16qam,'k-x');
plot(Span_num,0.02*ones(size(Span_num)),'r--','linewidth',2);
legend('QPSK','8QAM','16QAM','Ber����ֵ');
xlabel('���˶���');ylabel('ber');
title('���˳��� 80Km');
set(gca,'yscale','log');
grid on; box on;

function [Out] = qamdemod_8qam(In,mod)
switch mod
    case 8
        ModSymmbol = [1+1j -1+1j -1-1j 1-1j 1i*3 -1i*(1+3^(1/2)) (1+3^(1/2)) -(1+3^(1/2))];
        for a=1:length(In)
            Dis = abs(ModSymmbol-In(a));
            [~,Ind] = min(Dis);
            Out(a) = Ind-1;
        end
end
end

function out = span(sgl,One_span,Span_num,h,f,B,NF)

sgl_tmp=sgl;
for a=1:Span_num
    sgl_tmp=span_opt(sgl_tmp,One_span,h,f,B,NF);
end
out = sgl_tmp;
end

function out = span_opt(sgl,One_span,h,f,B,NF)
%������˵�˥����Gainֵ
Loss_dB = 3*(One_span/15000);
Loss = 10^(Loss_dB/10);
Gain = 1/Loss;
%������������
Pn = 2*pi*h*f*B*(NF+1/Gain);

InPow = mean(sgl.^2)/1e3;
Pno = 2*pi*h*f*B*NF;
%��������������
Pnnl = (InPow/(1e-3))^2*Pno*2/3;

%��������ź�
NLNoise = (randn(size(sgl))+randn(size(sgl))*1i)/sqrt(2)*sqrt(Pnnl);
In2 = sgl + NLNoise;

%����˥���ͷŴ�
Out1 = In2/Loss*Gain;

%�Ŵ�������
LnaNoise = (randn(size(sgl))+randn(size(sgl))*1i)/sqrt(2)*sqrt(Pn);
out = Out1+LnaNoise;
end
















