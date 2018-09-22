clear all;                                            
close all;
clc;
%假设传输10000个符号
N_symbol=10000;
%放大器自发辐射噪声 常数
h = 6.62606896*(1e-34);
f = 193e12;
B = 50e9;
NF = 4;

%符号单位功率
Transmit_Power  = 1e-3;
Transmit_Amp = sqrt(Transmit_Power);
Span_num = (10:10:250);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%80km为一跨%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
One_span_distance = 80;   %单位km 
for span_i=1:length(Span_num)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%qpsk原始信号%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %qpsk原始二进制位流
    si_qpsk=(1/(2^(1/2)))*2*(round(rand(1,N_symbol))-0.5);                      
    sq_qpsk=(1/(2^(1/2)))*2*(round(rand(1,N_symbol))-0.5);                                    
    s_qpsk=si_qpsk + 1i*sq_qpsk;  
    s_qpsk = s_qpsk*Transmit_Amp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%8qam原始信号%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dot_8qam = [1+1j -1+1j -1-1j 1-1j 1i*3 -1i*(1+3^(1/2)) (1+3^(1/2)) -(1+3^(1/2))];  %8qam星座映射点
    %产生符号序列% 映射16qam
    rand0_7 = randi(8,1,N_symbol);
    s_8qam = zeros(1,N_symbol);   
    for a=1:N_symbol
        s_8qam(1,a)=dot_8qam(1,rand0_7(1,a));
    end
    %8qam原始二进制位流
    s_8qam_oct=qamdemod_8qam(s_8qam,8); %此时解调出来的是8进制信号
    s_8qam_bit_stream=de2bi(s_8qam_oct,3,'left-msb');              %转化为对应的二进制比特流
    s_8qam_bit_stream=reshape(s_8qam_bit_stream.',numel(s_8qam_bit_stream),1');
    s_8qam = s_8qam*Transmit_Amp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%16qam原始信号%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dot_16qam = [1+1j -1+1j -1-1j 1-1j 3+1j 3+3j 1+3j -1+3j -3+3j -3+1j -3-1j -3-3j -1-3j 1-3j 3-3j 3-1j];  %16qam星座映射点
    %产生符号序列% 映射16qam
    rand0_15 = randi(16,1,N_symbol);
    s_16qam = zeros(1,N_symbol);   
    for a=1:N_symbol
        s_16qam(1,a)=dot_16qam(1,rand0_15(1,a));
    end
    %16qam原始二进制位流
    s_16qam_hex=demodulate(modem.qamdemod(16),s_16qam); %此时解调出来的是16进制信号
    s_16qam_bit_stream=de2bi(s_16qam_hex,4,'left-msb');              %转化为对应的二进制比特流
    s_16qam_bit_stream=reshape(s_16qam_bit_stream.',numel(s_16qam_bit_stream),1');
    s_16qam = s_16qam*Transmit_Amp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%经跨传输%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r_qpsk = span(s_qpsk,One_span_distance,Span_num(span_i),h,f,B,NF)/Transmit_Amp; 
    r_8qam = span(s_8qam,One_span_distance,Span_num(span_i),h,f,B,NF)/Transmit_Amp;
    r_16qam = span(s_16qam,One_span_distance,Span_num(span_i),h,f,B,NF)/Transmit_Amp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %经跨段传输之后的信号计算误码率
    
    %qpsk BER
    si__qpsk=sign(real(r_qpsk));                                
    sq__qpsk=sign(imag(r_qpsk));                               
    ber1_qpsk=(N_symbol-sum((si_qpsk/(1/(2^(1/2))))==si__qpsk))/N_symbol;                          
    ber2_qpsk=(N_symbol-sum((sq_qpsk/(1/(2^(1/2))))==sq__qpsk))/N_symbol;   
    ber_qpsk(span_i)=mean([ber1_qpsk ber2_qpsk]);        
    %8qam BER
    r_8qam_oct = qamdemod_8qam(r_8qam,8);             %此时解调出来的是8进制信号
    r_8qam_bit_stream=de2bi(r_8qam_oct,3,'left-msb');              %转化为对应的二进制比特流
    r_8qam_bit_stream=reshape(r_8qam_bit_stream.',numel(r_8qam_bit_stream),1');
    ber_8qam(span_i)=biterr(s_8qam_bit_stream,r_8qam_bit_stream)/(N_symbol*3);
    
    %16qam BER
    r_16qam_hex=demodulate(modem.qamdemod(16),r_16qam); %此时解调出来的是16进制信号
    r_16qam_bit_stream=de2bi(r_16qam_hex,4,'left-msb');              %转化为对应的二进制比特流
    r_16qam_bit_stream=reshape(r_16qam_bit_stream.',numel(r_16qam_bit_stream),1');
    ber_16qam(span_i)=biterr(s_16qam_bit_stream,r_16qam_bit_stream)/(N_symbol*4); 
end


figure(1);clf;hold on;
plot(Span_num,ber_qpsk,'r-*');
plot(Span_num,ber_8qam,'b-o');
plot(Span_num,ber_16qam,'k-x');
plot(Span_num,0.02*ones(size(Span_num)),'r--','linewidth',2);
legend('QPSK','8QAM','16QAM','Ber门限值');
xlabel('光纤段数');ylabel('ber');
title('光纤长度 80Km');
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
%计算光纤的衰减和Gain值
Loss_dB = 3*(One_span/15000);
Loss = 10^(Loss_dB/10);
Gain = 1/Loss;
%线性噪声功率
Pn = 2*pi*h*f*B*(NF+1/Gain);

InPow = mean(sgl.^2)/1e3;
Pno = 2*pi*h*f*B*NF;
%非线性噪声功率
Pnnl = (InPow/(1e-3))^2*Pno*2/3;

%生成入跨信号
NLNoise = (randn(size(sgl))+randn(size(sgl))*1i)/sqrt(2)*sqrt(Pnnl);
In2 = sgl + NLNoise;

%进行衰减和放大
Out1 = In2/Loss*Gain;

%放大器噪声
LnaNoise = (randn(size(sgl))+randn(size(sgl))*1i)/sqrt(2)*sqrt(Pn);
out = Out1+LnaNoise;
end
















