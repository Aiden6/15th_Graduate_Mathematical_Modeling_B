clear all;clc;
close all;

New_constellation=Genetic_algorithm();

%���贫��10000������
N_symbol=10000;

for snrdb=1:1:25
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%8qamԭʼ�ź�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dot_8qam = [1+1j -1+1j -1-1j 1-1j 1i*3 -1i*(1+3^(1/2)) (1+3^(1/2)) -(1+3^(1/2))];  %8qam����ӳ���
    %������������% ӳ��16qam
    rand0_7 = randi(8,1,N_symbol);
    s_8qam_old = zeros(1,N_symbol);   
    for a=1:N_symbol
        s_8qam_old(1,a)=dot_8qam(1,rand0_7(1,a));
    end
    %8qamԭʼ������λ��
    s_8qam_oct_old=qamdemod_8qam(s_8qam_old,8); %��ʱ�����������8�����ź�
    s_8qam_bit_stream_old=de2bi(s_8qam_oct_old,3,'left-msb');              %ת��Ϊ��Ӧ�Ķ����Ʊ�����
    s_8qam_bit_stream_old=reshape(s_8qam_bit_stream_old.',numel(s_8qam_bit_stream_old),1');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%8qamԭʼ�ź�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s_8qam_oct_new = randi(8,1,N_symbol)-1;
    %8qam_newԭʼ������λ�� 
    s_8qam_new = qammod_new_constellation(s_8qam_oct_new,New_constellation);%%%%%%%%%%%���µ�����ͼ�����ͷ���
    s_8qam_bit_stream_new=de2bi(s_8qam_oct_new,3,'left-msb');              %ת��Ϊ��Ӧ�Ķ����Ʊ�����
    s_8qam_bit_stream_new=reshape(s_8qam_bit_stream_new.',numel(s_8qam_bit_stream_new),1');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %8qam�Ӹ�˹������
    w_8qam_old = awgn(s_8qam_old,snrdb,'measured');
    r_8qam_old = w_8qam_old;
    
    %8qam�Ӹ�˹������
    w_8qam_new = awgn(s_8qam_new,snrdb,'measured');
    r_8qam_new = w_8qam_new;
    
     %8qam BER
    r_8qam_oct_old=qamdemod_8qam(r_8qam_old,8);                          %��ʱ�����������8�����ź�
    r_8qam_bit_stream_old=de2bi(r_8qam_oct_old,3,'left-msb');              %ת��Ϊ��Ӧ�Ķ����Ʊ�����
    r_8qam_bit_stream_old=reshape(r_8qam_bit_stream_old.',numel(r_8qam_bit_stream_old),1');
    ber_8qam_old(snrdb)=biterr(s_8qam_bit_stream_old,r_8qam_bit_stream_old)/(N_symbol*3);
    
    %8qam BER
    r_8qam_oct_new=qamdemod_new_constellation(r_8qam_new,New_constellation); %��ʱ�����������8�����ź�
    r_8qam_bit_stream_new=de2bi(r_8qam_oct_new,3,'left-msb');              %ת��Ϊ��Ӧ�Ķ����Ʊ�����
    r_8qam_bit_stream_new=reshape(r_8qam_bit_stream_new.',numel(r_8qam_bit_stream_new),1');
    ber_8qam_new(snrdb)=biterr(s_8qam_bit_stream_new,r_8qam_bit_stream_new)/(N_symbol*3);

end
%8qam_old����ͼ
hold on;
figure(4);clf;hold on;
subplot(1,2,1);plot(s_8qam_old,'b.');title('8QAM��������ͼ');axis equal;axis([-5 5 -5 5])
subplot(1,2,2);plot(r_8qam_old,'b.');title('8QAM��������ͼ(25db)');axis equal;axis([-5 5 -5 5])
%8qam_new����ͼ
hold on;
figure(2);clf;hold on;
subplot(1,2,1);plot(s_8qam_new,'b.');title('8QAM��������ͼ');axis equal;axis([-5 5 -5 5])
subplot(1,2,2);plot(r_8qam_new,'b.');title('8QAM��������ͼ(25db)');axis equal;axis([-5 5 -5 5])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%8qam BERvsSNR%%%%%%%%%%%%%%%%%%%%%%%
figure(5);
snrdb=1:1:25;
semilogy(snrdb,ber_8qam_old,'-bo',snrdb,ber_8qam_new,'-mh')
title('8QAM OLD 16QAM NEW  with awgn');
xlabel('Signal to noise ratio');
ylabel('Bit error rate'); 
legend('8QAM OLD','8QAM NEW');
axis([-2 25 -5 2]);
grid on;

function New_constellation = Genetic_algorithm()
%��Ⱥ��С
population_individual_num = 50;
%�����Ʊ��볤��
code_length=8;
%�������
p_crossover = 0.1;
%�������
p_mutae = 0.02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ʼ��Ⱥ
Initial_population = Init_population(population_individual_num,code_length);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�洢���������
save_fit_value =0;
New_constellation = [];

Evolution_popualtion = Initial_population;
for generation = 1:200000
    %������Ӧ��ֵ������ֵ��
    evaluate_value = cal_evaluate_value(Evolution_popualtion);
    fit_value = evaluate_value;
    %ѡ��
    new_population = select(Evolution_popualtion,fit_value);
    %����
    new_population = crossover(new_population,p_crossover);
    %����
    new_population = mutate(new_population,p_mutae);
    %����
    Evolution_popualtion = update_popualtion(new_population);
    
    %Ѱ�����Ž�
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %���������Ӧ��
    [bestindividual,best_fit_value] = Find_Best_Individual(Evolution_popualtion,fit_value);
    
    if (save_fit_value<best_fit_value)
       save_fit_value = best_fit_value;
       New_constellation = bestindividual.*3;
    end
    
    Fitness(generation) = (best_fit_value);
    if mod(generation,100) == 0
        figure(1);clf
        plot(New_constellation,'x');
        axis([-5 5 -5 5]);
        title(['��������Ϊn=' num2str(generation)]);
        pause(0.01);
    end
end

figure(2);clf
plot(New_constellation,'x');
figure(3);
plot(Fitness)
end
function Initial_population=Init_population(popsize,chromlength)
Initial_population = (rand(popsize,chromlength)*2-1)+(rand(popsize,chromlength)*2-1)*1i;
end

function [population] = update_popualtion(population)
population_real=real(population);
population_imag=imag(population);
population_real = max(min(population_real,1),-1);
population_imag = max(min(population_imag,1),-1);
population = population_real+population_imag*1i;
end


function [evaluate_value] = cal_evaluate_value(population)
[individual,len] = size(population);
for kk=1:individual
    Dis = inf(len,len);
    for ii=1:len
        for jj=ii+1:len
            Dis(ii,jj) = abs(population(kk,ii)-population(kk,jj));
        end
    end
    minDis = min(Dis(:));
    if (minDis~=0)
        population(kk,:) = population(kk,:)/minDis;
        evaluate_value(kk) = 1/(mean(abs(population(kk,:)).^2));
    else
        evaluate_value (kk)= 0;
    end
end

end
function [new_population] = select(population,individual_fitness)
[px,~] = size(population);
total_fitness = sum(individual_fitness);
p_fitness = individual_fitness/total_fitness;
p_fitness = cumsum(p_fitness);%�����������
ms = sort(rand(px,1));%��С��������
fitin = 1;
newin = 1;
while newin<=px
    if(ms(newin))<p_fitness(fitin)
        new_population(newin,:)=population(fitin,:);
        newin = newin+1;
    else
        fitin=fitin+1;
    end
end
end

function [new_population] = crossover(population,pc)
[population_x,population_y] = size(population);
new_population = ones(size(population));
for i = 1:2:population_x-1
    if(rand<pc)
        cross_num = round(rand*population_y);
        new_population(i,:) = [population(i,1:cross_num),population(i+1,cross_num+1:population_y)];
        new_population(i+1,:) = [population(i+1,1:cross_num),population(i,cross_num+1:population_y)];
    else
        new_population(i,:) = population(i,:);
        new_population(i+1,:) = population(i+1,:);
    end
end
end


function [new_populatin] = mutate(population,p_mutate)
[population_x,population_y] = size(population);
new_populatin = ones(size(population));
for i = 1:population_x
    mutate_num = ceil(rand*population_y);
    if(rand<p_mutate && ~isempty(mutate_num))        
        new_populatin(i,:) = population(i,:);
        new_populatin(i,mutate_num) =  new_populatin(i,mutate_num)+...
            (randn(1,length(mutate_num))+randn(1,length(mutate_num))*1i)*0.05;
    else
        new_populatin(i,:) = population(i,:);
    end
end
end

function [best_individual, best_fit_value] = Find_Best_Individual(popualtion,individual_fitness_value)
[population_x,~] = size(popualtion);
best_individual = popualtion(1,:);
best_fit_value = individual_fitness_value(1);
for i = 2:population_x
    if individual_fitness_value(i)>best_fit_value
        best_individual = popualtion(i,:);
        best_fit_value = individual_fitness_value(i);
    end
end
end


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

function Out = qammod_new_constellation(xIn,ModSymmbol)
Out = ModSymmbol(xIn+1);
end

function Out = qamdemod_new_constellation(In,ModSymmbol)
for i=1:length(In)
    Dis = abs(ModSymmbol-In(i));
    [~,Ind] = min(Dis);
    Out(i) = Ind-1;
end

end










