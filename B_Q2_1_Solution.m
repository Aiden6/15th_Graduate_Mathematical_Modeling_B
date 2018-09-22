clear;clc;close all;
%�Ŵ���������
%Link_Num = 16;
%P_cross = 0.3; %�������0.3-0.9
%P_mutate = 0.01; %�������0.01-0.2

%Link_Num = 33;
%P_cross = 0.3; %�������0.3-0.9
%P_mutate = 0.02; %�������0.01-0.2

%Initial_population_individual_num = 400; %��ʼ��Ⱥ�ĸ�������
%Max_generation = 500;  %������

%Ŀ��������
Link_Num = 16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������Ϣ
City_Num = 12;
%�ܿ�������
Link_Num_Total=(City_Num-1) * City_Num/2; 
City_Name = ["������","����&���", "֣��", "�Ϻ�", "�人",  "����&����", "����", "��³ľ��", "����", "�ɶ�", "����", "����"];
City_Population = [9.9526 34.8245 9.03 23.8043 10.12 23.3853 8.83 3.8 0.9025 16.0477 7.2131 8.56];
City_Locaiton=[45.8 126.53;39.92 116.46;34.75 113.66;31.23 121.47;30.6 114.3;23.13 113.27;34.27 108.93;43.82 87.62;29.97 91.11;30.67 104.07;25.05 102.72;29.57 106.55];
%City_coordinateΪ��γ��ת���������
City_coordinate=city_information(City_Locaiton);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�Ŵ��㷨����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Initial_population_individual_num = 400; %��ʼ��Ⱥ�ĸ�������
Max_generation = 500;  %������
P_crossover = 0.3; %�������0.3-0.9
P_mutate = 0.01; %�������0.01-0.2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%������ʼ��Ⱥ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Initial_population = zeros(Initial_population_individual_num,Link_Num_Total);
for i=1:Initial_population_individual_num
    %����tsp���룬����ÿ����·���б���
    Initial_population(i,:) = randperm(Link_Num_Total);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%������Ⱥ����Ӧ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,P_popualtion_distri]=population_fitness(Initial_population,City_coordinate,City_Population,Link_Num);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%�����׶�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
generation=1;
Evolution_popualtion = Initial_population;
Generations_of_best_individual=zeros(Initial_population_individual_num,Link_Num_Total);
Population_cross_new=zeros(Initial_population_individual_num,Link_Num_Total);
Population_mutate_new=zeros(Initial_population_individual_num,Link_Num_Total);


while generation < Max_generation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%��������Ⱥ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:2:Initial_population_individual_num
        %������Ⱥ���ʷֲ�����ѡ�����
        selected_num=select(P_popualtion_distri);
        %����
        Population_cross=crossover(Evolution_popualtion,selected_num,P_crossover);  
        Population_cross_new(j,:)=Population_cross(1,:);
        Population_cross_new(j+1,:)=Population_cross(2,:);
        %����
        Population_mutate_new(j,:)=mutate(Population_cross_new(j,:),P_mutate);         
        Population_mutate_new(j+1,:)=mutate(Population_cross_new(j+1,:),P_mutate);
    end
    Evolution_popualtion=Population_mutate_new;  
    %%%%%%%%%%%%%%%%%%%%%%%%��������Ⱥ�и������Ӧ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Individual_fitness,P_popualtion_distri]=population_fitness(Evolution_popualtion,City_coordinate,City_Population,Link_Num);  

    %��¼������Ѹ������Ӧ�Ⱥ͸�����Ⱥ���е�λ��
    [Best_fitness_individual_value,Best_fitness_individual_num]=max(Individual_fitness);          %�ҳ���Ӧ�����ĸ���
    %��¼��ǰ������Ѹ���
    Best_individual=Evolution_popualtion(Best_fitness_individual_num,:);
    %��¼ÿ����Ѹ���
    Generations_of_best_individual(generation,:)=Best_individual;
    %��¼ÿ����Ѹ������Ӧ�ȣ����������ֵ
    Generations_of_best_individual_fitness(generation) = Single_Individual_fitness(Best_individual(end,1:Link_Num),City_coordinate,City_Population);
    
    generation=generation+1; %������һ��
    generation
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�����ܼ�ֵ���������ӣ���������
[Individual_fitness,City_connection,City_Capacity] = Single_Individual_fitness(Generations_of_best_individual(end,1:Link_Num),City_coordinate,City_Population);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ͼ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
clf;hold on;
%�����������
plot(City_Locaiton(:,2),City_Locaiton(:,1),'ko');
%��ʾ��������
for i= 1:City_Num
    text(City_Locaiton(i,2),City_Locaiton(i,1),City_Name{i},'fontsize',18)
end
%�����м�����
for ii=1:City_Num
    for jj=1:City_Num
        if (City_connection(ii,jj)==1)
            plot([City_Locaiton(ii,2) City_Locaiton(jj,2)],...
                [City_Locaiton(ii,1) City_Locaiton(jj,1)],'r--');
        end
    end
end

title(['����Ϊ��' num2str(Individual_fitness) 'mTb/s'])
ylabel('γ��','fontsize',15);
xlabel('����','fontsize',15);
figure(2);
plot(Generations_of_best_individual_fitness);
title(['������Ϊ��' num2str(Link_Num)]);
ylabel('�����ֵ','fontsize',15);
xlabel('��������','fontsize',15);

%pause(0.01);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%����������Ⱥ����Ӧ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����ÿ���������Ӧ�ȼ����䱻ѡ��ĸ���
%P_individual_for_liveΪÿ��������ĸ���
%����ֵ��ÿ����������ʵĸ��ʷֲ�
function [Individual_fitness,P_individual_for_live_distri]=population_fitness(All_individual,City_coordinate,City_Population,Link_Num)

Individual_num = size(All_individual,1);                           %���ĳ�����������    

Individual_fitness=zeros(Individual_num,1);                    %��ʼ���洢ĳ�����и��������ܼ�ֵ����   

for i=1:Individual_num
    Individual_fitness(i)=Single_Individual_fitness(All_individual(i,1:Link_Num),City_coordinate,City_Population);  %�����i����Ⱥ�����㺯��ֵ����ÿ����Ⱥ��Ӧ�ȣ���Ӧ��Խ��Խ���״��
end
Individual_fitness=Individual_fitness';  %ת��
%%%%��Ӧ���ܺ�
All_individual_fitness_sum=0;
for i=1:Individual_num
    All_individual_fitness_sum=All_individual_fitness_sum+Individual_fitness(i)^10;% ����Ӧ��Խ�õĸ��屻ѡ�����Խ��   ָ��Խ�ߣ�ǿ��Խ���״��
end
%P_individual_for_live��ÿ��������ĸ��ʣ���Ӧ��Խ�ߣ�������Խ��
P_individual_for_live=zeros(Individual_num,1);
for i=1:Individual_num
    P_individual_for_live(i)=Individual_fitness(i)^10/All_individual_fitness_sum;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%�����ۻ�����%�����ʷֲ�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_individual_for_live_distri=zeros(Individual_num,1);
P_individual_for_live_distri(1)=P_individual_for_live(1);
for i=2:Individual_num
    P_individual_for_live_distri(i)=P_individual_for_live_distri(i-1)+P_individual_for_live(i);
end
P_individual_for_live_distri=P_individual_for_live_distri';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%���㵥���������Ӧ�Ⱥ���%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������������ܼ�ֵ
function [Network_total_value,City_connection,City_Capacity]=Single_Individual_fitness(Single_individual_gene,City_coordinate,City_Population)
%����
%��66�����������ҵ�Link_Num��16��33�������Ӳ�����Ϊ���м����Ӿ���
City_Num = size(City_coordinate,1);  %�������%�����и���
City_connection = 1-tril(ones(City_Num,City_Num),0);     %ȡ�����Ǿ���
City_connection = City_connection(:);                      %���һ������
Ind = find(City_connection==1);              %�����Ǿ���Ԫ��Ϊ1λ������
City_connection(Ind)=1:length(Ind);          %�������ǵ���
for ii=1:length(Single_individual_gene)
    City_connection(City_connection==(Single_individual_gene(ii)))=Inf;
end
City_connection=(City_connection==Inf);
City_connection = reshape(City_connection,City_Num,City_Num);         
City_connection = City_connection + City_connection';

%�ж��Ƿ����������н��
if (length(TSP(1,City_connection)) < City_Num)
    Network_total_value =0;
    return;
end
%�����ܵ������ֵ��ȨֵΪ1
Network_total_value = 0;
Max_trans_dist_Req = [Inf 3000 1200 600];
Capacity_Table = [0 8 16 32];
City_Capacity = zeros(City_Num,City_Num);
for ii =1:City_Num
    for jj=ii+1:City_Num
        if (City_connection(ii,jj)>=1)          %���������
           Dist= norm(City_coordinate(ii,:) - City_coordinate(jj,:));   %�����м����
           Ind = find(Max_trans_dist_Req>Dist);
           City_Capacity(ii,jj) = Capacity_Table(Ind(end));%ѡ��ĳ���е���������
           Network_total_value = Network_total_value + sqrt(City_Population(ii)*City_Population(jj))*City_Capacity(ii,jj);
           
        end        
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�������ɸ���
%���ݽ�������ж��Ƿ񽻲�
%���ݱ�������ж��Ƿ����

function P=pro(p)
table(1:100)=0;
p_100=round(100*p);
table(1:p_100)=1;
n=round(rand*99)+1;
P=table(n);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ѡ��
function selected_num=select(p)

selected_num=zeros(2,1);
%����Ⱥ��ѡ��������ͬ�ĸ���
%���ѡ��һ������
for i=1:2
    r=rand;  %����һ��0-1�����
    prand=p-r;
    j=1;
    while prand(j)<0
        j=j+1;
    end
    selected_num(i)=j; %ѡ�и�������
    %%����ͬ����ѡһ��
    if (i==2)&&(j==selected_num(i-1))    
        r=rand;  %����һ�������
        prand=p-r;
        j=1;
        while prand(j)<0
            j=j+1;
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����
function Population_crossover=crossover(Evolution_popualtion,selected_num,P_crossover)

bn=size(Evolution_popualtion,2);
P_crossover=pro(P_crossover);  %���ݽ�����ʾ����Ƿ���н��������1���ǣ�0���
Population_crossover(1,:)=Evolution_popualtion(selected_num(1),:);
Population_crossover(2,:)=Evolution_popualtion(selected_num(2),:);
if P_crossover==1
    c1=round(rand*(bn-2))+1;  %��[1,bn-1]��Χ���������һ������λ
    c2=round(rand*(bn-2))+1;
    chb1=min(c1,c2);
    chb2=max(c1,c2);
    middle=Population_crossover(1,chb1+1:chb2);
    Population_crossover(1,chb1+1:chb2)=Population_crossover(2,chb1+1:chb2);
    Population_crossover(2,chb1+1:chb2)=middle;
    for i=1:chb1 
        while find(Population_crossover(1,chb1+1:chb2)==Population_crossover(1,i))
            zhi=find(Population_crossover(1,chb1+1:chb2)==Population_crossover(1,i));
            y=Population_crossover(2,chb1+zhi);
            Population_crossover(1,i)=y;
        end
        while find(Population_crossover(2,chb1+1:chb2)==Population_crossover(2,i))
            zhi=find(Population_crossover(2,chb1+1:chb2)==Population_crossover(2,i));
            y=Population_crossover(1,chb1+zhi);
            Population_crossover(2,i)=y;
        end
    end
    
    for i=chb2+1:bn
        while find(Population_crossover(1,1:chb2)==Population_crossover(1,i))
            zhi=logical(Population_crossover(1,1:chb2)==Population_crossover(1,i));
            y=Population_crossover(2,zhi);
            Population_crossover(1,i)=y;
        end
        while find(Population_crossover(2,1:chb2)==Population_crossover(2,i))
            zhi=logical(Population_crossover(2,1:chb2)==Population_crossover(2,i));
            y=Population_crossover(1,zhi);
            Population_crossover(2,i)=y;
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����
function snnew=mutate(snew,pm)
bn=size(snew,2);
snnew=snew;
pmm=pro(pm);  %���ݱ�����ʾ����Ƿ���б��������1���ǣ�0���
if pmm==1
    c1=round(rand*(bn-2))+1;  %��[1,bn-1]��Χ���������һ������λ
    c2=round(rand*(bn-2))+1;
    chb1=min(c1,c2);
    chb2=max(c1,c2);
    x=snew(chb1+1:chb2);
    snnew(chb1+1:chb2)=fliplr(x);
end
end
% map coordinates
function [City_coordinate]=city_information(City_Locaiton)
axesm utm   %����ͶӰ��ʽ
dczone=utmzone(City_Locaiton);
setm(gca,'zone',dczone)
h = getm(gca);
for i=1:length(City_Locaiton)
    [x,y]= mfwdtran(h,City_Locaiton(i,1),City_Locaiton(i,2));
    City_coordinate(i,:)=[x;y]/1e3;         %���km
end
end


function result=TSP(start_Node,Graph)
[m n]=size(Graph);
nodelist=zeros(m,1);
queue=start_Node;
nodelist(start_Node)=1;
result=start_Node;
while isempty(queue)==false
    i=queue(1);
    queue(1)=[];
    for j=1:n
        if((Graph(i,j)>0)&&(nodelist(j)==0)&&(i~=j))
            queue=[queue;j];
            nodelist(j)=1;
            result=[result;j];
        end
    end
end
end






