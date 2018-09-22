clear;clc;close all;
%遗传函数参数
%Link_Num = 16;
%P_cross = 0.3; %交叉概率0.3-0.9
%P_mutate = 0.01; %变异概率0.01-0.2

%Link_Num = 33;
%P_cross = 0.3; %交叉概率0.3-0.9
%P_mutate = 0.02; %变异概率0.01-0.2

%Initial_population_individual_num = 400; %初始种群的个体数量
%Max_generation = 500;  %最大代数

%目标连接数
Link_Num = 16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%城市信息
City_Num = 12;
%总可连接数
Link_Num_Total=(City_Num-1) * City_Num/2; 
City_Name = ["哈尔滨","北京&天津", "郑州", "上海", "武汉",  "广州&深圳", "西安", "乌鲁木齐", "拉萨", "成都", "昆明", "重庆"];
City_Population = [9.9526 34.8245 9.03 23.8043 10.12 23.3853 8.83 3.8 0.9025 16.0477 7.2131 8.56];
City_Locaiton=[45.8 126.53;39.92 116.46;34.75 113.66;31.23 121.47;30.6 114.3;23.13 113.27;34.27 108.93;43.82 87.62;29.97 91.11;30.67 104.07;25.05 102.72;29.57 106.55];
%City_coordinate为经纬度转换后的坐标
City_coordinate=city_information(City_Locaiton);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%遗传算法参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Initial_population_individual_num = 400; %初始种群的个体数量
Max_generation = 500;  %最大代数
P_crossover = 0.3; %交叉概率0.3-0.9
P_mutate = 0.01; %变异概率0.01-0.2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%产生初始种群%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Initial_population = zeros(Initial_population_individual_num,Link_Num_Total);
for i=1:Initial_population_individual_num
    %采用tsp编码，即对每条线路进行编码
    Initial_population(i,:) = randperm(Link_Num_Total);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%计算种群的适应度%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,P_popualtion_distri]=population_fitness(Initial_population,City_coordinate,City_Population,Link_Num);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%进化阶段%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
generation=1;
Evolution_popualtion = Initial_population;
Generations_of_best_individual=zeros(Initial_population_individual_num,Link_Num_Total);
Population_cross_new=zeros(Initial_population_individual_num,Link_Num_Total);
Population_mutate_new=zeros(Initial_population_individual_num,Link_Num_Total);


while generation < Max_generation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%产生新种群%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:2:Initial_population_individual_num
        %根据种群概率分布进行选择操作
        selected_num=select(P_popualtion_distri);
        %交叉
        Population_cross=crossover(Evolution_popualtion,selected_num,P_crossover);  
        Population_cross_new(j,:)=Population_cross(1,:);
        Population_cross_new(j+1,:)=Population_cross(2,:);
        %变异
        Population_mutate_new(j,:)=mutate(Population_cross_new(j,:),P_mutate);         
        Population_mutate_new(j+1,:)=mutate(Population_cross_new(j+1,:),P_mutate);
    end
    Evolution_popualtion=Population_mutate_new;  
    %%%%%%%%%%%%%%%%%%%%%%%%计算新种群中个体的适应度%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Individual_fitness,P_popualtion_distri]=population_fitness(Evolution_popualtion,City_coordinate,City_Population,Link_Num);  

    %记录当代最佳个体的适应度和个体在群落中的位置
    [Best_fitness_individual_value,Best_fitness_individual_num]=max(Individual_fitness);          %找出适应度最大的个体
    %记录当前代的最佳个体
    Best_individual=Evolution_popualtion(Best_fitness_individual_num,:);
    %记录每代最佳个体
    Generations_of_best_individual(generation,:)=Best_individual;
    %记录每代最佳个体的适应度，即总网络价值
    Generations_of_best_individual_fitness(generation) = Single_Individual_fitness(Best_individual(end,1:Link_Num),City_coordinate,City_Population);
    
    generation=generation+1; %进入下一代
    generation
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%表达该网络%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算总价值，城市连接，城市容量
[Individual_fitness,City_connection,City_Capacity] = Single_Individual_fitness(Generations_of_best_individual(end,1:Link_Num),City_coordinate,City_Population);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%作图%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
clf;hold on;
%画城市坐标点
plot(City_Locaiton(:,2),City_Locaiton(:,1),'ko');
%显示城市名字
for i= 1:City_Num
    text(City_Locaiton(i,2),City_Locaiton(i,1),City_Name{i},'fontsize',18)
end
%画城市间连线
for ii=1:City_Num
    for jj=1:City_Num
        if (City_connection(ii,jj)==1)
            plot([City_Locaiton(ii,2) City_Locaiton(jj,2)],...
                [City_Locaiton(ii,1) City_Locaiton(jj,1)],'r--');
        end
    end
end

title(['容量为：' num2str(Individual_fitness) 'mTb/s'])
ylabel('纬度','fontsize',15);
xlabel('经度','fontsize',15);
figure(2);
plot(Generations_of_best_individual_fitness);
title(['连接数为：' num2str(Link_Num)]);
ylabel('网络价值','fontsize',15);
xlabel('迭代次数','fontsize',15);

%pause(0.01);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%计算所有种群的适应度%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%根据每个个体的适应度计算其被选择的概率
%P_individual_for_live为每个个体存活的概率
%返回值：每个个体存活概率的概率分布
function [Individual_fitness,P_individual_for_live_distri]=population_fitness(All_individual,City_coordinate,City_Population,Link_Num)

Individual_num = size(All_individual,1);                           %获得某代个体的数量    

Individual_fitness=zeros(Individual_num,1);                    %初始化存储某代所有个体网络总价值数据   

for i=1:Individual_num
    Individual_fitness(i)=Single_Individual_fitness(All_individual(i,1:Link_Num),City_coordinate,City_Population);  %输入第i个种群，计算函数值，即每个种群适应度，适应度越大越容易存活
end
Individual_fitness=Individual_fitness';  %转置
%%%%适应度总和
All_individual_fitness_sum=0;
for i=1:Individual_num
    All_individual_fitness_sum=All_individual_fitness_sum+Individual_fitness(i)^10;% 让适应度越好的个体被选择概率越高   指数越高，强者越容易存活
end
%P_individual_for_live放每个个体存活的概率，适应度越高，存活概率越大
P_individual_for_live=zeros(Individual_num,1);
for i=1:Individual_num
    P_individual_for_live(i)=Individual_fitness(i)^10/All_individual_fitness_sum;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%计算累积概率%即概率分布%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_individual_for_live_distri=zeros(Individual_num,1);
P_individual_for_live_distri(1)=P_individual_for_live(1);
for i=2:Individual_num
    P_individual_for_live_distri(i)=P_individual_for_live_distri(i-1)+P_individual_for_live(i);
end
P_individual_for_live_distri=P_individual_for_live_distri';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%计算单个个体的适应度函数%%%%%%%%%%%%%%%%%%%%%%%%%%%
%即计算网络的总价值
function [Network_total_value,City_connection,City_Capacity]=Single_Individual_fitness(Single_individual_gene,City_coordinate,City_Population)
%解码
%从66条连接里面找到Link_Num（16或33）个连接并解码为城市间连接矩阵
City_Num = size(City_coordinate,1);  %坐标个数%即城市个数
City_connection = 1-tril(ones(City_Num,City_Num),0);     %取下三角矩阵
City_connection = City_connection(:);                      %组成一个向量
Ind = find(City_connection==1);              %下三角矩阵元素为1位置索引
City_connection(Ind)=1:length(Ind);          %右上三角递增
for ii=1:length(Single_individual_gene)
    City_connection(City_connection==(Single_individual_gene(ii)))=Inf;
end
City_connection=(City_connection==Inf);
City_connection = reshape(City_connection,City_Num,City_Num);         
City_connection = City_connection + City_connection';

%判断是否连接了所有结点
if (length(TSP(1,City_connection)) < City_Num)
    Network_total_value =0;
    return;
end
%计算总的网络价值，权值为1
Network_total_value = 0;
Max_trans_dist_Req = [Inf 3000 1200 600];
Capacity_Table = [0 8 16 32];
City_Capacity = zeros(City_Num,City_Num);
for ii =1:City_Num
    for jj=ii+1:City_Num
        if (City_connection(ii,jj)>=1)          %如果有连接
           Dist= norm(City_coordinate(ii,:) - City_coordinate(jj,:));   %两城市间距离
           Ind = find(Max_trans_dist_Req>Dist);
           City_Capacity(ii,jj) = Capacity_Table(Ind(end));%选择某城市的网络容量
           Network_total_value = Network_total_value + sqrt(City_Population(ii)*City_Population(jj))*City_Capacity(ii,jj);
           
        end        
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%按次生成概率
%根据交叉概率判断是否交叉
%根据变异概率判断是否变异

function P=pro(p)
table(1:100)=0;
p_100=round(100*p);
table(1:p_100)=1;
n=round(rand*99)+1;
P=table(n);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%选择
function selected_num=select(p)

selected_num=zeros(2,1);
%从种群中选择两个不同的个体
%随机选择一个个体
for i=1:2
    r=rand;  %产生一个0-1随机数
    prand=p-r;
    j=1;
    while prand(j)<0
        j=j+1;
    end
    selected_num(i)=j; %选中个体的序号
    %%若相同就再选一次
    if (i==2)&&(j==selected_num(i-1))    
        r=rand;  %产生一个随机数
        prand=p-r;
        j=1;
        while prand(j)<0
            j=j+1;
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%交叉
function Population_crossover=crossover(Evolution_popualtion,selected_num,P_crossover)

bn=size(Evolution_popualtion,2);
P_crossover=pro(P_crossover);  %根据交叉概率决定是否进行交叉操作，1则是，0则否
Population_crossover(1,:)=Evolution_popualtion(selected_num(1),:);
Population_crossover(2,:)=Evolution_popualtion(selected_num(2),:);
if P_crossover==1
    c1=round(rand*(bn-2))+1;  %在[1,bn-1]范围内随机产生一个交叉位
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
%变异
function snnew=mutate(snew,pm)
bn=size(snew,2);
snnew=snew;
pmm=pro(pm);  %根据变异概率决定是否进行变异操作，1则是，0则否
if pmm==1
    c1=round(rand*(bn-2))+1;  %在[1,bn-1]范围内随机产生一个变异位
    c2=round(rand*(bn-2))+1;
    chb1=min(c1,c2);
    chb2=max(c1,c2);
    x=snew(chb1+1:chb2);
    snnew(chb1+1:chb2)=fliplr(x);
end
end
% map coordinates
function [City_coordinate]=city_information(City_Locaiton)
axesm utm   %设置投影方式
dczone=utmzone(City_Locaiton);
setm(gca,'zone',dczone)
h = getm(gca);
for i=1:length(City_Locaiton)
    [x,y]= mfwdtran(h,City_Locaiton(i,1),City_Locaiton(i,2));
    City_coordinate(i,:)=[x;y]/1e3;         %变成km
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






