clear ;clc;close all;
%Ŀ��������
Link_Num = 16;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%������Ϣ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
City_Num = 12;
Link_Num_Total=(City_Num-1)*City_Num/2;

City_Name = ["������","����&���", "֣��", "�Ϻ�", "�人",  "����&����", "����", "��³ľ��", "����", "�ɶ�", "����", "����"];
City_Population = [9.9526 34.8245 9.03 23.8043 10.12 23.3853 8.83 3.8 0.9025 16.0477 7.2131 8.56];
City_Locaiton=[45.8 126.53;39.92 116.46;34.75 113.66;31.23 121.47;30.6 114.3;23.13 113.27;34.27 108.93;43.82 87.62;29.97 91.11;30.67 104.07;25.05 102.72;29.57 106.55];
%City_coordinateΪ��γ��ת���������
City_coordinate=city_information(City_Locaiton);
City_Dist = zeros(City_Num,City_Num);
for ii=1:City_Num
    for jj=1:City_Num
        City_Dist(ii,jj) = norm(City_Locaiton(ii,:)-City_Locaiton(jj,:));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�Ŵ��㷨����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Initial_population_individual_num=50; %��ʼ��Ⱥ��С
Max_generation=500;  %������
P_crossover=0.3; %�������
P_mutate=0.02; %�������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%������ʼ��Ⱥ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_Link_num = 2;
Code_Num = Link_Num_Total+Link_Num_Total+City_Num*City_Num;
Initial_population=zeros(Initial_population_individual_num,Code_Num);
for i=1:Initial_population_individual_num
    si=randperm(Link_Num_Total);
    sid=randi([1 max_Link_num],1,Link_Num_Total);
    xx = rand(City_Num,City_Num);
    xx = xx/sum(xx(:));
    Initial_population(i,:) = Code(si,sid,xx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%������Ⱥ����Ӧ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,P_popualtion_distri]=population_fitness(Initial_population,City_coordinate,City_Dist,City_Population,Link_Num);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%�����׶�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
generation=1;
Evolution_popualtion = Initial_population;
Generations_of_best_individual=zeros(Initial_population_individual_num,Code_Num);
Population_cross_new=zeros(Initial_population_individual_num,Code_Num);
Population_mutate_new=zeros(Initial_population_individual_num,Code_Num);
while generation<Max_generation+1
    for j=1:2:Initial_population_individual_num
        %ѡ��
        selected_num=select(P_popualtion_distri);  
        %����
        Population_cross=crossover(Evolution_popualtion,selected_num,P_crossover,City_Num,Link_Num); 
        Population_cross_new(j,:)=Population_cross(1,:);
        Population_cross_new(j+1,:)=Population_cross(2,:);
        %����
        Population_mutate_new(j,:)=mutate(Population_cross_new(j,:),P_mutate,City_Num,Link_Num);  
        Population_mutate_new(j+1,:)=mutate(Population_cross_new(j+1,:),P_mutate,City_Num,Link_Num);
    end
    Evolution_popualtion=Population_mutate_new;  %�������µ���Ⱥ
    %%%%%%%%%%%%%%%%%%%%%%%%��������Ⱥ�и������Ӧ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Individual_fitness,P_popualtion_distri]=population_fitness(Evolution_popualtion,City_coordinate,City_Dist,City_Population,Link_Num);
    %��¼������Ѹ������Ӧ�Ⱥ͸�����Ⱥ���е�λ��
    [Best_fitness_individual_value,Best_fitness_individual_num]=max(Individual_fitness);
    %��¼��ǰ������Ѹ���
    Best_individual=Evolution_popualtion(Best_fitness_individual_num,:);
    %��¼ÿ����Ѹ���
    Generations_of_best_individual(generation,:)=Best_individual;
    
    [~,~,XFp,si2,sid2] = Decode(Best_individual,12,Link_Num);
    XId = [si2;sid2]';
    %��¼ÿ����Ѹ������Ӧ�ȣ����������ֵ
    Generations_of_best_individual_fitness(generation)=Single_Individual_fitness(XId,XFp,City_Dist,City_coordinate,City_Population);
    generation=generation+1;
    generation
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�����ܼ�ֵ���������ӣ���������
[~,~,XFp,si2,sid2] = Decode(Generations_of_best_individual(end,:),12,Link_Num);
XId = [si2;sid2]';
[Individual_fitness,City_connection,City_capacity]=Single_Individual_fitness(XId,XFp,City_Dist,City_coordinate,City_Population);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ͼ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
clf;hold on;
plot(City_Locaiton(:,2),City_Locaiton(:,1),'ko');
%��ʾ��������
for i= 1:City_Num
    text(City_Locaiton(i,2),City_Locaiton(i,1),City_Name{i},'fontsize',18)
end
%�����м�����
for ii=1:City_Num
    for jj=1:City_Num
        if (City_connection(ii,jj)>=1)
            plot([City_Locaiton(ii,2) City_Locaiton(jj,2)],...
                [City_Locaiton(ii,1) City_Locaiton(jj,1)],'r--','linewidth',double(City_connection(ii,jj)));
        end
    end
end
title(['����Ϊ' num2str(Individual_fitness) 'mTb/s'])
ylabel('γ��','fontsize',15);
xlabel('����','fontsize',15);
figure(2);
plot(Generations_of_best_individual_fitness);
title(['������Ϊ' num2str(Link_Num)]);
ylabel('�����ֵ','fontsize',15);
xlabel('��������','fontsize',15);



%����������Ⱥ����Ӧ��
function [Individual_fitness,P_popualtion_distri]=population_fitness(Evolution_popualtion,City_coordinate,City_Dist,City_Population,Link_Num)

City_num = size(City_coordinate,1);
Individual_num=size(Evolution_popualtion,1);  %��ȡ��Ⱥ��С
Individual_fitness=zeros(Individual_num,1);

for i=1:Individual_num
    
    [~,~,XFp,si2,sid2] = Decode(Evolution_popualtion(i,:),City_num,Link_Num);
    XId = [si2;sid2]';
    Individual_fitness(i)=Single_Individual_fitness(XId,XFp,City_Dist,City_coordinate,City_Population);  %���㺯��ֵ������Ӧ��
end
Individual_fitness=Individual_fitness'; 
%���ݸ������Ӧ�ȼ����䱻ѡ��ĸ���
All_individual_fitness_sum=0;
for i=1:Individual_num
    All_individual_fitness_sum=All_individual_fitness_sum+Individual_fitness(i)^15;% ����Ӧ��Խ�õĸ��屻ѡ�����Խ��
end
P_individual_for_live=zeros(Individual_num,1);
for i=1:Individual_num
    P_individual_for_live(i)=Individual_fitness(i)^15/All_individual_fitness_sum;
end

%�����ۻ�����
P_popualtion_distri=zeros(Individual_num,1);
P_popualtion_distri(1)=P_individual_for_live(1);
for i=2:Individual_num
    P_popualtion_distri(i)=P_popualtion_distri(i-1)+P_individual_for_live(i);
end
P_popualtion_distri=P_popualtion_distri';
end

%����
%���ݱ�������ж��Ƿ����
function P=pro(p)
table(1:100)=0;
a=round(100*p);
table(1:a)=1;
b=round(rand*99)+1;
P=table(b);
end

%--------------------------------------------------
%��ѡ�񡱲���
function selected_num=select(p)
selected_num=zeros(2,1);
%����Ⱥ��ѡ���������壬��ò�Ҫ����ѡ��ͬһ������
for i=1:2
    r=rand;  %����һ�������
    prand=p-r;
    j=1;
    while prand(j)<0
        j=j+1;
    end
    selected_num(i)=j; %ѡ�и�������
    if i==2&&j==selected_num(i-1)    %%����ͬ����ѡһ��
        r=rand;  %����һ�������
        prand=p-r;
        j=1;
        while prand(j)<0
            j=j+1;
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scro=crossover(s,seln,pc,citenum,linenum)
pcc=pro(pc);  %���ݽ�����ʾ����Ƿ���н��������1���ǣ�0���

s1 = s(seln(1),:);
s2 = s(seln(2),:);
scro(1,:) = s1;
scro(2,:) = s2;
if pcc==1
    [si1,sid1,X1,~,~] = Decode(s1,citenum,linenum);
    [si2,sid2,X2,~,~] = Decode(s2,citenum,linenum);
    sbn = length(si1);
    
    c1=round(rand*(sbn-2))+1;  %��[1,bn-1]��Χ���������һ������λ
    c2=round(rand*(sbn-2))+1;
    chb1=min(c1,c2);
    chb2=max(c1,c2);
    
    %�����ֱַ���н���
    middle=si1(chb1+1:chb2);
    si1(chb1+1:chb2)=si2(chb1+1:chb2);
    si2(chb1+1:chb2)=middle;
    middle=sid1(chb1+1:chb2);
    sid1(chb1+1:chb2)=sid2(chb1+1:chb2);
    sid2(chb1+1:chb2)=middle;
    
    cc1=round(rand*(citenum-2))+1;  %��[1,bn-1]��Χ���������һ������λ
    cc2=round(rand*(citenum-2))+1;
    cchb1=min(cc1,cc2);
    cchb2=max(cc1,cc2);
    middle=X1(cchb1+1:cchb2,:);
    X1(cchb1+1:cchb2,:)=X2(cchb1+1:cchb2,:);
    X1(cchb1+1:cchb2,:)=middle;
    
    for i=1:chb1 %�ƺ�������
        while find(si1(chb1+1:chb2)==si1(i))
            zhi=find(si1(chb1+1:chb2)==si1(i));
            si1(i)=si2(chb1+zhi);
            sid1(i)=sid2(chb1+zhi);
        end
        while find(si2(chb1+1:chb2)==si2(i))
            zhi=find(si2(chb1+1:chb2)==si2(i));
            si2(i)=si1(chb1+zhi);
            sid2(i)=sid1(chb1+zhi);
        end
    end
    for i=chb2+1:sbn
        while find(si1(1:chb2)==si1(i))
            zhi=logical(si1(1:chb2)==si1(i));
            si1(i)=si2(zhi);
            sid1(i)=sid2(zhi);
        end
        while find(si2(1:chb2)==si2(i))
            zhi=find(si2(1:chb2)==si2(i));
            si2(i)=si1(zhi);
            sid2(i)=sid1(zhi);
        end
    end
    XMask = 1-tril(ones(12,12),1);
    X1 = X1.*XMask;X1= X1+X1';
    X2 = X2.*XMask;X2= X2+X2';
    
    scro(1,:) = Code(si1,sid1,X1);
    scro(2,:) = Code(si2,sid2,X2);
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function snnew=mutate(snew,pm,citenum,linenum)

snnew=snew;

pmm=pro(pm);  %���ݱ�����ʾ����Ƿ���б��������1���ǣ�0���
if pmm==1
    [si1,sid1,X1,~,~] = Decode(snew,citenum,linenum);
    bn=size(si1,2);
    
    c1=round(rand*(bn-2))+1;  %��[1,bn-1]��Χ���������һ������λ
    c2=round(rand*(bn-2))+1;
    chb1=min(c1,c2);
    chb2=max(c1,c2);
    x=si1(chb1+1:chb2);
    si1(chb1+1:chb2)=fliplr(x);
    x=sid1(chb1+1:chb2);
    sid1(chb1+1:chb2)=fliplr(x);
    
    c1=round(rand*(citenum-2))+1;  %��[1,bn-1]��Χ���������һ������λ
    c2=round(rand*(citenum-2))+1;
    chb1=min(c1,c2); chb2=max(c1,c2);
    x=X1(chb1+1:chb2,:);
    X1(chb1+1:chb2,:)=fliplr(x')';
    
    bit = randperm(bn,round(rand*bn));
    randbit = randi([-1 1],1,length(bit));
    sid1(bit)=sid1(bit)+randbit;
    sid1 = min(max(sid1,1),4);
    
    bit = randperm(citenum,round(rand*citenum));
    randval = randn(1,length(bit),citenum)*0.1;
    X1(bit,:)=randval;X1 = max(X1,1);
    XMask = 1-tril(ones(citenum,citenum),1);
    X1 = X1.*XMask;X1= X1+X1';
    X1 = X1/sum(X1(:));
    
    snnew = Code(si1,sid1,X1);
    
end
end

function [Network_total_value,City_connection,Xload2]=Single_Individual_fitness(XId,XFp,City_coordinate,City_Locaiton,City_Population)
City_num = size(City_Locaiton,1);

%�������Ӿ���
City_connection = 1-tril(ones(City_num,City_num),0);
City_connection = City_connection(:);
Ind = find(City_connection==1);
City_connection(Ind)=1:length(Ind);
Xc = zeros(1,City_num*City_num);

for ii=1:length(XId)
    Xc(City_connection==XId(ii,1)) = 1;
    City_connection(City_connection==XId(ii,1))=XId(ii,2);    
    
end
City_connection(Xc==0)=0;
City_connection = reshape(City_connection,City_num,City_num);
City_connection = City_connection+City_connection';
%�ж��Ƿ����������н�������
Xd = City_connection>0;
result=TSP(1,Xd);
if (length(result)<City_num)
    Network_total_value =0;
    return;
end

XInt = City_connection>0;
City_coordinate= City_coordinate.*XInt;
City_coordinate(XInt==0) = Inf;

%·�����ʼ���
Xload = ones(City_num,City_num);
SnrReq = [Inf 3000 1200 600];
Capacity_Table = [0 8 16 32];
for ii =1:City_num
    for jj=1:City_num
        if (City_connection(ii,jj)>=1)
            Dist= norm(City_Locaiton(ii,:)-City_Locaiton(jj,:));
            Ind = find(SnrReq>Dist);
            ModRate = Capacity_Table(Ind(end));
            Xload(ii,jj) = City_connection(ii,jj)*ModRate;
        end
    end
end

%���㺬һ���м�������µ����ʷ���
Xload2 = zeros(City_num,City_num);
XFp = XFp*sum(Xload(:))/4;
for ii =1:City_num
    for jj=ii+1:City_num
        [shortestPath, totalCost] = Dijkstra(City_coordinate, ii, jj);
        
        if (totalCost==Inf || length(shortestPath)~=3)
            Xload2(ii,jj) = 0;
            Xload2(jj,ii) = 0;
        else
            Xload2(jj,ii) = XFp(ii,jj);
            Xload2(jj,ii) = XFp(ii,jj);           
            
            Xload(ii,shortestPath(2))=Xload(ii,shortestPath(2))-XFp(ii,jj);
            Xload(shortestPath(2),ii)=Xload(shortestPath(2),ii)-XFp(ii,jj);
            Xload(shortestPath(2),jj)=Xload(shortestPath(2),jj)-XFp(ii,jj);
            Xload(jj,shortestPath(2))=Xload(jj,shortestPath(2))-XFp(ii,jj);
        end
    end
end


if (min(Xload(:))<0)
    Network_total_value= 0;
    return;
else
    for ii =1:City_num
        for jj=ii+1:City_num
            Xload2(ii,jj)= Xload(ii,jj);
            Xload2(jj,ii)= Xload(jj,ii);
        end
    end
end

% �����ܵļ�ֵ
Network_total_value = 0;
for ii =1:City_num
    for jj=ii+1:City_num
        if (City_connection(ii,jj)>=1)
            Dist= norm(City_Locaiton(ii,:)-City_Locaiton(jj,:));
            Ind = find(SnrReq>Dist);
            ModRate = Capacity_Table(Ind(end));
            Network_total_value = Network_total_value+ sqrt(City_Population(ii)*City_Population(jj))*Xload2(ii,jj);
        end
    end
end


end

function [shortestPath, totalCost] = Dijkstra(netCostMatrix, s, d)

n = size(netCostMatrix,1);
for i = 1:n
    % initialize the farthest node to be itself;
    farthestPrevHop(i) = i; % used to compute the RTS/CTS range;
    farthestNextHop(i) = i;
end
% all the nodes are un-visited;
visited(1:n) = false;
distance(1:n) = inf;    % it stores the shortest distance between each node and the source node;
parent(1:n) = 0;

distance(s) = 0;
for i = 1:(n-1)
    temp = [];
    for h = 1:n
        if ~visited(h)  % in the tree;
            temp=[temp distance(h)];
        else
            temp=[temp inf];
        end
    end
    [~, u] = min(temp);      % it starts from node with the shortest distance to the source;
    visited(u) = true;         % mark it as visited;
    for v = 1:n                % for each neighbors of node u;
        if ( ( netCostMatrix(u, v) + distance(u)) < distance(v) )
            distance(v) = distance(u) + netCostMatrix(u, v);   % update the shortest distance when a shorter shortestPath is found;
            parent(v) = u;     % update its parent;
        end
    end
end

shortestPath = [];
if parent(d) ~= 0   % if there is a shortestPath!
    t = d;
    shortestPath = [d];
    while t ~= s
        p = parent(t);
        shortestPath = [p shortestPath];
        
        if netCostMatrix(t, farthestPrevHop(t)) < netCostMatrix(t, p)
            farthestPrevHop(t) = p;
        end
        if netCostMatrix(p, farthestNextHop(p)) < netCostMatrix(p, t)
            farthestNextHop(p) = t;
        end
        
        t = p;
    end
end

totalCost = distance(d);
end


function result=TSP(startNode,Graph)

% �����������
% Graph ͼ��ͨ����
[m n]=size(Graph);
nodelist=zeros(m,1);
queue=startNode;
nodelist(startNode)=1;
result=startNode;
while isempty(queue)==false
    i=queue(1);
    queue(1)=[];
    for j=1:n
        if(Graph(i,j)>0&&nodelist(j)==0&&i~=j)
            queue=[queue;j];
            nodelist(j)=1;
            result=[result;j];
        end
    end
end
end
function y = Code(si,sid,X)
y = [si sid X(:)'];
end


function [si,sid,X,si2,sid2] = Decode(s,cn,linenum)

City_num = cn*(cn-1)/2;
si = s(1:City_num);
sid = s(City_num+1:2*City_num);
X = s(2*City_num+1:end);
X = reshape(X,cn,cn);

kk = cumsum(sid);
[~,Ind] =find(kk>=linenum,1);
si2 = si(1:Ind);
kk(Ind) = linenum;
sid2 = diff([0 kk(1:Ind)]);
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
