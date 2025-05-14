plt.rcParams['font.sans-serif]=['STSong'] class data_m():
def       nit     (self,load_max,time_max,car_v,parm_p,cl1,c12,c21,c22,popsize=50): self.load_max   =load_max   #  载重
self.time_max   =time_max

self.car_v=car_v                      %速度
self.s1=parm_p[0]       %alpha
self.s2=parm_p[1]       %beta
self.pi=parm_p[2]%c 单位时间惩罚成本
self.cll=cl1
self.c12=c12
self.c21=c21
self.c22=c22
self.popsize=popsize
def          distance(self,lat1,lon1,lat2,lon2):#  纬度，经度
%将角度转换为弧度
d=((lat1-lat2)**2+(lon1-lon2)**2)**0.5
%计算结果单位为公里
return d
def information(self):
data   =pd.read_excel('customer.xlsx')
location_xy=data[['x','y']].values %x 为坐标，y 为坐标
demand=data['demand'].values
time_window=data[['ET','LT']]
service_time   =data['Service_time'].values
distance  =[]





for i in range(len(location_xy)):
dis  =[]

for jin  range(len(location_xy)):
d                  self.distance(location_xy[i][1],        location_xy[i][0],         location_xy[i][1], location_xy[j][0])
dis.append(d)
distance.append(dis)
return      np.array(location_xy),np.array(demand),np.array(time_window),np.array(service_time), np.array(
distance)
def   draw_initial(self,location_xy):
plt.figure()
plt.scatter(location_xy[0][0],location_xy[0][1],color='red')
for i in range(1,len(location_xy)):
plt.plot([location_xy[i][0],     location_xy[0][0]],      [location_xy[i][1],     location_xy[0][0]], color='black')
plt.scatter(location_xy[i][0],location_xy[i][1],color='blue')
plt.annotate('客户%.0f%'(i),xy=(location_xy[i][0],location_xy[i][1]),
xytext=(location_xy[i][0]-0.3,location_xy[i][1]+0.5),fontsize=10,
weight='heavy')
plt.show(block=True)
def  decode(self,road,distance,demand):
car_road,car_p,car_s,car_time=[],[],[,[
#car_road,car_p,car_s,car_time=[],[],[,[
time,car_code,car_demand,car_window,car_distance,car_un=[],[],[,[],[],[
signal=-1      %车辆标号
window_low=[]
cost11,cost12,cost21,cost22=0,0,0,0
total_distance=0
total_time=0

for i in range(len(road)):
forjin   range(len(road[i])-1):
d=distance[road[i][j]][road[i][j+1]]
time_transport=d/self.car_v
cost11+=time_transport*self.cl1
cost12+=demand[road[i][j+1]]*service_time[road[i][j+1]]

cost21+=d*self.c21
total_distance+=d
total_time+=time_transport+service_time[road[i][j+1]]
cost22+=self.c22

cost=cost11+cost12+cost21+cost22
car_num=len(road)
return cost,total_distance,total_time,car_num
def        distance_saving(self,location_xy,demand,time_window,service_time,distance,d_sort): customers=copy.deepcopy(demand).tolist()
car_road,car_p,car_s,car_time         =[],[],[],[]
time,car_code,car_demand,car_window,car_distance,car_un=[],[],[],[],[],[]
signal=-1%车辆标号
car_load_left=[]
d_saving=copy.deepcopy(d_sort)
A=False
AA=[i for i in range(1,28)]
while  len(d_saving[0])!=0:
for i in range(len(d_saving[0])):
ifi<1 or A:
ifi>0:
car_road[signal].append(0)
car_load=self.load_max
car_t=self.time_max
car_road.append([0])

signal+=1
A=False
if d_saving[1][i]not in car_road[signal]:
if car_load>demand[d_saving[1][i]]:
car_road[signal].append(d_saving[1][i])
customers.remove(d_saving[1][i])
car_load   -=demand[d_saving[1][i]]
elif   car_load>min(demand[list(set(d_saving[1]+d_saving[2]))]):
pass
else:
A=True
if d_saving[2][i]not in car_road[signal]:
if car_load>demand[d_saving[2][i]]:
car_road[signal].append(d_saving[2][i])
customers.remove(d_saving[2][i])
car_load  -=demand[d_saving[2][i]]
elif       car_load>min(demand[list(set(d_saving[1]+d_saving[2]))]):
pass
else:
A=True
ifA==True:
index=[]
forjin   range(len(d_saving[0])):
if   d_saving[1][j]    not    in   car_road[signal]    and      d_saving[2][j]not      in
car_road[signal]:
index.append(j)
re=[[],[,[]]
for jin range(len(index)):
re[0].append(d_saving[0][index[j]])
re[1].append(d_saving[1][index[j]])

re[2].append(d_saving[2][index[j]]) d_saving=copy.deepcopy(re)
break
car_road[signal].append(0)
if i==len(re[0])-1:
break
car_load_left.append(car_load) print(car_road)
print(car_load_left)
return car_road
def create_road(self.road):
roadl=copy.deepcopy(road)
for i in range(len(road)):
a,b=[],[0]
forjin    range(1,len(road[i])-1):
a.append(road[i][j])
np.random.shuffle(a)
a.append(0)
b=b+a
road1[i]=b
return  roadl
def    roulette(self,total_chrom,total_answer):
population=[]
fitness     =(1/np.array(total_answer)).tolist() sum_fitness   =sum(fitness)
pk,qk=[].[]
for  i  in  range(self.popsize*2):
pk.append(fitness[i]/sum_fitness) for  i  in  range(self.popsize*2):
cumulative  =0

forjin       range(0,i+1):
cumulative  =cumulative  +pk[j]
qk.append(cumulative)
selection_rand=[np.random.rand()for  i  in  range(self.popsize)] for i in range(self.popsize):
if   selection_rand[i]<=qk[0]:
population.append(copy.deepcopy(total_chrom[0])) else:
forj    in    range(0,self.popsize*2-1):
if       qk[j]<selection_rand[i]<=qk[j+1]:
population.append(copy.deepcopy(total_chrom[j+1])) break
return population
def GA(self,road_init,distance,demand):
population=[]
Total_road1=[]
answer=[]
fit_every  =[]
for gen in range(100):
if(gen<1):
for i in range(50):
road=self.create_road(road_init,)
cost,_,_,_=self.decode(road,distance,demand)
population.append(road)
answer.append(cost)
best_index=answer.index(min(answer))
fit_every.append(min(answer)
Tbest=min(answer)
Tbest_road=population[best_index]
print(Tbest)

answer1=[]
population1=[]
for  i  in range(0,50,2):
chroml=copy.deepcopy(population[i])
chrom2=copy.deepcopy(population[i+1])
chroml=self.create_road(chroml)
chrom2=self.create_road(chrom2)
cost,_,_,_=self.decode(chrom1,distance,demand)
answer1.append(cost)
cost,_,_,_=self.decode(chrom2,distance,demand)
answer1.append(cost)
population1.append(chrom1)
population1.append(chrom2)
total_population=copy.deepcopy(population)+copy.deepcopy(populationl) total_answer=answer+answer1
population         =self.roulette(total_population,total_answer) Tbest_now=min(total_answer)
if     Tbest_now<Tbest:
Tbest=Tbest_now

Tbest_road=copy.deepcopy(total_population[total_answer.index(min(total_answer)))
fit_every.append(Tbest)
print('当前迭代次数为%.0f,最小成本为%.3f,最优路径为%s'(gen,Tbest,Tbest_road))
return    Tbest_road,fit_every
def draw_change(self,fit_every):
x=[i  for  i  in  range(len(fit_every))] plt.plot(x,fit_every)
plt.xlabel ('迭代次数',fontsize=20) 
plt.ylabel ('总费用',fontsize=20)
plt.tick_params(labelsize=15)
plt.show()



def draw(self,car_road):
.plt.figure()
for i in range(len(car_road)):
x=[location_xy[i][0]for jin car_road[i]]
y=[location_xy[j][1]for jin car_road[i]
plt.scatter(x,y,c="red")
if(i<len(car_road)-1):
plt.plot(x,y)
else:
plt.plot(x,y)
for k in range(len(car_road[i])):
plt.annotate(car_road[i][k],xy=(x[k],y[k]),xytext=(x[k],y[k]),fontsize=20,weight='heavy') plt.legend(prop={'family':['STSong'],'size':16})#标签字体大小，可以修改 font1={'weight:'bold,'size':22}#汉字字体大小，可以修改
plt.xlabel('横坐标',font1)
plt.title('最优配送方案路线图',font1)
plt.ylabel('纵坐标',font1)
plt.axis()
#plt.plot(x,y,label='目标函数值为：%.2f%(Z))
plt.show()
load_max,time_max=9,8       #  载重  car_v=50  # 速度20/s 转化为km/h alpha,beta,c=1,2,10
parm_p=[alpha,beta,c]#         惩罚成本各参数，s1,s2  和违反约束比例 customers=28      #100个需求点
c11,c12,c21,c22=20,25,100,50
a=data_m(load_max,time_max,car_v,parm_p,cl1,c12,c21,c22)#             解码模块 location_xy,demand,time_window,service_time,distance=a.information()
a.draw_initial(location_xy)
#print(distance)


d_sort1.d_sort2,d_sort3=[],[],[] for i in range(2,len(distance)):
forjin     range(1,i):
d_sort1.append(distance[i][0]+distance[0][j]-distance[i][j])
d_sort2.append(i)
d_sort3.append(j)
print(d_sort1)
index =np.argsort(d_sort1)[::-1]# 降序排列 d_sort  =[]
d_sort.append(np.array(d_sort1)[index].tolist())
d_sort.append(np.array(d_sort2)[index].tolist())


for i in range(len(d_sor[0])):
print('客户%.0f与客户%.0f之间节约里程为%0.3f%(d_sor[1][i],d_sort[2][i],d_sort[0][i])) road_r=[[0,9,3,4,5,20,22,23,0],[0,2,1,8,10,21,24,28,0],[0,6,7,12,11,13,25,26,0],[0,14,15,16,17,18,19,2  7,0]]
r_cost,r_total_distance,r_total_time,r_car_num=a.decode(road_r,distance,demand)
print( '原始方案',road_r,  '原始方案成本：% . 3f, 原始方案总距离：% . 3f, 原始方案总时间：% . 3f, 原 始方案车辆数：% . 3f%(r_cost,r_total_distance,r_total_time,r_car_num))
a.draw(road_r)
car_road=a.distance_saving(location_xy,demand,time_window,service_time,distance,d_sort)
print(car_road)
cost,total_distance,total_time,car_num=a.decode(car_road,distance,demand)
print('初始方案’,car_road,  '初始方案成本：%.3f, 初始方案总距离：% . 3f,初始方案总时间：% . 3f,初 始方案总车辆数：% . 3f%(cost,total_distance,total_time,car_num))
a.draw(car_road)
Tbest_road,fit_every=a.GA(car_road,distance,demand)

Tbest_cost,Tbest_total_distance,Tbest_total_time,Tbest_car_num=a.decode(Tbest_road,distance,dema nd)
print(  '初始方案',Tbest_road,   '最优方案成本：%.3f, 最优方案总距离：% . 3f, 最优方案总时间：% . 3f,
68-
中国知网   https://www.cnki.net

最优方案总车辆数：%.3f%(Tbest_cost,Tbest_total_distance,Tbest_total_time,Tbest_car_num))
a.draw_change(fit_every)
a.draw(Tbest_road)