纯量
name <- "XiaoMing"; L <- TRUE
1+2*3+4/5+6^2
sqrt(2)
log(10)
3==5
x <- 3
x<-"jjjjj"
y <- as.character(x); y
is.numeric (y)
向量
x <- c(10.4,5.6,3.1,6.4,21.7)
y <- c(x,0,x)
1:10
10:1
y <- c(8,3,5,7,6,2,8,9);y > 5
all (y > 5) 
any(y>5) 
which(y>5) 
seq(0,1,length.out=11)
seq(1,9,by=2)
seq(1,9,by=pi)
seq(1,6,by=3)
seq(10)
seq(0,1,along.with=rnorm(11))
rep(1:4,times=2);rep(1:4,2)
rep(1:4,length.out=10)
rep(1:4,each=2)
rep(1:4,c(1,2,2,3))
x<-c(1,4,7);x[2]
x
x[1:5]
x[c(1,2,3,2,1)]
c("a","b","c")[rep(c(2,1,3),times=3)]
x<-10:20;x[-(1:5)]
x <- c(1,4,7); x[x<5]
(ages <- c(Li=33,Zhang=29,Liu=18))
因子
data<-c(1,2,3,3,1,2,2,3,1,3,2,1)
(fdata<-factor(data))
(rdata<-factor(data,labels=c("I","II","III")))
gl(3,5,labels=paste0("A",1:3))
gl(5,1,length=15,labels=paste0("B",1:5))
矩阵
(mdata<-matrix(c(1,2,3,11,12,13),nrow=2,ncol=3,byrow=TRUE,dimnames = list(c("row1","row2"),c("C.1","C.2","C.3"))))
mdata
A<-matrix(1:15,nrow=3,ncol=5);A
A<-matrix(1:15,nrow=3)
A<-matrix(1:15,ncol=5)
B<-matrix(nr=2,nc=3)
B[1,1]<-1;B[1,3]<-0;B[2,2]<-3;B
X<-1:12;dim(X)<-c(3,4);X
X1<-rbind(1:2,101:102);X1
X2<-cbind(1:2,101:102);X2
cbind(X1,X2)
rbind(X1,X2)
dim(A)
nrow(A)
ncol(A)
as.vector(A)
A[1,2];A[1,2]<-102
A[c(1,3),2:4]
A[-3,-2];A[-1,];A[,-2]
A[c(1,3),];A[2,]<-201:205;A
数组
X<-array(1:20,dim=c(4,5));X
Y<-array(1:24,dim=c(3,4,2));Y
Y<-1:24
dim(Y)<-c(3,4,2)
a<-1:24
dim(a)<-c(2,3,4);a
a[2,1,2]
a[1,2:3,2:3]
a[1,,]
a[,2,]
a[1,1,]
a[]
a[]<-0;a
列表
Lst<-list(name="Fred",wife="Mary",no.children=3,child.ages=c(4,7,9))
Lst
Lst[[2]]
Lst[[4]][2]
Lst[["name"]]
Lst[["child.ages"]]
Lst$name;Lst$wife;Lst$child.ages
Lst$name<-"John"
Lst$income<-c(1980,1600)
Lst$income<-NULL
数据框
df<-data.frame(
     Name=c("Alice","Becka","James","Jeffrey","John"),
     Sex=c("F","F","M","M","M"),
     Age=c(13,13,12,13,12),
     Height=c(56.5,65.3,57.3,62.5,59.0),
     Weight=c(84.0,98.0,83.0,84.0,99.5)
);df
df[1:2,3:5]
df[["Height"]]
df$Weight
attach(df)
r<-Height/Weight;r
df$r<-Height/Weight
detach()
df$r<-with(df,Height/Weight)
读写数据文件
读写纯文本文件
rt<-read.table("houses.data");rt
is.data.frame(rt)
class(rt)
rt<-read.table("houses.data",header=TRUE);rt
w<-scan("weight.data");w
is.vector(w)
inp<-scan("h_w.data",list(height=0,weight=0));inp
x<-scan()
names<-scan(what="")
读写Excel表格数据
read.delim("educ_scores.txt")
read.csv("educ_scores.csv")
install.packages("RODBC")
library(RODBC)
con<-odbcConnectExcel("educ_scores.xls")
tabls<-sqlTables(con)
sh1<-sqlFetch(con,tbls$TABLE_NAME[1])
qry<-paste("selet*from[",tbls$TABLE_NAME[1],"]",sep="")
sh2<-sqlQuery(con,qry);sh2
数据集的读取
data()
data(infert)
str(infert)#查看
data(package="cluster")
data(agriculture,package="cluster")
library("cluster")
data()
data(agriculture)
str(agriculture)#
写数据文件
X<-matrix(1:12,ncol=6);X
write(X,file="Xdata.txt")
df<-data.frame(
  Name=c("Alice","Becka","James","Jeffrey","John"),
  Sex=c("F","F","M","M","M"),
  Age=c(13,13,12,13,13),
  Height=c(56.5,65.3,57.3,62.5,59.0),
  Weight=c(84.0,98.0,83.0,84.0,99.5)
)
write.table(df,file="foo.txt")
write.csv(df,file="foo.csv")
控制流
if(any(x<=0))y<-log(1+x) else y<-log(x)
y<-if(any(x<=0))log(1+x) else log(x)
n<-4;x<-array(0,dim=c(n,n))
for(i in 1:n){
  for(j in 1:n){
    x[i,j]<-1/(i+j-1)
  }
}
x
f<-c(1,1);i<-1
while(f[i]+f[i+1]<1000){
  f[i+2]<-f[i]+f[i+1]
  i<-i+1;
}
f
f<-c(1,1);i<-1
repeat{
  f[i+2]<-f[i]+f[i+1]
  i<-i+1
  if(f[i]+f[i+1]>=1000)break
}
f
R语言的程序设计
X<-scan("sample.data",nlines=2)
Y<-scan("sample.data",skip=2)
source("t.stat.R");t.stat(X,Y)
t.stat(x=X,y=Y);t.stat(y=Y,x=X)
source("moment.R");moment(X,k=2);moment(X,k=2,mean=mean(X))
source("fac1.R");fac1(4)
source("fac2.R");fac1(4)
prod(seq(4))
gamma(5)
factorial(4)
