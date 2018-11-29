【包】
library(zoo)            #时间格式预处理
library(xts)            #同上
library(timeSeires)      #同上
library(urca)           #进行单位根检验
library(tseries)         #arma模型
library(fUnitRoots)     #进行单位根检验
library(FinTS)         #调用其中的自回归检验函数
library(fGarch)        #GARCH模型
library(nlme)          #调用其中的gls函数
library(fArma)        #进行拟合和检验

【基本函数】
数学函数
abs，sqrt：绝对值，平方根 log, log10, log2 , exp：对数与指数函数 sin，cos，tan，asin，acos，atan，atan2：三角函数 sinh，cosh，tanh，asinh，acosh，atanh：双曲函数 
简单统计量
sum, mean, var, sd, min, max, range, median, IQR（四分位间距）等为统计量，sort，order，rank与排序有关，其它还有ave，fivenum，mad，quantile，stem等。


【数据处理】
#具体说明见文档1
#转成时间序列类型
x = rnorm(2)
charvec = c(“2010-01-01”,”2010-02-01”)
zoo(x,as.Date(charvec))     #包zoo
xts(x, as.Date(charvec))     #包xts
timeSeries(x,as.Date(charvec))  #包timeSeries
#规则的时间序列，数据在规定的时间间隔内出现
tm = ts(x,start = c(2010,1), frequency=12 )  #12为按月份，4为按季度，1为按年度
zm = zooreg(x,start = c(2010,1), frequency=12 )  #包zoo
xm = as.xts(tm)     #包xts
sm = as.timeSeries(tm) #包timeSeries
#判断是否为规则时间序列
is.regular(x)

#排序
zoo()和xts()会强制变换为正序（按照时间名称）
timeSeries不会强制排序；其结果可以根据sort函数排序，也可以采用rev()函数进行逆序；参数recordIDs，可以给每个元素（行）标记一个ID，从而可以找回原来的顺序

#预设的时间有重复的时间点时
zoo会报错
xts按照升序排列
timeSeries把重复部分放置在尾部； 

#行合并和列合并
#都是按照列名进行合并，列名不同的部分用NA代替
cbind()
rbind()
merge() 列合并

#取子集
xts()默认将向量做成了矩阵；其他与常规向量或者矩阵没有差别

#缺失值处理
na.omit(x) 
x[is.na(x)] = 0
x[is.na(x)] = mean(x,na.rm=TRUE)
x[is.na(x)] = median(x,na.rm=TRUE)
na.approx(x)  #对缺失值进行线性插值
na.spline(x)   #对缺失值进行样条插值
na.locf(x)     #末次观测值结转法
na.trim(x, sides=”left” )  #去掉最后一个缺失值
#对timeSreies数据
na.omit(x, “ir” )  #去掉首末位置的缺失值
na.omit(x, “iz” )  #用替换首末位置的缺失值
na.omit(x, “ie” )  #对首末位置的缺失值进行插值
na.omit(x, method=“ie”, interp= c(“before”,”linear”,”after”) ) #可以选择插值方法，before末次观测值法，after下次观测结转法

as.contiguous(x)  #返回x中最长的连续无缺失值的序列片段，如果有两个等长的序列片段，则返回第一个。

#时间序列数据的显示
#zoo和xts都只能按照原来的格式显示，timeSeries可以设置显示格式
print(x, format= “%m/%d/%y %H:%M”)  #%m表示月，%d表示天，%y表示年，%H表示时，%M表示分钟，%A表示星期，%j表示天的序号
      #timeSeries也可以按照ts的格式显示
      print(x, style=”ts”)
      print(x, style=”ts”, by=”quarter”)
      
      【图形展示】
      plot.zoo(x)
      plot.xts(x)
      plot.zoo(x, plot.type=”single”) #支持多个时间序列数据在一个图中展示
      plot(x, plot.type=”single”) #支持多个时间序列数据在一个图中展示，仅对xts不行
      
      【基本统计运算】
      1、自相关系数、偏自相关系数等
      例题2.1 
      d=scan("sha.csv")
      sha=ts(d,start=1964,freq=1)
      plot.ts(sha)   #绘制时序图
      acf(sha,22)   #绘制自相关图，滞后期数22
      pacf(sha,22)  #绘制偏自相关图，滞后期数22
      corr=acf(sha,22)   #保存相关系数
      cov=acf(sha,22,type = "covariance")   #保存协方差
      
      2、同时绘制两组数据的时序图
      d=read.csv("double.csv",header=F)
      double=ts(d,start=1964,freq=1)
      plot(double, plot.type = "multiple")   #两组数据两个图
      plot(double, plot.type = "single")     #两组数据一个图 
      plot(double, plot.type = "single",col=c("red","green"),lty=c(1,2)) #设置每组数据图的颜色、曲线类型)
      
      3、纯随机性检验
      例题2.3续
      d=scan("temp.csv")
      temp=ts(d,freq=1,start=c(1949))
      Box.test(temp, type="Ljung-Box",lag=6)
      
      4、差分运算和滞后运算
      diff
      lag
      
      5、模拟ARIMA模型的结果
      arima.sim(n = 100, list(ar = 0.8))
      plot.ts(arima.sim(n = 100, list(ar = 0.8)))   #会随机产生一个包含100个随机数的时序图
      plot.ts(arima.sim(n = 100, list(ar = -1.1)))   #非平稳，无法得到时序图。
      plot.ts(arima.sim(n = 100, list(ar = c(1,-0.5))))
      plot.ts(arima.sim(n = 100, list(ar = c(1,0.5))))
      arima.sim(n = 1000, list(ar = 0.5, ma = -0.8))
      acf(arima.sim(n = 1000, list(ar = 0.5, ma = -0.8)),20)
      pacf(arima.sim(n = 1000, list(ar = 0.5, ma = -0.8)),20)
      
      【单位根检验】
      #方法1
      b=ts(read.csv("6_1.csv",header=T)) 
      x=b[,1]
      y=b[,1]
      summary(ur.df(x,type="trend",selectlags="AIC"))
      #方法2：单位根检验更好的函数，加了画图的功能 
      library(fUnitRoots)
      urdfTest(x)
      #方法3：ADF检验的一个自编函数
      library(urca)
      #...
      ur.df.01=function(x,lags=8){    
        #将三种ADF检验形式汇总的函数（结果和EVIEWS不一致）
        res=matrix(0,5,3)
        colnames(res)=c("无","含常数项","含常数项和趋势项")
        rownames(res)=c("tau统计量","1%临界值","5%临界值",
                        "10%临界值","是否稳定（1/0）")
        types=c("none","drift","trend")
        for(i in 1:3){
          x.adf=ur.df(x,type=types[i],lags=lags,selectlags="AIC")
          x.adf.1=x.adf@teststat  #统计量
          x.adf.2=x.adf@cval      #临界值
          res[1,i]  =x.adf.1[1]
          res[2:4,i]=x.adf.2[1,]
          res[5,i]=if( abs(res[1,i]) > abs(res[3,i]) ) 1 else 0
        }
        return(res)
      }
      #...
      ur.df.01(x)              #对原序列进行判断
      
      【一般的ARIMA模型】
      d=scan("a1.5.txt")               #导入数据
      prop=ts(d,start=1950,freq=1)      #转化为时间序列数据
      plot(prop)                   #作时序图
      acf(prop,12)                 #作自相关图，拖尾
      pacf(prop,12)                #作偏自相关图，1阶截尾
      Box.test(prop, type="Ljung-Box",lag=6)  
      #纯随机性检验,p值小于5%，序列为非白噪声
      Box.test(prop, type="Ljung-Box",lag=12)
      ( m1=arima(prop, order = c(1,0,0),method="ML") )    #用AR(1)模型拟合，如参数method="CSS"，估计方法为条件最小二乘法，用条件最小二乘法时，不显示AIC。
      ( m2=arima(prop, order = c(1,0,0),method="ML", include.mean = F) ) #用AR(1)模型拟合，不含截距项。
      tsdiag(m1)  #对估计进行诊断，判断残差是否为白噪声
      summary(m1)
      r=m1$residuals  #用r来保存残差
      Box.test(r,type="Ljung-Box",lag=6, fitdf=1)#对残差进行纯随机性检验，fitdf表示残差减少的自由度
      AutocorTest(m1$resid)                    #加载FinTS包，进行自相关检验
      prop.fore = predict(m1, n.ahead =5)  #将未来5期预测值保存在prop.fore变量中
      U = prop.fore$pred + 1.96* prop.fore$se  #会自动产生方差
      L = prop.fore$pred – 1.96* prop.fore$se   #算出95%置信区间
      ts.plot(prop, prop.fore$pred, col=1:2)      #作时序图，含预测。
      lines(U, col="blue", lty="dashed")
      lines(L, col="blue", lty="dashed")#在时序图中作出95%置信区间
      
      ——说明：运行命令arima(prop, order = c(1,0,0),method="ML")之后，显示：
      Call:
        arima(x = prop, order = c(1, 0, 0), method = "ML")
      Coefficients:
        ar1    intercept
      0.6914    81.5509
      s.e.   0.0989     1.7453
      sigma^2 estimated as 15.51:  log likelihood = -137.02,  aic = 280.05
      注 意：intercept下面的81.5509是均值，而不是截距！虽然intercept是截距的意思，这里如果用mean会更好。（the mean and the intercept are the same only when there is no AR term，均值和截距是相同的，只有在没有AR项的时候）
      如果想得到截距，利用公式计算。int=(1-0.6914)*81.5509= 25.16661。
      
      ——说明：Box.test(r,type="Ljung-Box",lag=6，fitdf=1)
      fitdf表示p+q，number of degrees of freedom to be subtracted if x is a series of residuals，当检验的序列是残差到时候，需要加上命令fitdf，表示减去的自由度。
      运行Box.test(r,type="Ljung-Box",lag=6,fitdf=1)后，显示的结果：
      Box.test(r,type="Ljung-Box",lag=6,fitdf=1)
      Box-Ljung test
      data:  r 
      X-squared = 5.8661, df = 5, p-value = 0.3195
      “df = 5”表示自由度为5，由于参数lag=6，所以是滞后6期的检验。
      
      #另一个参数估计与检验的方法（加载fArma程序包）
      ue=ts(scan("unemployment.txt"),start=1962,f=4) #读取数据
      due=diff(ue)
      ddue=diff(due,lag=4)
      fit2=armaFit(~arima(4,0,0),include.mean=F,data=ddue,method="ML")  #另一种拟合函数
      summary(fit2)
      fit3=armaFit(~arima(4,0,0),data=ddue,transform.pars=F,fixed=c(NA,0,0,NA),include.mean=F,method="CSS")
      summary(fit3)
      
      【一些特殊的模型】
      #固定某些系数的值
      arima(dw,order=c(4,0,0),fixed=c(NA,0,0,NA,0),method="CSS")
      
      #乘积季节模型
      wue=ts(scan("wue.txt"),start=1948,f=12)
      arima(wue,order=c(1,1,1),seasonal=list(order=c(0,1,1),period=12),include.mean=F,method="CSS")
      
      #拟合自回归模型，因变量关于时间的回归模型
      eg1=ts(scan("582.txt"))
      ts.plot(eg1)
      fit.gls=gls(eg1~-1+time(eg1), correlation=corARMA(p=1), method="ML") #看nlme包
      summary(fit.gls2) 
      #或
      fit=arima(eg1,c(1,0,0),xreg=time(eg1),include.mean=F,method="ML")
      AutocorTest(fit$resid)    #残差白噪声检验 
      
      #延迟因变量回归模型
      leg1=lag(eg1,-1)
      y=cbind(eg1,leg1)
      fit=arima(y[,1],c(0,0,0),xreg=y[,2],include.mean=F)
      
      #拟合GARCH模型
      library(tseries)
      library(fGarch)
      library(FinTS)
      a=ts(scan("583.txt"))
      ts.plot(a)
      fit=lm(a~-1+time(a))
      r=resid(fit)
      summary(fit)
      pacf(r^2)
      acf(r)
      acf(r^2)
      AutocorTest(r)  #残差是否存在序列相关
      ArchTest(r)     #是否存在ARCH效应
      fit1=garchFit(~arma(2,0)+garch(1,1), data=r, algorithm="nlminb+nm", 
                    trace=F, include.mean=F)
      summary(fit1)
      
      #协整检验   
      fit=arima(b[,2],xreg=b[,1],method="CSS")
      r=resid(fit)
      summary(ur.df(r,type="drift",lag=1))
      Box.test(r,lag=6,fitdf=1)
      
      
      
      【自动运行的自编函数】
      acf.3(x)    #同时绘制3个相关图，acf函数的扩展
      ur.df.01(x)  #进行单位根检验，得到更加舒服的结果
      tsdiag2(x)  #返回x的
      arma.choose(x,ari=3,mai=3)  #选择合适的AR和MA，基于包tseries的arma函数
      
      
      
      #########################附属自编函数
      #...
      acf.3=function(x,lag.max=10,…){
        ol=par(mfrow=c(3,1),mar=c(2,4,1,1))
        acf(x,lag.max=lag.max,type="correlation")
        acf(x,lag.max= lag.max,type="covariance")
        acf(x,lag.max= lag.max,type="partial")
        par(ol)
      }
      #...
      #...类似于tsgiag函数的扩展
      tsdiag2=function(xx.model,fitdf=0,testlag=10){
        t1=xx.arma$residuals
        t2=acf(na.omit(t1),plot=F)
        t3=sapply(1:testlag,
                  function(x,r,fitdf){
                    Box.test(r,type="Ljung-Box",lag=x, fitdf=fitdf)
                  },
                  r=t1,fitdf=fitdf)
        par(mfrow=c(3,1))
        plot(t1,type="b",ylab="",main="残差走势")
        lines(c(0,length(t1)*2),c(0,0),col=2,lty=2)
        plot(t2,type="h",ylab="ACF",main="残差的自相关系数")
        plot(do.call("c",t3[3,]),type="p",ylab="P-value",pch=16,col=4,
             ylim=c(0,1),main="残差的Ljung-Box检验")
        lines(c(0,attr(t1,"tsp")[2]),c(0.05,0.05),lty=2,col=2)
      }
      #...
      ur.df.01=function(x,lags=8){    
        #将三种ADF检验形式汇总的函数（结果和EVIEWS不一致）
        res=matrix(0,5,3)
        colnames(res)=c("无","含常数项","含常数项和趋势项")
        rownames(res)=c("tau统计量","1%临界值","5%临界值",
                        "10%临界值","是否稳定（1/0）")
        types=c("none","drift","trend")
        for(i in 1:3){
          x.adf=ur.df(x,type=types[i],lags=lags,selectlags="AIC")
          x.adf.1=x.adf@teststat  #统计量
          x.adf.2=x.adf@cval      #临界值
          res[1,i]  =x.adf.1[1]
          res[2:4,i]=x.adf.2[1,]
          res[5,i]=if( !is.nan(res[1,i]) & abs(res[1,i]) > abs(res[3,i]) ) 1 else 0
        }
        return(res)
      }
      #...
      #...
      arma.choose.02=function(x){
        #二进制进位运算，以矩阵形式,x=c(0,1,0,1,...)
        n=length(x)
        if( all(!as.logical(x-rep(1,n))) ) stop("已不能再加1!")
        x[1]=x[1]+1
        for(i in 1:(n-1)) if(x[i]>1){ x[i]=0;x[i+1]=x[i+1]+1 }
        return(x)
      }
      arma.choose.01=function(ti){
        #把ti变换成所有可能的ti个0或1的组合
        if(ti<0)  stop("ti要大于0！")
        if(ti==0) return(0)
        if(ti%%1!=0) stop("ti要整数！")
        res=matrix(0,2^ti,ti)
        for(i in 2:2^ti) res[i,]=arma.choose.02(res[i-1,])
        return(res)
      }
      arma.choose.03=function(t0){
        gsub(", ",".",toString(t0,sep=""))
      }
      arma.choose.04=function(i,ari,tti){
        #ari是最大滞后期,tti由ari生成
        ar.lag=((1:ari)*tti[i,])
        ar.lag=ar.lag[ar.lag!=0]
        ar.lag
      }
      arma.choose=function(x,ari=3,mai=3,...){
        tti=arma.choose.01(ari)
        ttj=arma.choose.01(mai)
        ti=2^ari;tj=2^mai
        res.aic=matrix(Inf,ti,tj)     #保存所有组合的AIC
        rownames(res.aic)=paste("AR",apply(tti,1,arma.choose.03),sep=".")
        colnames(res.aic)=paste("MA",apply(ttj,1,arma.choose.03),sep=".")
        res.rss=matrix(Inf,ti,tj)     #保存所有组合的RSS
        rownames(res.rss)=paste("AR",apply(tti,1,arma.choose.03),sep=".")
        colnames(res.rss)=paste("MA",apply(ttj,1,arma.choose.03),sep=".")
        for(i in 2:ti){
          j=1
          ar.lag=arma.choose.04(i,ari,tti)
          x.arma=arma(x,lag=list(ar=ar.lag),...)
          ss=summary(x.arma)
          res.aic[i,j]=ss$aic
          res.rss[i,j]=sum(ss$residuals^2)
        }
        for(j in 2:tj){
          i=1
          ma.lag=arma.choose.04(j,mai,ttj)
          x.arma=arma(x,lag=list(ma=ma.lag),...)
          ss=summary(x.arma)
          res.aic[i,j]=ss$aic
          res.rss[i,j]=sum(ss$residuals^2)
        }
        for(i in 2:ti){for(j in 2:tj){
          ar.lag=arma.choose.04(i,ari,tti)
          ma.lag=arma.choose.04(j,mai,ttj)
          x.arma=arma(x,lag=list(ar=ar.lag,ma=ma.lag),...)
          ss=summary(x.arma)
          res.aic[i,j]=ss$aic
          res.rss[i,j]=sum(ss$residuals^2)
        }}
        res=list()
        res[["tt.ar"]]=tti
        res[["tt.ma"]]=ttj
        temp1=which.min(res.aic)   #找到最小的位置，把res.aic当做按列排的向量
        temp2=temp1 %% ti          #ti是行数，取余以后就是（temp2）行号
        #AR可以直接被arma调用，MA同理
        res[["AR"]]=if(temp2==0) arma.choose.04(ti,ari,tti) else arma.choose.04(temp2,ari,tti)
        res[["MA"]]=arma.choose.04( ceiling( temp1 / ti ), mai,ttj)
        res[["aic"]]=res.aic
        res[["rss"]]=res.rss
        res
      }
      
      #.
      
      出处：http://www.xiaowanxue.com/up_files/2012919161514.html
      
      