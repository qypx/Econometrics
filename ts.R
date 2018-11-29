### 21.3 检验数据集mpyr.dta中各主要变量是否含有单位根
#### 读入数据
library(foreign)
setwd("E://新建文件夹 (2)/计量经济")
mpyr=read.dta("mpyr.dta")
attach(mpyr)
library(tseries)

#### ①ADF检验（原假设为存在单位根）
adf.test(logm1)
adf.test(logp)
adf.test(logy)
adf.test(r)
adf.test(logmr)
adf.test(logv)
#### 由以上结果可见p值均大于0.05，不拒绝原假设，认为各变量均存在单位根

#### ②PP检验（原假设为存在单位根）
pp.test(logm1)
pp.test(logp)
pp.test(logy)
pp.test(r)
pp.test(logmr)
pp.test(logv)
#### 由以上结果可见p值均大于0.05，不拒绝原假设，认为各变量均存在单位根

#### ③KPSS检验（原假设为时间序列为平稳序列）
kpss.test(logm1)
kpss.test(logp)
kpss.test(logy)
kpss.test(r)
kpss.test(logmr)
kpss.test(logv)
#### 由以上结果可见p值均小于0.05，拒绝原假设，认为各变量均存在单位根

detach(mpyr)



### 22.2 数据集dow1.dta包含了1953-1990年美国道琼斯股股指的收盘价。计算道琼斯股指的日收益率，画其时间趋势图，并参照本章实例进行ARCH/GARCH估计。
#### 读入数据
setwd("E://新建文件夹 (2)/计量经济")
dow<-read.dta("dow1.dta")
attach(dow)

#### 计算日收益率：
r<-c()
for(i in 1:nrow(dow)-1){
  r[i]=(dowclose[i+1]-dowclose[i])/dowclose[i]
}


detach(dow)

#### 日收益率的时间趋势图
ts.plot(r)
#### 从图中可以看到，存在波动性聚集

#### 考虑自回归模型
fit1<-ar(r,method = "ols");fit1
#### 选择AR(6)

#### 提取AR模型残差
e<-fit1$resid
e=e[-(1:6)]


#### 检验是否存在ARCH效应
library(FinTS)
ArchTest(e)
#### 拒绝原假设，认为有ARCH效应

#### Ljung-Box Q检验
Box.test(e^2,type="Ljung")
#### 拒绝原假设，认为有ARCH效应

#### 画残差平方的自相关图与偏自相关图
par(mfrow=c(1,2))
acf(e^2)
pacf(e^2)
#### 由图可见，残差平方序列存在自相关，故扰动项存在条件异方差

#### 下面考察GARCH(1,1)模型
library(rugarch)
spec=ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                mean.model = list(armaOrder=c(6,0),include.mean=F,arfima=F))
fit=ugarchfit(spec=spec,r);fit
#### 结果显示，ARCH(1)与GARCH(1)项均很显著（程序中alpha1代表ARCH(1),beta1代表GARCH(1)）

