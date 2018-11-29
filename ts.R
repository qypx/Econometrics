### 21.3 �������ݼ�mpyr.dta�и���Ҫ�����Ƿ��е�λ��
#### ��������
library(foreign)
setwd("E://�½��ļ��� (2)/��������")
mpyr=read.dta("mpyr.dta")
attach(mpyr)
library(tseries)

#### ��ADF���飨ԭ����Ϊ���ڵ�λ����
adf.test(logm1)
adf.test(logp)
adf.test(logy)
adf.test(r)
adf.test(logmr)
adf.test(logv)
#### �����Ͻ���ɼ�pֵ������0.05�����ܾ�ԭ���裬��Ϊ�����������ڵ�λ��

#### ��PP���飨ԭ����Ϊ���ڵ�λ����
pp.test(logm1)
pp.test(logp)
pp.test(logy)
pp.test(r)
pp.test(logmr)
pp.test(logv)
#### �����Ͻ���ɼ�pֵ������0.05�����ܾ�ԭ���裬��Ϊ�����������ڵ�λ��

#### ��KPSS���飨ԭ����Ϊʱ������Ϊƽ�����У�
kpss.test(logm1)
kpss.test(logp)
kpss.test(logy)
kpss.test(r)
kpss.test(logmr)
kpss.test(logv)
#### �����Ͻ���ɼ�pֵ��С��0.05���ܾ�ԭ���裬��Ϊ�����������ڵ�λ��

detach(mpyr)



### 22.2 ���ݼ�dow1.dta������1953-1990����������˹�ɹ�ָ�����̼ۡ��������˹��ָ���������ʣ�����ʱ������ͼ�������ձ���ʵ������ARCH/GARCH���ơ�
#### ��������
setwd("E://�½��ļ��� (2)/��������")
dow<-read.dta("dow1.dta")
attach(dow)

#### �����������ʣ�
r<-c()
for(i in 1:nrow(dow)-1){
  r[i]=(dowclose[i+1]-dowclose[i])/dowclose[i]
}


detach(dow)

#### �������ʵ�ʱ������ͼ
ts.plot(r)
#### ��ͼ�п��Կ��������ڲ����Ծۼ�

#### �����Իع�ģ��
fit1<-ar(r,method = "ols");fit1
#### ѡ��AR(6)

#### ��ȡARģ�Ͳв�
e<-fit1$resid
e=e[-(1:6)]


#### �����Ƿ����ARCHЧӦ
library(FinTS)
ArchTest(e)
#### �ܾ�ԭ���裬��Ϊ��ARCHЧӦ

#### Ljung-Box Q����
Box.test(e^2,type="Ljung")
#### �ܾ�ԭ���裬��Ϊ��ARCHЧӦ

#### ���в�ƽ���������ͼ��ƫ�����ͼ
par(mfrow=c(1,2))
acf(e^2)
pacf(e^2)
#### ��ͼ�ɼ����в�ƽ�����д�������أ����Ŷ�����������췽��

#### ���濼��GARCH(1,1)ģ��
library(rugarch)
spec=ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                mean.model = list(armaOrder=c(6,0),include.mean=F,arfima=F))
fit=ugarchfit(spec=spec,r);fit
#### �����ʾ��ARCH(1)��GARCH(1)�����������������alpha1����ARCH(1),beta1����GARCH(1)��
