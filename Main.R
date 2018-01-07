setwd("C:/bujji/SEM2/MATH265/Data")
library(fGarch)
library(fBasics)
library(psych)
library(forecast)
library(PerformanceAnalytics)
library(quantmod)
library(rugarch)
library(car)
library(FinTS)
library(rmgarch)
library(MASS)
options(digits=10)
symbol.vec = c("EBAY", "AMZN", "WMT")
getSymbols(symbol.vec, from ="2011-01-01", to = "2014-12-31")
colnames(EBAY)
start(EBAY)
end(EBAY)
colnames(AMZN)
start(AMZN)
end(AMZN)
colnames(WMT)
start(WMT)
end(WMT)
save(EBAY,file="EBAY.Rdata")
save(AMZN,file="AMZN.Rdata")
save(WMT,file="WMT.Rdata")
# extract adjusted closing prices
EBAY = EBAY[, "EBAY.Adjusted", drop=F]
AMZN = AMZN[, "AMZN.Adjusted", drop=F]
WMT = WMT[, "WMT.Adjusted", drop=F]
# calculate log-returns for GARCH analysis
EBAY.ret = CalculateReturns(EBAY, method="log")
AMZN.ret = CalculateReturns(AMZN, method="log")
WMT.ret = CalculateReturns(WMT, method="log")
save(EBAY.ret,file="EBAYreturn.Rdata")
save(AMZN.ret,file="AMZNreturn.Rdata")
save(WMT.ret,file="WMTreturn.Rdata")
# remove first NA observation
EBAY.ret = EBAY.ret[-1,]
AMZN.ret = AMZN.ret[-1,]
WMT.ret = WMT.ret[-1,]
colnames(EBAY.ret) ="EBAY"
colnames(AMZN.ret) = "AMZN"
colnames(WMT.ret) = "WMT"
# plot the log returns in a single 1x3 plot
par(mfcol=c(1,3))
# plot log returns
plot(EBAY.ret)
plot(AMZN.ret)
plot(WMT.ret)
stock.ret = merge(EBAY.ret, AMZN.ret, WMT.ret)
basicStats(stock.ret)
t.test(stock.ret$EBAY)
t.test(stock.ret$AMZN)
t.test(stock.ret$WMT)
## Testing for skewness ##
skew = abs((sapply(stock.ret,skewness))/(sqrt(6/nrow(stock.ret))))
skew
pvalue = 2*(1-(sapply(abs(skew),pnorm)))
pvalue
kurt = (sapply(stock.ret,kurtosis))/(sqrt(24/nrow(stock.ret)))
kurt
pval = 2*(1-(sapply(abs(kurt),pnorm)))
pval
# plot returns
par(mfcol=c(1,3))
# print density plots
densityPlot(as.timeSeries(EBAY.ret))
densityPlot(as.timeSeries(AMZN.ret))
densityPlot(as.timeSeries(WMT.ret))
# print QQ-norm plots to assess normal distribution
qqnormPlot(as.timeSeries(EBAY.ret))
mtext("EBAY-Gaussian",side=3,cex=0.7,padj=-0.5)
qqnormPlot(as.timeSeries(AMZN.ret))
mtext("AMZN-Gaussian",side=3,cex=0.7,padj=-0.5)
qqnormPlot(as.timeSeries(WMT.ret))
mtext("WMT-Gaussian",side=3,cex=0.7,padj=-0.5)
# print QQ-GHT plots to assess Generalized Hyperbolic Student-t distribution
qqghtPlot(as.timeSeries(EBAY.ret))
mtext("EBAY-Student-t Dist",side=3,cex=0.7,padj=-0.5)
qqghtPlot(as.timeSeries(AMZN.ret))
mtext("AMZN-Student-t Dist",side=3,cex=0.7,padj=-0.5)
qqghtPlot(as.timeSeries(WMT.ret))
mtext("WMT-Student-t Dist",side=3,cex=0.7,padj=-0.5)
# print QQ-GLD plots to assess generalized lambda distribution
qqgldPlot(as.timeSeries(EBAY.ret))
mtext("EBAY-Lambda Dist",side=3,cex=0.7,padj=-0.5)
qqgldPlot(as.timeSeries(AMZN.ret))
mtext("AMZN-Lambda Dist",side=3,cex=0.7,padj=-0.5)
qqgldPlot(as.timeSeries(WMT.ret))
mtext("WMT-Lambda Dist",side=3,cex=0.7,padj=-0.5)

# Use a more robust method for cleaning Outliers from the dataset


EBAY.ret.clean = Return.clean(EBAY.ret, method="boudt")
par(mfrow=c(2,1))
plot(EBAY.ret, main="Raw EBAY Returns", ylab="EBAY")
plot(EBAY.ret.clean, main="Cleaned EBAY Returns", ylab="EBAY")

### There are some outliers.
AMZN.ret.clean = Return.clean(AMZN.ret, method="boudt")
par(mfrow=c(2,1))
plot(AMZN.ret, main="Raw AMZN Returns", ylab="AMZN")
plot(AMZN.ret.clean, main="Cleaned AMZN Returns", ylab="AMZN")

### There are no outliers.
WMT.ret.clean = Return.clean(WMT.ret, method="boudt")
par(mfrow=c(2,1))
plot(WMT.ret, main="Raw WMT Returns", ylab="WMT")
plot(WMT.ret.clean, main="Cleaned WMT Returns", ylab="WMT")

# The BoxCox.lambda() function will choose a value of lambda for you.
lambda <- BoxCox.lambda(AMZN.ret.clean)
#lambda = 0.794
lambda <- BoxCox.lambda(WMT.ret)
#lambda = 
lambda <- BoxCox.lambda(EBAY.ret)
#lambda = 
###There is no good transformation.

# Model identification
par(mfcol=c(1,3))
acf(EBAY.ret)
acf(AMZN.ret.clean)
acf(WMT.ret)

pacf(EBAY.ret)
pacf(AMZN.ret.clean)
pacf(WMT.ret)

e1auto=ar(as.ts(EBAY.ret),method="mle")
e1auto$order
e1auto

a1auto=ar(as.ts(AMZN.ret),method="mle")
a1auto$order
a1auto

w1auto=ar(as.ts(WMT.ret),method="mle")
w1auto$order
w1auto

# series. Include details of each step of the process, and support your final model selection for
# each series.
e1auto <-auto.arima(EBAY.ret)
summary(e1auto)
par(mfcol=c(1,1))
tsdiag(e1auto)
mtext("Diagnostics for EBAY log Returns",side=3,cex=1.0,padj=-2)
sqrt(e1auto$sigma2) # Calculate the residual standard error =
Box.test(e1auto$resid,lag=24,type='Ljung')

a2=arima(AMZN.ret.clean, order=c(0,0,1),include.mean=F)
summary(a2)
a1=arima(AMZN.ret.clean, order=c(3,0,3),include.mean=F)
summary(a1)

a1auto <-auto.arima(AMZN.ret.clean)
summary(a1auto)
coef(a1auto)
tsdiag(a1auto)
mtext("Diagnostics for AMZN log Returns",side=3,cex=1.0,padj=-2)
sqrt(a1auto$sigma2) # Calculate the residual standard error =
Box.test(a1auto$resid,lag=24,type='Ljung')

w1auto <-auto.arima(WMT.ret)
summary(w1auto)
coef(w1auto)
tsdiag(w1auto)
mtext("Diagnostics for WMT log Returns",side=3,cex=1.0,padj=-2)
sqrt(w1auto$sigma2) # Calculate the residual standard error =
Box.test(w1auto$resid,lag=24,type='Ljung')

#g
e1.forecast <- forecast(e1auto, h=21,)
a1.forecast <- forecast(a1auto, h=21,)
w1.forecast <- forecast(w1auto, h=21,)
e1.forecast
a1.forecast
w1.forecast
par(mfcol=c(3,1))
plot(e1.forecast, include=100,main="Forecast of EBAY Log Return Data for January 2015")
plot(a1.forecast, include=100, main="Forecast of WMT Log Return Data for January 2015")
plot(w1.forecast,include=100, main="Forecast of Microsoft Log Return Data for January 2015")

#h)
Box.test(e1auto$resid^2,lag=24,type='Ljung')
Box.test(a1auto$resid^2,lag=24,type='Ljung')
Box.test(w1auto$resid^2,lag=24,type='Ljung')

Box.test(EBAY.ret^2,lag=12,type='Ljung')
Box.test(AMZN.ret.clean^2,lag=24,type='Ljung')
Box.test(WMT.ret^2,lag=24,type='Ljung')

ArchTest(EBAY.ret^2, lags=24)
ArchTest(AMZN.ret.clean^2, lags=24)
ArchTest(WMT.ret^2, lags=24)

par(mfcol=c(1,3))
acf(EBAY.ret^2)
acf(AMZN.ret^2)
acf(WMT.ret^2)

pacf(EBAY.ret^2)
pacf(AMZN.ret^2)
pacf(WMT.ret^2)

par(mfcol=c(1,2))
acf(AMZN.ret.clean^2)
pacf(AMZN.ret.clean^2)


squared.res.e1=e1auto$resid^2
squared.res.a1=a1auto$resid^2
squared.res.w1=w1auto$resid^2
par(mfcol=c(3,1))
plot(squared.res.e1,main='Squared Residuals')
plot(squared.res.a1,main='Squared Residuals')
plot(squared.res.w1,main='Squared Residuals')


p)
dcc.fcst = dccforecast(dcc.fit, n.ahead=21)
class(dcc.fcst)
slotNames(dcc.fcst)
class(dcc.fcst@mforecast)
names(dcc.fcst@mforecast)
plot(dcc.fcst)
#########################################################################################
# #
# Plot the DCC Forecasted Conditional Covariances. #
# Create 1 plot for each set of pairs. #
# #
#########################################################################################
plot(dcc.fcst, which=3, series=c(1,2))
plot(dcc.fcst, which=3, series=c(1,3))
plot(dcc.fcst, which=3, series=c(2,3))

plot(dcc.fcst, which=3, series=c(1,2))
plot(dcc.fcst, which=3, series=c(1,3))
plot(dcc.fcst, which=3, series=c(2,3))

# firstmonth in 2015's 1-step ahead forecasts of conditional covariance and correlation. #
# #
# #
####################################################################################################
#
# forecasting conditional volatility and correlations
#
#########################################################################################
# #
# Plotting workaround for DCC Forecasted Conditional Correlation Plot #
# Create 1 plot for each set of pairs. #
# #
#########################################################################################
# Extract the forecasted conditional correlations into a dataframe.
corr.dcc.fcst = as.data.frame(rcor(dcc.fcst))
names(corr.dcc.fcst)
row.names(corr.dcc.fcst)
corr.dcc.fcst.df <- as.data.frame(t(corr.dcc.fcst)) # Transpose the rows and columns in the dataframe
names(corr.dcc.fcst.df)
# Create binary flags for subsetting the dataframe to extract desired pairs.
# The dataframe repeats the rows of the correlation matrix: EBAY, AMZN< WMT for each day.
# To extract the desired pair, create a flag=1 to identify the desired row.
corr.dcc.fcst.df$X2014.12.31.AMZN <- as.vector(rep(c(0,1,0),times=21)) # Forcast horizon was 21 days for Jan 2015.
corr.dcc.fcst.df$X2014.12.31.WMT <- as.vector(rep(c(0,0,1),times=21))
# subset the dataframe extract the desired correlation pairs from the forcast dataframe
plot.EBAY.AMZN.df <- subset(corr.dcc.fcst.df,corr.dcc.fcst.df$X2014.12.31.AMZN==1, select=EBAY )
plot.EBAY.WMT.df <- subset(corr.dcc.fcst.df,corr.dcc.fcst.df$X2014.12.31.WMT==1, select=EBAY )
plot.AMZN.WMT.df <- subset(corr.dcc.fcst.df,corr.dcc.fcst.df$X2014.12.31.WMT==1, select=AMZN )

# Extract the fitted conditional correlations into a dataframe.
corr.dcc.fit = as.data.frame(rcor(dcc.fit))
corr.dcc.fit.df <- as.data.frame(t(corr.dcc.fit)) # transposte the rows and columns
corr.dcc.fit.df[1:10,] # verify the transposed columns
corr.dcc.fit.df$X2014.12.31.AMZN <- as.vector(rep(c(0,1,0),times=1005)) # Fitted data contained 1005 observations.
corr.dcc.fit.df$X2014.12.31.WMT <- as.vector(rep(c(0,0,1),times=1005))
# subset the dataframe to extract the desired conditional correlation pairs from the fitted dataframe.
plot.EBAY.AMZN.fit.df <- subset(corr.dcc.fit.df,corr.dcc.fit.df$X2014.12.31.AMZN==1, select=EBAY )
plot.EBAY.WMT.fit.df <- subset(corr.dcc.fit.df,corr.dcc.fit.df$X2014.12.31.WMT==1, select=EBAY )
plot.AMZN.WMT.fit.df <- subset(corr.dcc.fit.df,corr.dcc.fit.df$X2014.12.31.WMT==1, select=AMZN )

# Combine the fitted and forecasted conditional correlations into a single dataframe.
plot.EBAY.AMZN.fit.df <- rbind(plot.EBAY.AMZN.fit.df,plot.EBAY.AMZN.df)
plot.EBAY.WMT.fit.df <- rbind(plot.EBAY.WMT.fit.df,plot.EBAY.WMT.df)
plot.AMZN.WMT.fit.df <- rbind(plot.AMZN.WMT.fit.df,plot.AMZN.WMT.df)
library(ggplot2)
# Generate Correlation Plots
# First set-up a dataframe for printing the plot with line colors then generage the plot
#
# Plot Ebay-AMZNle Correlation
#
row.names(plot.EBAY.AMZN.fit.df)[940:950]
x.start <- 942 # begin plotting in the 3rd quarter of 2014 (October 1,2014)
x.end <- nrow(plot.EBAY.AMZN.fit.df)
fit.length <- 1005 - x.start +1
forcast.length <- nrow(plot.EBAY.AMZN.df)
line.color <- c(rep("darkgrey",fit.length), # 1005 - days of fitted observations
                rep("red", forcast.length)) # 21 - days in the forecast
df <- data.frame(
  x = x.start:x.end, # 1:nrow(plot.EBAY.AMZN.fit.df),
  y = plot.EBAY.AMZN.fit.df$EBAY[x.start:x.end]
  col = line.color)
p1 <- ggplot(df, aes(x=x, y=y)) +
  geom_line(aes(colour=col, group=1)) +
  scale_colour_identity() +
  labs(y="Correlation", x = "Days") +
  ggtitle("Ebay-Amazon Forecasted Conditional Correlation\n October 2014 through January 2015")
p1 + theme_bw() + theme(plot.title = element_text(color="black",size=20,face="bold"))
#
# Plot Ebay-Microsoft Correlation
#
df <- merge(
  x = x.start:x.end,
  y = plot.EBAY.WMT.fit.df$EBAY[x.start:x.end],
  col = line.color)
p2 <- ggplot(df, aes(x=x, y=y)) +
  geom_line(aes(colour=col, group=1)) +
  scale_colour_identity() +
  labs(y="Correlation", x = "Days") +
  ggtitle("Ebay-Microsoft Forecasted Conditional Correlation\n October 2014 through January 2015")
p2 + theme_bw() + theme(plot.title = element_text(color="black",size=20,face="bold"))
#
# Plot AMZNle-Microsoft Correlation
#
df <- data.frame(
  x = x.start:x.end,
  y = plot.AMZN.WMT.fit.df$AMZN[x.start:x.end],
  col =line.color)
p3 <- ggplot(df, aes(x=x, y=y)) +
  geom_line(aes(colour=col, group=1)) +
  scale_colour_identity() +
  labs(y="Correlation", x = "Days") +
  ggtitle("AMZNle-Microsoft Forecasted Conditional Correlation\n October 2014 through January 2015")
p3 + theme_bw() + theme(plot.title = element_text(color="black",size=20,face="bold"))


# (q) Compare your mean and volatility forecast from part (p) with parts (g) and (k). Which is the
# best model?
#
###################################################################
# #
# Using the student distribution for the garch specification #
# #
###################################################################
garchstd.spec=ugarchspec(mean.model = list(armaOrder = c(0,0)),
                         variance.model=list(garchOrder = c(1,1), model="sGARCH"), distribution.model="std")
# dcc specification - GARCH(1,1) for conditional correlations
dcc.garchstd.spec = dccspec(uspec=multispec(replicate(3,garchstd.spec)), dccOrder = c(1,1), distribution="mvt")
dcc.garchstd.spec
dccstd.fit = dccfit(dcc.garchstd.spec, stock.ret)
show(dccstd.fit)
coef(dccstd.fit)
likelihood(dccstd.fit)




accuracy(e1.forecast)




