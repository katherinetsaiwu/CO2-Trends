# Group Project I: Trends in Atmospheric Carbon Dioxide (#5)
# STP 429: Dassanayake - T/TH 1:30-2:45
# Team: Danny Ober-Reynolds, Jimin Nam, Katherine Wu, TJ Radigan

# Part I: The Data

# Set Directory (different for everyone)
# setwd("C:/Users/dober_000/Documents/School 2016 Spring/STP 429/Project1")
setwd("C:/Users/dober_000/Documents/School 2016 Spring/STP 429/Project1/")

# Import Data
# Source: http://www.esrl.noaa.gov/gmd/ccgg/trends/
# Sets: Mauna Loa CO2 monthly mean data and Mauna Loa CO2 annual mean data
monthlydta = read.table("co2_mm_mlo.txt")
annualdta = read.table("co2_annmean_mlo.txt")

# Rename Columns
names(monthlydta) <- c("year", "month", "decimaldate", "average", "interpolated", "trend", "ndays")
names(annualdta) <- c("year", "mean", "unc")

# Replace missing observations with "NA"
# -99.99 and -1 denote missing
monthlydta[monthlydta == -99.99] <- NA
monthlydta[monthlydta == -1] <- NA


# Part II: The Model

# 1) simple linear regression model on the annual mean dataset with annual mean CO2
simplereg = lm(mean ~ year, data=annualdta)
# residuals for model 1
plot(fitted(simplereg), residuals(simplereg), main = "SLR: mean ~ year, Fitted vs. Residuals", xlab = "Fitted", ylab = "Residuals")

# Diagnostic graphs (saves to working directory)
png('SimpleRegDiagnosticGraphs.png',
width = 640, height = 640, units = "px")
par(mfrow=c(2,2))
plot(simplereg)
dev.off()

# Normal QQ plot (saves to working direectory)
png('Graph1a_SimpleRegQQPlot.png',
width = 640, height = 640, units = "px")
qqnorm(residuals(simplereg), main = "Graph 1(a): Normal QQ Plot")
qqline(residuals(simplereg))
dev.off()
# Test normality assumption
shapiro.test(residuals(simplereg))

# We see a pattern in residual vs fitted plot. We'll try transforming the independent variable

# standardize year
annualdta$stdYear = (annualdta$year - mean(annualdta$year))/sd(annualdta$year)
# square
annualdta$yearSqrd = annualdta$year * annualdta$year
# other year variables (staring from 1958 = 0)
annualdta$yearFrom1958 = annualdta$year - 1958
annualdta$yearFrom1958Sqrd = annualdta$yearFrom1958 * annualdta$yearFrom1958
annualdta$yearFrom1958Sqrt = sqrt(annualdta$yearFrom1958)
annualdta$yearFrom1958log = log(annualdta$yearFrom1958)
annualdta$stYearFrom1958 = (annualdta$yearFrom1958 - mean(annualdta$yearFrom1958))/sd(annualdta$yearFrom1958)
annualdta$invYearFrom1958 = 1/annualdta$yearFrom1958
#Tranformed mean (dependent variable)
annualdta$sqrtMean = sqrt(annualdta$mean)
annualdta$sqrdMean = annualdta$mean * annualdta$mean
annualdta$logMean = log(annualdta$mean)
annualdta$meanPowNeg2.75 = annualdta$mean^(-2.75)


# Alternative SLR models
# Standardized year
simpleregDeMeanedYr <- lm(mean ~ stdYear, data = annualdta)
par(mfrow=c(2,2))
plot(simpleregDeMeanedYr)
shapiro.test(residuals(simpleregDeMeanedYr))

# Year squared
simpleregSqrdYr <- lm(mean ~ yearSqrd, data = annualdta)
par(mfrow=c(2,2))
plot(simpleregSqrdYr)
shapiro.test(residuals(simpleregSqrdYr))

# Year from 1958 squared (THIS MODEL HAS THE BEST NORMALITY ASSUMPTION)
simpleregYrFm1958Sqrd <- lm(mean ~ yearFrom1958Sqrd, data = annualdta)
par(mfrow=c(2,2))
plot(simpleregYrFm1958Sqrd)
shapiro.test(residuals(simpleregYrFm1958Sqrd))

# Year from 1959 sqrt
simpleregYrFm1958Sqrt <- lm(mean ~ yearFrom1958Sqrt, data = annualdta)
par(mfrow=c(2,2))
plot(simpleregYrFm1958Sqrt)
shapiro.test(residuals(simpleregYrFm1958Sqrt))

# Log year from 1958
simpleregLogYrFrom1958 <- lm(mean ~ yearFrom1958log, data = annualdta)
par(mfrow=c(2,2))
plot(simpleregLogYrFrom1958)
shapiro.test(residuals(simpleregLogYrFrom1958))

# Standardized year from 1958
simpleregStYrFrom1958 <- lm(mean ~ stYearFrom1958, data = annualdta)
par(mfrow=c(2,2))
plot(simpleregStYrFrom1958)
shapiro.test(residuals(simpleregStYrFrom1958))

# Inverse year from 1958
simpleregInvYearFrom1958 <- lm(mean ~ invYearFrom1958, data = annualdta)
par(mfrow=c(2,2))
plot(simpleregInvYearFrom1958)
shapiro.test(residuals(simpleregInvYearFrom1958))

# Sqrt mean (dependent variable)
simpleregSqrtMean <- lm(sqrtMean ~ year, data = annualdta)
par(mfrow=c(2,2))
plot(simpleregSqrtMean)
shapiro.test(residuals(simpleregSqrtMean))

# Sqrd mean (dependent variable)
simpleregSqrdMean <- lm(sqrdMean ~ year, data = annualdta)
par(mfrow=c(2,2))
plot(simpleregSqrdMean)
shapiro.test(residuals(simpleregSqrdMean))

# Log mean (dependent variable)
simpleregLogMean <- lm(logMean ~ year, data = annualdta)
par(mfrow=c(2,2))
plot(simpleregLogMean)
shapiro.test(residuals(simpleregLogMean))

# Double log model
simpleregDoubleLog <- lm(logMean ~ yearFrom1958log, data = annualdta)
par(mfrow=c(2,2))
plot(simpleregDoubleLog)
shapiro.test(residuals(simpleregDoubleLog))

# sqrt mean, log year from 1958
simpleregSqrtMeanLogYr <- lm(sqrtMean ~ yearFrom1958log, data = annualdta)
par(mfrow=c(2,2))
plot(simpleregSqrtMeanLogYr)
shapiro.test(residuals(simpleregSqrtMeanLogYr))

# mean^-2.75, year
simpleregMeanPowNeg275 <- lm(meanPowNeg2.75 ~ year, data = annualdta)
par(mfrow=c(2,2))
plot(simpleregMeanPowNeg275)
shapiro.test(residuals(simpleregMeanPowNeg275)) # results: p-value 0.2475

# Save graphs from mean^-2.75 model
png('SimpleRegPowModelDiagnosticGraphs.png',
width = 640, height = 640, units = "px")
par(mfrow=c(2,2))
plot(simpleregMeanPowNeg275)
dev.off()

# Normal QQ plot (saves to working direectory)
png('Graph1b_SimpleRegPowModelQQPlot.png',
width = 640, height = 640, units = "px")
qqnorm(residuals(simpleregMeanPowNeg275), main = "Graph 1(b): Normal QQ Plot")
qqline(residuals(simpleregMeanPowNeg275))
dev.off()



# The last model, simpleregMeanPowNeg275, appears to have the least issues with the normality assumption. However, there still appears to be a sort of pattern in the residuals vs. fitted, indicating heteroskedasticity. Perhaps including the months as a factor variable will help.


# 2) multiple regression model on the monthly mean dataset with monthly mean CO2
multiplereg = lm(average ~ year + factor(month), data=monthlydta)
# check residuals
par(mfrow=c(2,2))
plot(multiplereg)
# get residuals - use qqnorm

# We still see the same issue as before, using mean rather than mean^-2.75.

monthlydta$averagePowNeg2.75 = monthlydta$average^(-2.75)

# 2a) multiple regression model on the monthly mean dataset with monthly mean CO2 and mean^-2.75 as response variable
multipleregMeanPowNeg275 <- lm(averagePowNeg2.75 ~ year + factor(month), data = monthlydta)
par(mfrom=c(2,2))
plot(multipleregMeanPowNeg275)
shapiro.test(residuals(multipleregMeanPowNeg275)) # results: p-value 0.01407
summary(multipleregMeanPowNeg275)


# Part III: Graphs
# 1) For the simple linear regression model: scatter plot of the annual mean CO2 data, with the predicted line and its 95% confidence band
attach(annualdta)
predictlineSimpleReg <- predict(simplereg, interval="confidence")
png('SimpleRegScatterAndPrediction.png',
width = 640, height = 640, units = "px")
plot (mean ~ year, data=annualdta, main = "Graph 2(a): Simple Pre-Transformation Model")
lines(year[order(year)], predictlineSimpleReg[order(year), 1])
lines(year[order(year)], predictlineSimpleReg[order(year), 2], col="blue", lty=2)
lines(year[order(year)], predictlineSimpleReg[order(year), 3], col="blue", lty=2)
dev.off() # save graph (export)

# 1a) For the fixed simple linear regression model (simpleregMeanPowNeg275): scatter plot of the annual mean CO2 data, with the predicted line and its 95% confidence bands
predictlineSimpleRegPowModel <- predict(simpleregMeanPowNeg275, interval = "confidence")
png('SimpleRegPowModelScatterAndPrediction.png',
    width = 640, height = 640, units = "px")
plot (meanPowNeg2.75 ~ year, main = "Graph 2(b): Simple Model (y^-2.75)")
lines(year[order(year)], predictlineSimpleRegPowModel[order(year), 1])
lines(year[order(year)], predictlineSimpleRegPowModel[order(year), 2], col="blue", lty=2)
lines(year[order(year)], predictlineSimpleRegPowModel[order(year), 3], col="blue", lty=2)
dev.off()
detach(annualdta)

# 2) For the multiple regression model: connected line plot of the monthly mean CO2 data to show the trend
fixYearMonthWrong <- function(yearMonthWrong) {
  if (nchar(yearMonthWrong) < 6){
    year <- substr(yearMonthWrong,0,4)
    month <- substr(yearMonthWrong,5,5)
    result <- as.integer(paste(year,"0",month, sep = ""))
  }
  else {
    result <- as.integer(yearMonthWrong)
  }
  return (result)
}

monthlydta$yearMonthWrong = paste(monthlydta$year, monthlydta$month, sep = "")
monthlydta$yearMonth = unlist(lapply(monthlydta$yearMonthWrong, fixYearMonthWrong))
monthlydta <- monthlydta[order(monthlydta$yearMonth),]

# Plot Monthly Mean CO2
png('MontlyMeanCO2.png',
    width = 640, height = 640, units = "px")
plot(0,0,type = "n", xlim = c(0,696), ylim = c(min(monthlydta$average, na.rm = TRUE), max(monthlydta$average, na.rm = TRUE)), axes = FALSE, ann = FALSE)
lines(monthlydta$average, type = "l", col = "blue")
title(main = "Monthly Mean CO2", xlab = "Month", ylab = "Mean CO2")
axis(1, at = c(seq(from = 1, to = 700, by = 70)), labels = c("1958 Mar","1964 Jan","1969 Nov","1975 Sep","1981 Jul","1987 May","1993 Mar","1999 Jan","2004 Nov","2010 Sep"))
axis(2, at = c(seq(310,410,by=10)))
dev.off()


## Bin and histogram CO2 values annually
binsAnnually <- 10
cutpointsAnnually <- quantile(annualdta$mean, (0:binsAnnually)/binsAnnually, na.rm = TRUE)
binnedAnnually <- cut(annualdta$mean, cutpointsAnnually, include.lowest = TRUE, na.rm = TRUE)
summary(binnedAnnually)
png('AnnualCO2Distribution.png', width = 640, height = 640, units = "px")
plot(binnedAnnually, main = "Annual CO2 Distribution")
dev.off()

## Bin and histogram CO2 values monthly
binsMonthly <- 10
cutpointsMonthly <- quantile(monthlydta$average, (0:binsMonthly)/binsMonthly, na.rm = TRUE)
binnedMonthly <- cut(monthlydta$average, cutpointsMonthly, include.lowest =TRUE, na.rm = TRUE)
summary(binnedMonthly)
png('MonthlyCO2Distribution.png', width = 640, height = 640, units = "px")
plot(binnedMonthly, main = "Monthly CO2 Distribution")
dev.off()

# Calculate predicution values from the monthly model. Coefficients:
#  (Intercept)            year  factor(month)2  factor(month)3  factor(month)4  factor(month)5  factor(month)6 
#2.451319e-06   -1.181988e-09   -5.999201e-10   -1.342449e-09   -2.340740e-09   -2.795695e-09   -2.318805e-09 
#factor(month)7  factor(month)8  factor(month)9 factor(month)10 factor(month)11 factor(month)12 
#-1.193374e-09    4.507518e-10    1.763915e-09    1.805585e-09    7.221231e-10   -3.499770e-10 
jan = (2.451319e-06 + -1.181988e-09*2016)^(-1/2.75)
feb = (2.451319e-06 + -1.181988e-09*2016 - 5.999201e-10)^(-1/2.75)
mar = (2.451319e-06 + -1.181988e-09*2016 - 1.342449e-09)^(-1/2.75)
apr = (2.451319e-06 + -1.181988e-09*2016 - 2.340740e-09)^(-1/2.75)
may = (2.451319e-06 + -1.181988e-09*2016 - 2.795695e-09)^(-1/2.75)
jun = (2.451319e-06 + -1.181988e-09*2016 - 2.318805e-09)^(-1/2.75)
jul = (2.451319e-06 + -1.181988e-09*2016 - 1.193374e-09)^(-1/2.75)
aug = (2.451319e-06 + -1.181988e-09*2016 + 4.507518e-10)^(-1/2.75)
sep = (2.451319e-06 + -1.181988e-09*2016 + 1.763915e-09)^(-1/2.75)
oct = (2.451319e-06 + -1.181988e-09*2016 + 1.805585e-09)^(-1/2.75)
nov = (2.451319e-06 + -1.181988e-09*2016 + 7.221231e-10)^(-1/2.75)
dec = (2.451319e-06 + -1.181988e-09*2016 - 3.499770e-10)^(-1/2.75)

