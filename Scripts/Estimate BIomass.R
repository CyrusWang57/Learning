
#### predict biomass ####
#########################
setwd("E:/wangc/Documents/R/Data/Estimate Biomass") #must use "/" instead of "\"
install.packages("readxl")
library(readxl)
coexist <- read_excel('coexist.xlsx') # or read.table('coexist.txt',header = T)
plot(bio ~ height, data=coexist)
plot(bio ~ diameter, data=coexist)
lm_co <- lm(bio ~ height + diameter, data=coexist)#lm is used to fit linear models.
plot(lm_co)
summary(lm_co)
lrt_co <- drop1(lm_co, test="Chi")
lrt_co

newdata=data.frame(height=c(13, 15, 17, 18), diameter=c(3, 4, 5, 6))
newdata
predict(lm_co, newdata=newdata, type = "response")

