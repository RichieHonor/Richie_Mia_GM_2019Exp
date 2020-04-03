
#Accounting for amoung date sampling variation. 
```{r}

#Removing Runs that i know are bad. 

#Removing all data from runs on Dec 6 2019 and Dec 8 2019 
GF<-GF[GF$DataFile!="Dat_Dec_6_2019.txt"& GF$DataFile!="Dat_Dec_8_2019.txt",]

#Removing only glucosinolate and flavonoid runs from Dat_Jan_30_2020_PM.txt, Dat_Jan_31_2020_AM.txt, Dat_Jan_31_2020_PM.txt
GF[GF$DataFile=="Dat_Jan_30_2020_PM.txt" | GF$DataFile=="Dat_Jan_31_2020_AM.txt"| GF$DataFile=="Dat_Jan_31_2020_PM.txt","Conc"]<-NA

#Removing Flavonoid data from Mabels Run on Feb 3rd 2020 AM.
GF$Conc[GF$DataFile=="Dat_Feb_3_2020_AM.txt" & GF$Compound=="flav"]<-NA



#Glucosinolate
gluc<-GF[GF$Compound=="gluc",]

#Flavonoid
flav<-GF[GF$Compound=="flav",]

#Making glucosinolate specific column names
colnames(gluc)<-paste("gluc",colnames(gluc),sep="_")

#Making flavonoid specific column names
colnames(flav)<-paste("flav",colnames(flav),sep="_")

#Removing controls from chlorophyll
ChlorNoCont<-ChlorAB[ChlorAB$ID!="|cont|",]

#ISolating genotypes we care about for the fertilizer experiment and all of the pools from the dataset. 

#flavonoid Dataframe.
flav<-flav[grepl("SMITH1-3-80|MSMID1-2-0|KVEDG1-8-60|pool",flav$flav_Tag),]

#Glucosinolate Dataframe.
gluc<-gluc[grepl("SMITH1-3-80|MSMID1-2-0|KVEDG1-8-60|pool",gluc$gluc_Tag),]

#Chlorophyll Dataframe.
ChlorNoCont<-ChlorNoCont[grepl("SMITH1-3-80|MSMID1-2-0|KVEDG1-8-60|pool",ChlorNoCont$Tag),]

#Assign family column
prefamily<-gsub("*.\\|","",ChlorNoCont$Tag)
ChlorNoCont$Family<-gsub("\\-.*","",prefamily)
#Visualizing variation amoung leaves. There can be alot. 

#Accounting for Variation........
#---------------------------------------------------

#Assigning pool collumns for visualization
#GlucPools
gluc$Pool[gluc$gluc_Tag=="pool"]<-"y"
gluc$Pool[is.na(gluc$Pool)]<-"n"
#FlavPools
flav$Pool[flav$flav_Tag=="pool"]<-"y"
flav$Pool[is.na(flav$Pool)]<-"n"
#ChlorPools
ChlorNoCont$Pool[ChlorNoCont$Tag=="pool"]<-"y"
ChlorNoCont$Pool[is.na(ChlorNoCont$Pool)]<-"n"


#Glucosinolates!!!
#---------------------------------------------------

#As you can see from these graphs, there is variation in the pools, but there is more variation in the samples, which is good. 
ggplot(gluc)+
  geom_point(aes(x=gluc$Pool,y=gluc_Conc))

ggplot(flav)+
  geom_point(aes(x=flav$Pool,y=flav_Conc))

#I will now account for the variation in sampling date, using the data from the pools
glucPool<-gluc[gluc$Pool=="y",]
mean(glucPool$gluc_Conc,na.rm=T)
median(glucPool$gluc_Conc,na.rm=T)
#The mean (Red) and median (blue) are approximately equal. The difference between each sample date mean (of pools) and the grand mean (of pools) will be used to determine the correction that will be applied to all samples from that date. 
glucPool

ggplot(glucPool)+
  geom_histogram(aes(x=gluc_Conc),bins=20)+
  geom_vline(aes(xintercept=1.107457,colour="red"))+
  geom_vline(aes(xintercept=1.105041,colour="blue"))

#Accounting for the sampling day variation.  
glucPoolMean<-glucPool %>% group_by(gluc_DataFile) %>% summarize(gluc_Conc=mean(gluc_Conc))

#Determining Adjustment values by the difference from the mean of the sample date to the median of the sample.
glucPoolMean$DateAdjust<-glucPoolMean$gluc_Conc-1.105041

#rejoining the adjustment data with the gluc dataset. 

gluc<-gluc %>% left_join(select(glucPoolMean,gluc_DataFile,DateAdjust),by="gluc_DataFile")


###Flavonols!
#---------------------------------------------------

#I will now account for the variation in sampling date, using the data from the pools
flavPool<-flav[flav$Pool=="y",]
mean(flavPool$flav_Conc,na.rm=T)
median(flavPool$flav_Conc,na.rm=T)
#The mean (Red) and median (blue) are approximately equal. The difference between each sample date mean (of pools) and the grand mean (of pools) will be used to determine the correction that will be applied to all samples from that date. 
flavPool

ggplot(flavPool)+
  geom_histogram(aes(x=flav_Conc),bins=25)+
  geom_vline(aes(xintercept=1.570148,colour="red"))+
  geom_vline(aes(xintercept=1.570296,colour="blue"))

#Accounting for the sampling day variation.  
flavPoolMean<-flavPool %>% group_by(flav_DataFile) %>% summarize(flav_Conc=mean(flav_Conc))

#Determining Adjustment values by the difference from the mean of the sample date to the median of the sample.
flavPoolMean$DateAdjust<-flavPoolMean$flav_Conc-1.105041

#rejoining the adjustment data with the flav dataset. 

flav<-flav %>% left_join(select(flavPoolMean,flav_DataFile,DateAdjust),by="flav_DataFile")

gluc$gluc_Conc_Adj<-gluc$gluc_Conc+gluc$DateAdjust
flav$flav_Conc_Adj<-flav$flav_Conc+flav$DateAdjust



#Visualizing how accounting for date with pool data affect the dispersal of variation in the dataset.  
```{r}
#Removing pools 
gluc<-gluc[gluc$Pool=="n",]
flav<-flav[flav$Pool=="n",]

#visualizing again. 
ggplot(gluc)+
  geom_boxplot(aes(x=gluc_Sample,y=gluc_Conc_Adj))
ggplot(flav)+
  geom_boxplot(aes(x=flav_Sample,y=flav_Conc_Adj))

glucFit<-lmer(gluc_Conc_Adj~1+(1|gluc_Tag/gluc_Sample/gluc_DataFile),data=gluc)
flavFit<-lmer(flav_Conc_Adj~1+(1|flav_Tag/flav_Sample/flav_DataFile),data=flav)
summary(glucFit)
summary(flavFit)

glucVar<-data.frame(Variance=c(0.041324, 0.014003,0.024288,0.008855),Variable=c("AmoungDate","AmoungLeaf", "AmoungPlant", "TechnicalReplicate"))
glucVar$Percent<-glucVar$Variance/sum(glucVar$Variance)*100

flavVar<-data.frame(Variance=c(0.042491,0.017512,0.012571,0.007988),Variable=c("AmoungDate","AmoungLeaf", "AmoungPlant", "TechnicalReplicate"))
flavVar$Percent<-flavVar$Variance/sum(flavVar$Variance)*100

ggplot(glucVar)+
  geom_col(aes(x=Variable,y=Percent))+labs(title="Glucosinolate Variance")
ggplot(flavVar)+
  geom_col(aes(x=Variable,y=Percent))+labs(title="Flavonol Variance")


table(FatDat$Sample)
flav[flav$flav_Sample=="i27_1",]

gluc[order(gluc$gluc_Sample),]

#Clearly there is no point of doing that. 

```


one<-rep(1,9)
one
t(one)
one %*% one

t(one) %*% one
  

t(one) %*% t(one)

one %*% t(one)

t(one)

growth<-c(12,10,8,11,6,7,2,3,3)  
tannin<-c(0,1,2,3,4,5,6,7,8)  

sum(growth*tannin)


t(growth)%*%tannin

# Sum of products
growth%*%t(tannin)

#sum of squares 

growth %*% growth

tannin %*% tannin

a<-cbind(growth,tannin)
a

t(a)%*%a

  

t(a)
a  




#-------------------------------------------------------------------
GROWTHRATE
#-------------------------------------------------------------------




#########################Growth rate analysis 

#Repeat this analysis on RGR1
```{r}
rm(fit.1,fit.2,fit.3,fit.4,fit.5,fit.5.Whole,fit.6,fitfull,fitfull0,fitfull2)

#reading in growthrate data. 
rgr<-read.csv("GrowthRateIntegrated.csv")

#Merging growth rate data with all other greenhouse data. I am doing this with preference to the growth rate data, which will remove those without a final measurement, therefore excluding those that died, and without glucosinolate data, which would make the analysis impossible anyways. 


GrowthDat<-rgr %>% select("Tag"=Genotype,RGR1,RGR2,RGR3,RGR4)  %>% right_join(dat2,by="Tag")

#Creating dataframe with mean values. 
GrowthDatMean<-rgr %>% select("Tag"=Genotype,RGR1,RGR2,RGR3,RGR4,Family,treatment)  %>% group_by(Family,treatment) %>% summarise_if(is.numeric,mean) 

levels(GrowthDatMean$treatment)<-c("a","gm","m")

GrowthDatMean<-left_join(dat3,GrowthDatMean,by=c("Family","treatment"))


```


#Modelling relative growth rate 1. Pathogens and fern abundance was removed from this analysis because they did not appear this early in the experiment. 
```{r}

#Modelling fixed effects. Ignoring flavonoids for the time being. 

fitfull3<-lmer(RGR1~treatment*gluc_Conc+(1|Family:treatment)+(1|gh_bench/gh_col), data=GrowthDat)


#Removing three way interaction
fit.1<-update(fitfull3, ~.-gluc_Conc:treatment)
anova(fitfull3,fit.1) #Good to remove

fit.2<-update(fit.1,~.-gluc_Conc)
anova(fit.2,fit.1) #Good to remove. 

fit.3<-update(fit.2,~.-treatment)
anova(fit.2,fit.3) #Significant
```
Treament was the only predictor for growth rate 1. 

#Modelling on growth rate 2. 
```{r}
#Modelling fixed effects. Ignoring flavonoids for the time being. 

fitfull3<-lmer(RGR2~treatment*gluc_Conc+(1|Family:treatment)+(1|gh_bench/gh_col), data=GrowthDat)


#Removing three way interaction
fit.1<-update(fitfull3, ~.-gluc_Conc:treatment)
anova(fitfull3,fit.1) #Good to remove

fit.2<-update(fit.1,~.-gluc_Conc)
anova(fit.2,fit.1) #Good to remove. 

fit.3<-update(fit.2,~.-treatment)
anova(fit.2,fit.3) #Significant

```
Treament was the only predictor for growth rate 2.


#Modelling on growth rate 3. 
```{r}
#Modelling fixed effects. Ignoring flavonoids for the time being. 

fitfull3<-lmer(RGR3~treatment*gluc_Conc+Fern+WhiteFungDam+BlackPathDam+ThripsDam+(1|Family:treatment)+(1|gh_bench/gh_col), data=GrowthDat)


#Removing three way interaction
fit.1<-update(fitfull3, ~.-gluc_Conc:treatment)
anova(fitfull3,fit.1) #Not Significant

fit.2<-update(fit.1,~.-gluc_Conc)
anova(fit.2,fit.1) # Not Significant

fit.3<-update(fit.2,~.-treatment)
anova(fit.2,fit.3) #Not Significant

fit.4<-update(fit.3,~.-WhiteFungDam)
anova(fit.3,fit.4) #Not Significant

fit.5<-update(fit.4,~.-BlackPathDam)
anova(fit.4,fit.3) #Not Significant

#Rebuilding model to be the same size as that without thrips dam (missing data)
fit.5<-lmer(RGR3 ~ Fern + ThripsDam + (1 | Family:treatment) + (1 | gh_bench/gh_col),data=GrowthDat[!is.na(GrowthDat$ThripsDam),])

fit.6<-update(fit.5,~.-ThripsDam)
anova(fit.5,fit.6) #Not Significant

fit.7<-update(fit.6,~.-Fern)
anova(fit.6,fit.7) #Not Significant

#No significant predictors of this data. 

```


#Modelling on growth rate 4. This one include the other variables as they were now involved in the experiement around this time. 
```{r}
#Modelling fixed effects. Ignoring flavonoids for the time being. 


fitfull3<-lmer(RGR4~treatment*gluc_Conc+Fern+WhiteFungDam+BlackPathDam+ThripsDam+(1|Family:treatment)+(1|gh_bench/gh_col), data=GrowthDat)


#Removing three way interaction
fit.1<-update(fitfull3, ~.-gluc_Conc:treatment)
anova(fitfull3,fit.1) #This interaction was significant for growth rate 4... interesting.

fit.1<-update(fitfull3,~.-Fern)
anova(fitfull3,fit.1) #Fern is significant

fit.1<-update(fitfull3,~.-BlackPathDam)
anova(fitfull3,fit.1) #Black path dam is Significant

fit.1<-update(fitfull3,~.-WhiteFungDam)
anova(fitfull3,fit.1) #WhiteFungDam not Significant

fit.2<-update(fit.1,~.-ThripsDam)
anova(fit.1,fit.2) #ThripsDam not Significant

#Same conclusion as for total body mass. This means that the trends in total body mass we are seeing likely came about towards the end of the experiment. 

library(lme4)
fit.2<-lmer(RGR4~treatment*gluc_Conc+Fern+BlackPathDam+(1|Family:treatment), data=GrowthDat)

fit.3<-update(fit.2,~.-treatment:gluc_Conc)
anova(fit.2,fit.3) #ThripsDam not Significant


summary(fit.2)
plot(fit.2)
confint(fit.2) #The same trend is observed for RGR4 as for final body mass. 


ggplot(GrowthDat)+
  geom_point(aes(y=RGR4,x=gluc_Conc,colour=treatment))


ggplot(GrowthDatMean)+
  geom_point(aes(y=exp(RGR4),x=gluc_Conc,colour=treatment))+
  geom_abline(intercept=-1.594e-07,slope=1.715e-07,colour="red")+
  geom_abline(intercept=-1.594e-07+2.462e-07,slope=1.715e-07-2.346e-07,colour="green")+
  geom_abline(intercept=-1.594e-07+2.626e-07,slope=1.715e-07-3.034e-07,colour="blue")+theme_classic()

#Running the analysis without the random effects leads to maple treatment having a positive slope....
summary(lm(RGR4~gluc_Conc,data = GrowthDatMean %>% filter(treatment=="m"))) #When you look at means
summary(lm(RGR4~gluc_Conc,data = GrowthDat %>% filter(treatment=="m"))) #But not when you dont. This makes me wonder if my analyses are done properly. 

lm(RGR4~gluc_Conc,data = GrowthDatMean %>% filter(treatment=="a"))





```

















