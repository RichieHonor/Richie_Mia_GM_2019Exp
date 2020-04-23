
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





####Body Size integration garbage. 



#First Remove those without bodydata from the final data frame. 

#Full data set without samples (leaf data) lacking final measurement.
FullDataSetLeafData<-final %>% select(Tag,Sample) %>% filter(!Sample %in% Remainder$Sample[!is.na(Remainder$Sample)])

LeavesWithGenotypeInFinal<-final %>% select(Tag,Sample) %>% filter(Sample %in% Remainder$Sample[!is.na(Remainder$Sample)])

#Full data set without genotypes lacking final measurement.
FullDataSetNoLeafData<-final %>% select(Tag,Sample) %>% filter(!Sample %in% Remainder$Sample[is.na(Remainder$Sample)])

#Determine which samples (leaf data) in the remainder do not have genotypes in the final data set.
Remainder %>% filter(!Tag %in% FullDataSetLeafData$Tag)
#These individuals do not have genotypes in the final data set and the average of the leaf data and glucosinolate data can be taken for these leaves. 

#However, there are 12 leaf which do have body mass and pathogen data available at the level of the plant. These ill will fill in with the average pathogen data available for the plant, as well as the final body mass and number of leaves. 
Remainder %>% filter(Tag %in% UnJoined$Genotype)

#Code Here

#Now all genotypes in the remainder that have bodymass data available have been filled with that data. However, are there still genotypes that we have bodymass data on that are not in the final data frame?



#If the genotype is not in the final data frame, i want to add in the average of the pathogen data for that genotype. 







#Modelling -- Poisson -What influences White pathogen damage?
```{r}
#Rounding to integers for white fung dam. 
dat$WhiteFungDam<-ceiling(dat$WhiteFungDam)

fit_full<-glmer(WhiteFungDam~treatment*gluc_Conc+treatment*flav_Conc+(1|Family/Tag)+(1|gh_bench),family=poisson,data=dat)

fit.1<-update(fit_full,~.-flav_Conc:treatment)
anova(fit.1,fit_full) #Fit1 is a better model


fit.2<-update(fit.1,~.-gluc_Conc:treatment)
anova(fit.2,fit.1) #models are the same... eliminating interaction.

fit.3<-update(fit.2,~.-gluc_Conc)
anova(fit.3,fit.2) #The model including glucConc is better. 

fit.2<-update(fit.2,data=dat[!is.na(dat$flav_Conc),]) #Using a model with a data subset.
fit.3<-update(fit.2,~.-flav_Conc) 
anova(fit.3,fit.2)#flav_Conc is unimportant, even though the AIC and BIC are lower.

#Therefore, the best model is:
summary(glmer(WhiteFungDam~treatment+gluc_Conc+(1|Family/Tag)+(1|gh_bench),family=poisson,data=dat))

#the inclusion of GH_Bench does not affect the outcome of the model, thereofre, I am excluding it as the other model was saturated and took alot time to compute. 
summary(glmer(WhiteFungDam~treatment+gluc_Conc+(1|Family/Tag),family=poisson,data=dat))
#Using a Permutation test to determine the probability of getting a significant glucosinolate concentration predictor, if the data is ramdomly sampled. 

datPerTest<-dat
zStore<-c()
#for(i in 1:1000){
#Randomize glucosinolate concentration.
datPerTest$gluc_Conc<-sample(dat$gluc_Conc,length(dat$gluc_Conc),replace = F)

#Run Model with randomized glucConc and extract p-value
glucZval<-summary(glmer(WhiteFungDam~treatment+gluc_Conc+(1|Family/Tag),family=poisson,data=datPerTest))$coef[4,3]

#Store P value.
zStore[i]<-glucZval
#}

#How many samples are less than or equal to the Z score that i got? (the p value)
sum(zStore<=-3.262)/length(zStore)

#The p value is actually around .10, which is not significant. It was estimated to be 0.001. 

#Given that 64% of the models return a significant p value, I would say that this is a very very very bad model. 
```

#Visuallizing -- white fungal damage is influenced by gluc and treatmnet. 
```{r}

PoisSlope=function(x,int){
  y=exp(-2.32292*x+int)
  return(y)
}

glucplot=seq(min(dat$gluc_Conc),max(dat$gluc_Conc),length.out = 695)

FungA<-PoisSlope(glucplot,-0.15701)
FungM<-PoisSlope(glucplot,-0.15701-0.82060)
FungGM<-PoisSlope(glucplot,-0.15701+0.02483)

fit.g<-glm(WhiteFungDam~gluc_Conc+treatment,family=poisson,data=dat)
summary(fit.g)

predict(fit.g,type="response",newdata=data.frame("treatment"=rep("a",695),"gluc_Conc"=glucplot))

datv<-dat[!is.na(dat$WhiteFungDam),]

#tiff("Defence_Figures/GlucosinolateFungi.tiff", units="in", width=10, height=6, res=300)
ggplot(datv[datv$WhiteFungDam<30,])+
  geom_point(aes(y=WhiteFungDam,x=gluc_Conc,colour=treatment))+
  scale_colour_manual(values=c("#009E73","#56B4E9","#E69F00","#009E73","#56B4E9","#E69F00"),labels=c("Alone","Garlic Mustard","Maple"))+theme_simple()+
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12))+
  geom_path(x=glucplot,y=FungA,colour="#009E73")+
  geom_path(x=glucplot,y=FungM,colour="#E69F00")+
  ylab("Powdery Mildew Damage")+
  xlab(bquote(bold("[Total Glucosinolate] " (mg/ml))))
#dev.off()

```



#Modelling -- What influences Black Pathogen Damage? 
```{r}
dat$BlackPathDam
dat$BlackPathDam<-ceiling(dat$BlackPathDam)
#Glucosinolate and flavonoids are correlated, but the R2 is only 0.39, so i dont think it is a very big deal. 
summary(lm(gluc_Conc~flav_Conc,dat=dat))

fit_full<-glmer(BlackPathDam~treatment*gluc_Conc+treatment*flav_Conc+(1|Family/Tag)+(1|gh_bench),family=poisson,data=dat)

fit.1<-update(fit_full,~.-flav_Conc:treatment)
anova(fit.1,fit_full) #There is a significant interaction between flav concentration and treatment. 

fit.2<-update(fit_full,~.-gluc_Conc:treatment)
anova(fit.2,fit_full) #There is also a significant interaction between gluc Conc and treatment. 

#Diagnostic


plot(fit_full)


#Permutation test. 
datPerTest<-dat
newframe<-NA
for(i in 1:300){
  #Randomize glucosinolate concentration.
  datPerTest$flav_Conc<-sample(dat$flav_Conc,length(dat$flav_Conc),replace = F)
  datPerTest$gluc_Conc<-sample(dat$gluc_Conc,length(dat$gluc_Conc),replace = F)
  datPerTest$treatment<-sample(dat$treatment,length(dat$treatment),replace = F)
  
  #Run Model with randomized glucConc and extract p-value
  ModelZval<-summary(glmer(BlackPathDam~treatment*gluc_Conc+treatment*flav_Conc+(1|Family/Tag),family=poisson,data=datPerTest))$coef[,3]
  
  #Saving all of the Z values in a dataframe. 
  newframe<-rbind(newframe,ModelZval)
}

newframe<-data.frame(newframe)
summary(fit_full)
newframe
summary(fit_full)$coef[1,3]
sum(newframe[,1]<=-1.35,na.rm = T)/length(newframe[,1])

#Flavonoid concentration interaction is not significant. Must be less than 0.25. This is a poor model for this data. 
sum(newframe$flav_Conc<=-5.220,na.rm = T)/length(newframe$flav_Conc)

#This as well is not significant. 
sum(newframe$gluc_Conc>=3.605 ,na.rm = T)/length(newframe$gluc_Conc)

#The p value calculated for this not actually close to zero as in the model, but is around 0.1533333. This is not significant. 
#Approximitely half of the p values are less than 0.05, this suggests that this model is quite horrible. 

```








#Visuallizing -- Black Pathogen damage is influenced by flavonoids. 
```{r}


PoisSlope=function(x){
  y=exp(-1.1150*x+1.6816)
  return(y)
}


flavplot=seq(min(dat$flav_Conc,na.rm = T),max(dat$flav_Conc,na.rm=T),length.out = 687)

flavy<-PoisSlope(flavplot)

#tiff("Defence_Figures/FlavonoidBlackPath.tiff", units="in", width=10, height=6, res=300)
ggplot(dat[!is.na(dat$flav_Conc)&!is.na(dat$flav_Conc),])+
  geom_point(aes(y=BlackPathDam,x=flav_Conc))+theme_simple()+
  geom_path(x=flavplot,y=flavy,size=1,colour="#999999")+
  
  scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+
  ylab("Black Pathogen Damage")+
  xlab(bquote(bold("[Total Flavonoid] " (mg/ml))))
#dev.off()

```







#What influences Thrips damage? 
```{r}

fit_full<-glmer(ThripsDam~treatment+gluc_Conc+flav_Conc+(1|Family/Tag)+(1|gh_bench),family=poisson,data=dat)
fit_full_0.1<-glmer(ThripsDam~treatment+gluc_Conc+flav_Conc+(1|Family/Tag),family=poisson,data=dat)
anova(fit_full,fit_full_0.1) #Bench is an important predictor in this analysis. 

fit.1<-update(fit_full,~.-gluc_Conc)
anova(fit.1,fit_full) #Gluc Conc is a very significant predictor. 


fit.1<-update(fit_full,~.-treatment)
anova(fit.1,fit_full)  #treatment is not a significant predictor. 


fit.full_f<-update(fit.1,data=dat[!is.na(dat$flav_Conc),])#Changing data set to account for missing flavonoid data. 
fit.2<-update(fit.full_f,~.-flav_Conc)
anova(fit.full_f,fit.2) #Flav conc is also a very significant predictor. 


#Now i will test for interactions. 

fit.full_f_int<-update(fit.full_f,~.+flav_Conc:gluc_Conc)
anova(fit.full_f,fit.full_f_int)

#There is a signiicant interaction between gluc conc and flav conc. I am dubious, however, at how robust this conclusion is with the data modelled with a oisson, i will use a permutation test to see. 

summary(fit.full_f_int)

#Permutation test. 
dat$ThripsDam<-ceiling(dat$ThripsDam)
datPerTest<-dat
zStore<-c()
for(i in 1:1000){
  #Randomize glucosinolate concentration.
  datPerTest$flav_Conc<-sample(dat$flav_Conc,length(dat$flav_Conc),replace = F)
  
  #Run Model with randomized glucConc and extract p-value
  flavZval<-summary(glmer(ThripsDam~flav_Conc+(1|Family/Tag),family=poisson,data=datPerTest))$coef[2,3]
  
  #Store P value.
  zStore[i]<-flavZval
}

sum(zStore<=-18.71)/length(zStore)

hist(zStore)
#The p value for this model is still about zero according to the permutation test, which might suggest that it is a better model than the others. 

plot(fit.4)


sum(dat$ThripsDam==0,na.rm = T)/length(dat$ThripsDam)
sum(dat$WhiteFungDam==0,na.rm = T)/length(dat$WhiteFungDam)
sum(dat$BlackPathDam==0,na.rm = T)/length(dat$BlackPathDam)
sum(dat$Fern==0,na.rm = T)/length(dat$Fern)

#42 percent of thrips damage is zero. 
#0.76 % of white fungal damage is zero. 
#60% of black path damage is zero. 
#86% of fern data is zero. 



```





#BlackPathDam -- logistic regression--
```{r}
#Because of how zero inflated white pathogen damage is, i will use a logistic regression to model it.
dat$BlackPathLogis<-NA
dat$BlackPathLogis[dat$BlackPathDam==0]<-0
dat$BlackPathLogis[dat$BlackPathDam>0]<-1

#Models fail to converge with interactions and with gh_bench included in the analysis. 
fit_full_g<-glmer(BlackPathLogis~treatment+gluc_Conc+flav_Conc+(1|Family/Tag),family=binomial,data=dat[!is.na(dat$flav_Conc),])


fit.1<-update(fit_full_g,~.-flav_Conc)
anova(fit_full_g,fit.1) #Flavonoid Concentration is very close to significance (p=0.06 with a AIC lower in the model that included flavonoid concentration)


fit.1<-glmer(BlackPathLogis~treatment+gluc_Conc+(1|Family/Tag),family=binomial,data=dat[!is.na(dat$gluc_Conc),])

fit.2<-update(fit.1,~.-gluc_Conc)
anova(fit.2,fit.1) #Gluc conc is also not significant although it also has a lower AIC value. 

fit.2<-glmer(BlackPathLogis~treatment+(1|Family/Tag),family=binomial,data=dat)
fit.3<-update(fit.2,~.-treatment)
anova(fit.3,fit.2) #Treatment is not a significant predictor. 

#Lets see if flavonoids are significant without glucosinolates in the model, as they are correlated. 
fit.3<-glmer(BlackPathLogis~1+(1|Family/Tag),family=binomial,data=dat[!is.na(dat$flav_Conc),])

fit.4<-update(fit.3,~.+flav_Conc)
anova(fit.3,fit.4) #Without glucosinolates in the equation, flavonoids are very significant. This is the best model.  

summary(fit.4)

```






#Modelling : Negative Binomial White fungal Damage. 
```{r}

dat$WhiteFungDam<-ceiling(dat$WhiteFungDam)


fit_full<-glmmTMB(WhiteFungDam~treatment+gluc_Conc+flav_Conc+(1|Family/Tag)+(1|gh_bench),family=nbinom2,data=dat)

fit_full_0.1<-glmmTMB(WhiteFungDam~treatment+gluc_Conc+flav_Conc+(1|Family/Tag),family=nbinom2,data=dat)

anova(fit_full,fit_full_0.1)#Bench is an important predictor to keep. 

summary(fit_full)

fit_full<-glmmTMB(WhiteFungDam~treatment+gluc_Conc+flav_Conc+(1|Family/Tag)+(1|gh_bench),family=nbinom2,data=dat[!is.na(dat$flav_Conc),])
fit_1<-glmmTMB(WhiteFungDam~treatment+gluc_Conc+(1|Family/Tag)+(1|gh_bench),family=nbinom2,data=dat[!is.na(dat$flav_Conc),])
anova(fit_full,fit_1) #flav_Conc not important. 

fit_1<-glmmTMB(WhiteFungDam~treatment+gluc_Conc+(1|Family/Tag)+(1|gh_bench),family=nbinom2,data=dat[!is.na(dat$gluc_Conc),])
fit_2<-glmmTMB(WhiteFungDam~treatment+(1|Family/Tag)+(1|gh_bench),family=nbinom2,data=dat[!is.na(dat$gluc_Conc),])
anova(fit_1,fit_2) #Gluc_Conc is not important. 


fit_1<-glmmTMB(WhiteFungDam~treatment+(1|Family/Tag)+(1|gh_bench),family=nbinom2,data=dat)
fit_2<-glmmTMB(WhiteFungDam~1+(1|Family/Tag)+(1|gh_bench),family=nbinom2,data=dat)
anova(fit_1,fit_2) #Treatment is also unimportant predictor. (even though it has a lower AIC.)

summary(fit_full)
plot(resid(fit_full))


```












