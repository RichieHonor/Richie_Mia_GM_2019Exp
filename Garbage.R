
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




MAPLE STUFF 



#Redoing model selection with different variable. 
```{r}
#What affects maple growth? secondary compounds, initial leaf size (to control for starting leaf area size), Fern, and maple Stem height (signifies the year of the plant).  


fitFull<-lmer(Maple_TotalLeafArea_End~gluc_Conc*MapleDamage+MapleDamage*flav_Conc+maple_height_0+Fern+maple_leaf_length_0+(1|Family)+(1|gh_bench),data=Maple)

summary(fitFull)

Maple$Family


#Removing interaction
fit.1<-update(fitFull,~.-MapleDamage:flav_Conc)
anova(fitFull,fit.1) #This interaction is NOT significant... removing it. 


fit.2<-update(fit.1,~.-MapleDamage:gluc_Conc)
anova(fit.2,fit.1) #The interaction is also not significant -- Removing this. 

fit.2<-update(fit.2,data=Maple[!is.na(Maple$flav_Conc),])#Using flav_conc subset (accounting for missing values. )
fit.3<-update(fit.2,~.-flav_Conc)
anova(fit.3,fit.2) #Flavonoids are not significant... removing from data set.  

fit.3<-update(fit.3,data=Maple[!is.na(Maple$gluc_Conc),])#Using gluc_conc subset (accounting for missing values. )
fit.4<-update(fit.3,~.-gluc_Conc)
anova(fit.4,fit.3) #The AIC is the same and the BIC is even lower, it is not a significant predictor, but it would belong in the candidate model set. 


fit.5<-update(fit.4,~.-Fern)
anova(fit.5,fit.4) #Fern is not significant. 





fit.6<-update(fit.5,~.-maple_height_0)
anova(fit.6,fit.5)  # Maple stem height is NOT a significant predictor.

fit.7<-update(fit.6,~.-maple_leaf_length_0)
anova(fit.7,fit.6)  # Initial Maple leaf length  is also not a significant predictor

fit.8<-update(fit.7,~.-MapleDamage)
anova(fit.8,fit.7)  # Damage is also a  NOT a significant predictor


#Best model is fit 8. ... Is family a significant predictor? 
fit.9<-update(fit.8,~.-(1|Family))
anova(fit.9,fit.8)  # Family is not a significant predictor. Is bench?


fit.10<-update(fit.8,~.-(1|gh_bench))
anova(fit.10,fit.8)  #Bench isnt either really. 


fit.11<-update(fit.10,~.+(1|gh_bench/gh_col),)
anova(fit.11,fit.8)  #neither is gh column....


summary(fit.8)

```


#If non of the random effects are significant, can i just do a regular linear regression?
```{r}
#Testing over flavonoid data set.... are flavonoids significant?

#Interaction with damage
fit.f1<-lm(Maple_TotalLeafArea_End~MapleDamage*flav_Conc+maple_height_0+maple_leaf_length_0,data=Maple[!is.na(Maple$flav_Conc),])

#no interaction with damage
fit.f2<-lm(Maple_TotalLeafArea_End~MapleDamage+flav_Conc+maple_height_0+maple_leaf_length_0,data=Maple[!is.na(Maple$flav_Conc),])

#no flavonoid at all. 
fit.f3<-lm(Maple_TotalLeafArea_End~flav_Conc+MapleDamage+maple_height_0+Fern+maple_leaf_length_0,data=Maple[!is.na(Maple$flav_Conc),])

AIC(fit.f1,fit.f2,fit.f3)

#The best model is one with flavonoid, this model is just as good as the interaction, meaning the interaction probably isnt a significant predictor. 

#Testing over gluc data set. 
fit.g1<-lm(Maple_TotalLeafArea_End~gluc_Conc*MapleDamage+maple_height_0+maple_leaf_length_0,data=Maple[!is.na(Maple$gluc_Conc),])

fit.g2<-lm(Maple_TotalLeafArea_End~gluc_Conc+MapleDamage+maple_height_0+maple_leaf_length_0,data=Maple[!is.na(Maple$gluc_Conc),])

fit.g3<-lm(Maple_TotalLeafArea_End~MapleDamage+maple_height_0+maple_leaf_length_0,data=Maple[!is.na(Maple$gluc_Conc),])

AIC(fit.g1,fit.g2,fit.g3)

#The model with and without glucosinolates are just as good, meaning the model with glucosinolate is likely not good. 

#Effects of other things... 
fit.o1<-lm(Maple_TotalLeafArea_End~MapleDamage+maple_height_0+maple_leaf_length_0,data=Maple)

fit.o2<-lm(Maple_TotalLeafArea_End~maple_height_0+maple_leaf_length_0,data=Maple)

fit.o3<-lm(Maple_TotalLeafArea_End~MapleDamage+maple_leaf_length_0,data=Maple)

fit.o4<-lm(Maple_TotalLeafArea_End~MapleDamage+maple_height_0,data=Maple)

AIC(fit.o1,fit.o2,fit.o3,fit.o4)

#Model without damage is better and maple stem height and leaf length initial are important. Thereofore the best model is. 


fit.final<-lm(Maple_TotalLeafArea_End~flav_Conc+maple_height_0+maple_leaf_length_0,data=Maple[!is.na(Maple$flav_Conc),])
#Testing the significance of flav concentration again... 
fit.final.2<-lm(Maple_TotalLeafArea_End~maple_height_0+maple_leaf_length_0,data=Maple[!is.na(Maple$flav_Conc),])


AIC(fit.final)
AIC(fit.final.2)

anova(fit.final,fit.final.2)

#Flav conc is not a good predictor. 
summary(fit.final.2)

```
```











#visualization 
```{r}
ggplot(Maple)+
  geom_point(aes(y=Maple_TotalLeafArea_End,x=gluc_Conc,colour=MapleDamage))

ggplot(Maple)+
  geom_point(aes(y=Maple_TotalLeafArea_End,x=flav_Conc,colour=MapleDamage))


#this should not be analyzed as linear because there is a limit at zero. 
ggplot(Maple)+
  geom_point(aes(y=Maple_TotalLeafArea_End,x=Maple_StemHeight_0))


Maple


ggplot(Maple)+
  geom_point(aes(y=Maple_TotalLeafArea_End,x=maple_leaf_length_0))

summary(lm(maple_leaf_length_0~gluc_Conc,data=Maple)) #Adjusted r squared is only 0.1. 


fit.1<-lmer(Maple_TotalLeafArea_End~maple_leaf_length_0+maple_height_0+gluc_Conc*MapleDamage+MapleDamage*flav_Conc+(1|Family)+(1|gh_bench/gh_col),data=Maple)

summary(fit.1)





#Maple effect on fitness. 


```{r}
#There is an odd distribution of leaf size, but it is fairly normal. 
hist(MapleModel$GM_TotalLeaf_Area,breaks=25)

#What effect does maple have on garlic mustard performance? 
MapleModel
#Standardizing maple size, as it is now a predictor variable is has too large of a SD and was messing up the model. 
MapleModel$Maple_TotalLeafArea_End<-Standardize(MapleModel$Maple_TotalLeafArea_End)



#MaxModel (I predict possible interactions between pathogens, because the infection of one pathogen could make others worst and decrease fitness more when infected with more than one pathogen, however i expect competitive effects of fern and maple total leaf area to only be additive because elimiting nutrients would not have any forseeable compounding effects. It should be noted that glucoisnolate concentration is not in the model, as explained above)

#Modelling random effects 
fit<-lmer(GM_TotalLeaf_Area~Maple_TotalLeafArea_End+Fern+ThripsDam*BlackPathDam*WhiteFungLogis+(1|Family)+(1|gh_bench/gh_col),data=MapleModel)

fit2<-lmer(GM_TotalLeaf_Area~Maple_TotalLeafArea_End+Fern+ThripsDam*BlackPathDam*WhiteFungLogis+(1|Family)+(1|gh_bench),data=MapleModel)

fit3<-lmer(GM_TotalLeaf_Area~Maple_TotalLeafArea_End+Fern+ThripsDam*BlackPathDam*WhiteFungLogis+(1|gh_bench/gh_col),data=MapleModel)

#Is collumn important?
anova(fit,fit2) #Yes it is. 

#Is family important? 
anova(fit,fit3) #No, it is just on the verge of significance (p=0.058)


#Modeling fixed effects: 

fit<-lmer(GM_TotalLeaf_Area~Maple_TotalLeafArea_End+Fern+ThripsDam*BlackPathDam*WhiteFungLogis+(1|gh_bench/gh_col),data=MapleModel)

fit2<-update(fit,~.-ThripsDam:BlackPathDam:WhiteFungLogis)
anova(fit,fit2) #Not a significant three way interaction although it is close (p=0.06)

fit3<-update(fit2,~.-BlackPathDam:WhiteFungLogis)
anova(fit3,fit2) #There is not a significant two way interaction 

fit4<-update(fit3,~.-BlackPathDam:ThripsDam)
anova(fit4,fit3)
#There is a significant interaction between black path dam and thrips dam. 

fit5<-update(fit3,~.-WhiteFungLogis:ThripsDam)
anova(fit5,fit3) #There is no interaction between these two. 

summary(fit5)

#is White path Damage significant alone? 
fit6<-update(fit5,~.-WhiteFungLogis)
anova(fit6,fit5) #White Fung logis is not important, although it has an identical AIC.   

summary(fit6)

fit7<-update(fit6,~.-BlackPathDam)
anova(fit7,fit6)

summary(fit7)

#The best model is one with thrips damage and and interaction between thrips and black pathogen damage. Oddly, the coefficent for thrips damage is positive, which seems like thrips damage increases performance, however this is likely just a reflection of there being more thrips damage and black pathogen damage on larger leaves, and of course plants with larger total leaf area will have larger leaves. Therefore, what i need to do is correct for leaf size when using pathogen data and have something like % infection on leaf. spots per leaf area? Otherwise, i am afraid that this model and all pathogen models are not informative. 

plot(ChlorA~GM_TotalLeaf_Area,data=dat2)
summary(lm(ChlorA~GM_TotalLeaf_Area,data=dat2))
#Interestinly, chlrophyllA is not associated with garlic mustard total leaf area -- our measure of performance. bizarre. Could this be because smaller leaves have much more chlorophyll? 

plot(ChlorA~GM_Leaf_Len,data=dat)
dat
#Not really no



```

Fertilizer figures 
#Effect of treatment
ggplot(Final2)+
  geom_boxplot(aes(x=treatment,y=GM_Shoot_Mass,fill=Fert))+theme_simple()+ylab("Shoot Mass (g)")+
  scale_y_continuous(breaks = seq(0,5,0.25))+
  scale_x_discrete(name="",labels=c("Alone","Garlic Mustard","Maple"))+
  scale_fill_manual(values=c("#0072B2","#F0E442"),labels=c("No amendment","Fertilizer"))+
  theme(axis.title.y =  element_text(color = "black", size = 16, face = "bold",margin=margin(3,20,3,0)))
#dev.off()


#Plasticity notes 
#Interstingly there is significant interaction between chlorophyll a and treatment with glucosinolate concentration. What this means exacly is unclear, but the relationship between chlorophyll a and glucosinolate concentration is higher in the maple and garlic mustard treatment. 


#In terms of the fixed effects, the amount of glucosinolates still decreases in the maple and garlic mustard treatment, suggesting that the plasticity we see is actually true plasticity, and is not just occuring because the plants perform worse and so make less glucosinolates, but are actively producing less glucosinolates. This could also just represent a steeper decrease in glucosinolates than in chlorophyll. It would make sense for individuals to reduce the amount of glucosinolate they are producing before they reduce the amount of chlorophyll. This is in line with what is found for thrips damage, for some reason, glucosinolate concentration is reduced with increasing thrips damage after controlling for performance. This does not change, whether thrips damage is weighed by leaf size or not -- the effect is still negative. 


#What predicts flavonoid concentration? (Apparent Plasticity)
```{r}
hist(dat$flav_Conc) 
#The distribution looks faily normal, so i am going to use a linear mixed model. # My maximum model is one with possible interactions between pathogens, but no other interactions

#------------------
#Modelling Random Effects 

#Maximum model
fit<-lmer(flav_Conc~treatment+WhiteFungLogis*BlackPathDamW*ThripsDamW+Fern+GM_Leaf_Area+I(GM_Leaf_Area^2)+(1|Family/Tag)+(1|gh_bench),data=dat)

#Is Tag important?
fit2<-lmer(flav_Conc~treatment+WhiteFungLogis*BlackPathDamW*ThripsDamW+Fern+GM_Leaf_Area+I(GM_Leaf_Area^2)+(1|Family)+(1|gh_bench),data=dat)
anova(fit,fit2)#The model with Tag is not significant,although p=0.06, thereofer i will leave it in the model to avoid psuedoreplication

#Is greenhouse bench important?
fit3<-lmer(flav_Conc~treatment+WhiteFungLogis*BlackPathDamW*ThripsDamW+Fern+GM_Leaf_Area+I(GM_Leaf_Area^2)+(1|Family/Tag),data=dat)
anova(fit,fit3) #gh_Bench is very significant. Keeping this in the analysis 

#Is family important? (Genotype effects)
fit4<-lmer(flav_Conc~treatment+WhiteFungLogis*BlackPathDamW*ThripsDamW+Fern+GM_Leaf_Area+I(GM_Leaf_Area^2)+(1|Tag)+(1|gh_bench),data=dat)
anova(fit,fit4) #Family is not important. 

#Is there a GXE interaction? 
fit5<-lmer(flav_Conc~treatment+WhiteFungLogis*BlackPathDamW*ThripsDamW+Fern+GM_Leaf_Area+I(GM_Leaf_Area^2)+(1|Family:treatment)+(1|gh_bench),data=dat)
anova(fit2,fit5) #There is no GXE interaction

#Double check family is not important because there was no significant effect of Tag. 
fit6<-lmer(flav_Conc~treatment+WhiteFungLogis*BlackPathDamW*ThripsDamW+Fern+GM_Leaf_Area+I(GM_Leaf_Area^2)+(1|gh_bench),data=dat)
anova(fit2,fit6) #Family is not important. 

#There is no significant family within family genetic variation for flavonoid concentration. 

#------------------
#Modelling fixed effects. 

#Maximum model:
fit<-lmer(flav_Conc~treatment+WhiteFungLogis*ThripsDamW*BlackPathDamW+Fern+GM_Leaf_Area+I(GM_Leaf_Area^2)+(1|Tag)+(1|gh_bench),data=dat)

fit2<-update(fit,~.-WhiteFungLogis:ThripsDamW:BlackPathDamW)
anova(fit,fit2)

fit3<-update(fit2,~.-WhiteFungLogis:ThripsDamW)
anova(fit2,fit3) # There is a significant interaction between white fungi presense and thripsdamage. 

fit4<-update(fit2,~.-ThripsDamW:BlackPathDamW)
anova(fit2,fit4) #no relationship here. 

fit5<-update(fit4,~.-WhiteFungLogis:BlackPathDamW)
anova(fit5,fit4)#There are no significant interactions. 


fit4<-lmer(flav_Conc~treatment+BlackPathDamW+WhiteFungLogis+ThripsDamW+Fern+GM_Leaf_Area+WhiteFungLogis:ThripsDamW+I(GM_Leaf_Area^2)+(1|gh_bench),data=dat)
summary(fit4)
#Removing black Path Dam

fit4<-lmer(flav_Conc~treatment+WhiteFungLogis+ThripsDamW+Fern+GM_Leaf_Area+WhiteFungLogis:ThripsDamW+I(GM_Leaf_Area^2)+(1|gh_bench),data=dat)
summary(fit4)#Revoving Fern (although it is close to significance (p=0.054))

fit5<-lmer(flav_Conc~treatment+WhiteFungLogis+ThripsDamW+GM_Leaf_Area+WhiteFungLogis:ThripsDamW+I(GM_Leaf_Area^2)+(1|gh_bench),data=dat)
summary(fit5) #Best Model


summary(fit5)#This is the best model. 
plot(fit5)
qqnorm(residuals(fit5))
plot(residuals(fit5))

#Treatment (Maple and garlic mustard) and an interaction between whitefungal abundance and thripsdamage reduced flavonoid concentration (apparent plasticity)
```


Modelling Random Effects 

#Maximum model
fit<-lmer(flav_Conc~treatment+WhiteFungLogis*BlackPathDamW*ThripsDamW+Fern+GM_Leaf_Area+I(GM_Leaf_Area^2)+ChlorA+treatment:ChlorA+(1|Family/Tag)+(1|gh_bench),data=dat)

#Is Tag important?
fit2<-lmer(flav_Conc~treatment+WhiteFungLogis*BlackPathDamW*ThripsDamW+Fern+GM_Leaf_Area+I(GM_Leaf_Area^2)+ChlorA+treatment:ChlorA+(1|Family)+(1|gh_bench),data=dat)
anova(fit,fit2)#Tag is unimportant (i.e. there is no significant genetic variation within an individual.. leaves are very different) 

#Is greenhouse bench important?
fit3<-lmer(flav_Conc~treatment+WhiteFungLogis*BlackPathDamW*ThripsDamW+Fern+GM_Leaf_Area+I(GM_Leaf_Area^2)+ChlorA+treatment:ChlorA+(1|Family/Tag),data=dat)
anova(fit,fit3) #greenhouse bench is very significant. I will keep this.

#Is family important? (Genotype effects)
fit4<-lmer(flav_Conc~treatment+WhiteFungLogis*BlackPathDamW*ThripsDamW+Fern+GM_Leaf_Area+I(GM_Leaf_Area^2)+ChlorA+treatment:ChlorA+(1|Tag)+(1|gh_bench),data=dat)
anova(fit,fit4) #Family is not important.

#Is there a GXE interaction? 
fit5<-lmer(flav_Conc~treatment+WhiteFungLogis*BlackPathDamW*ThripsDamW+Fern+GM_Leaf_Area+I(GM_Leaf_Area^2)+ChlorA+treatment:ChlorA+(1|Family:treatment)+(1|gh_bench),data=dat)
anova(fit2,fit5) #There is no GXE interaction

#Double check that family isnt important because Tag was not significant.
fit6<-lmer(flav_Conc~treatment+WhiteFungLogis*BlackPathDamW*ThripsDamW+Fern+GM_Leaf_Area+I(GM_Leaf_Area^2)+ChlorA+treatment:ChlorA+(1|gh_bench),data=dat)
anova(fit6,fit2) #Family is not important, although it was close to significant (p=0.068)

#Modelling Fixed Effects. 

#Maximum model
fit<-lmer(flav_Conc~treatment+WhiteFungLogis*BlackPathDamW*ThripsDamW+Fern+GM_Leaf_Area+I(GM_Leaf_Area^2)+ChlorA+treatment:ChlorA+(1|gh_bench),data=dat)

fit1<-update(fit,~.- WhiteFungLogis:BlackPathDamW:ThripsDamW)
anova(fit,fit1) #No three way interaction

fit2<-update(fit1,~.- WhiteFungLogis:BlackPathDamW)
anova(fit2,fit1)#No interaction

fit3<-update(fit2,~.- WhiteFungLogis:ThripsDamW)
anova(fit3,fit2) #No interaction

fit4<-update(fit3,~.- BlackPathDamW:ThripsDamW)
anova(fit4,fit3) #No interaction

summary(fit4) #Polynomial of leaf size isn't significant. 

fit4<-lmer(flav_Conc~treatment+WhiteFungLogis+BlackPathDamW+ThripsDamW+Fern+GM_Leaf_Area+ChlorA+treatment:ChlorA+(1|gh_bench),data=dat)

fit5<-update(fit4,~.- ChlorA:treatment)
anova(fit5,fit4) #the chlorophyll A interaction is not important. 

fit5<-lmer(flav_Conc~treatment+WhiteFungLogis+BlackPathDamW+ThripsDamW+Fern+GM_Leaf_Area+ChlorA+(1|gh_bench),data=dat)
summary(fit5) #Fern not significant... but just barely (p=0.055)

fit6<-lmer(flav_Conc~treatment+WhiteFungLogis+BlackPathDamW+ThripsDamW+GM_Leaf_Area+ChlorA+(1|gh_bench),data=dat)
summary(fit6) #WhiteFung Dam is not significant


fit7<-lmer(flav_Conc~treatment+BlackPathDamW+ThripsDamW+GM_Leaf_Area+ChlorA+(1|gh_bench),data=dat)
summary(fit7)#Black Path Dam sint significant. 

fit7<-lmer(flav_Conc~treatment+ThripsDamW+GM_Leaf_Area+ChlorA+(1|gh_bench),data=dat)

#is treatment significant? 
fit8<-update(fit7,~.-treatment)
anova(fit7,fit8) #Treatment is highly significant. 

#Therefore, the best model is fit7, with predictors of thrips damage, black pathogen damage, treatment and leaf area all affecting the true plasticity of flavonoid compounds in this experiement. 

#However, the surprising result is that the amount of flavonoids actually decreases with increasing thrips and blackpathogen damage, even after chlorophyll is accounted for. This suggests that when the plant is suffering some sort of fitness consequence, it reduces the amount of glucosinolates and flavonoids before it reduces the amount of chlorophyll, which makes sense. 

plot(fit7)
qqnorm(residuals(fit7))
plot(residuals(fit7))


```
Even aft



#What predicts flavonoid concentration after controlling for performance? (True Plasticity). 
```{r}
#Modelling Random Effects 

#Maximum model
fit<-lmer(flav_Conc~treatment+WhiteFungLogis*BlackPathDamW*ThripsDamW+Fern+GM_Leaf_Area+I(GM_Leaf_Area^2)+ChlorA+treatment:ChlorA+(1|Family/Tag)+(1|gh_bench),data=dat)

#Is Tag important?
fit2<-lmer(flav_Conc~treatment+WhiteFungLogis*BlackPathDamW*ThripsDamW+Fern+GM_Leaf_Area+I(GM_Leaf_Area^2)+ChlorA+treatment:ChlorA+(1|Family)+(1|gh_bench),data=dat)
anova(fit,fit2)#Tag is unimportant (i.e. there is no significant genetic variation within an individual.. leaves are very different) 

#Is greenhouse bench important?
fit3<-lmer(flav_Conc~treatment+WhiteFungLogis*BlackPathDamW*ThripsDamW+Fern+GM_Leaf_Area+I(GM_Leaf_Area^2)+ChlorA+treatment:ChlorA+(1|Family/Tag),data=dat)
anova(fit,fit3) #greenhouse bench is very significant. I will keep this.

#Is family important? (Genotype effects)
fit4<-lmer(flav_Conc~treatment+WhiteFungLogis*BlackPathDamW*ThripsDamW+Fern+GM_Leaf_Area+I(GM_Leaf_Area^2)+ChlorA+treatment:ChlorA+(1|Tag)+(1|gh_bench),data=dat)
anova(fit,fit4) #Family is not important.

#Is there a GXE interaction? 
fit5<-lmer(flav_Conc~treatment+WhiteFungLogis*BlackPathDamW*ThripsDamW+Fern+GM_Leaf_Area+I(GM_Leaf_Area^2)+ChlorA+treatment:ChlorA+(1|Family:treatment)+(1|gh_bench),data=dat)
anova(fit2,fit5) #There is no GXE interaction

#Double check that family isnt important because Tag was not significant.
fit6<-lmer(flav_Conc~treatment+WhiteFungLogis*BlackPathDamW*ThripsDamW+Fern+GM_Leaf_Area+I(GM_Leaf_Area^2)+ChlorA+treatment:ChlorA+(1|gh_bench),data=dat)
anova(fit6,fit2) #Family is not important, although it was close to significant (p=0.068)

#Modelling Fixed Effects. 

#Maximum model
fit<-lmer(flav_Conc~treatment+WhiteFungLogis*BlackPathDamW*ThripsDamW+Fern+GM_Leaf_Area+I(GM_Leaf_Area^2)+ChlorA+treatment:ChlorA+(1|gh_bench),data=dat)

fit1<-update(fit,~.- WhiteFungLogis:BlackPathDamW:ThripsDamW)
anova(fit,fit1) #No three way interaction

fit2<-update(fit1,~.- WhiteFungLogis:BlackPathDamW)
anova(fit2,fit1)#No interaction

fit3<-update(fit2,~.- WhiteFungLogis:ThripsDamW)
anova(fit3,fit2) #No interaction

fit4<-update(fit3,~.- BlackPathDamW:ThripsDamW)
anova(fit4,fit3) #No interaction

summary(fit4) #Polynomial of leaf size isn't significant. 

fit4<-lmer(flav_Conc~treatment+WhiteFungLogis+BlackPathDamW+ThripsDamW+Fern+GM_Leaf_Area+ChlorA+treatment:ChlorA+(1|gh_bench),data=dat)

fit5<-update(fit4,~.- ChlorA:treatment)
anova(fit5,fit4) #the chlorophyll A interaction is not important. 

fit5<-lmer(flav_Conc~treatment+WhiteFungLogis+BlackPathDamW+ThripsDamW+Fern+GM_Leaf_Area+ChlorA+(1|gh_bench),data=dat)
summary(fit5) #Fern not significant... but just barely (p=0.055)

fit6<-lmer(flav_Conc~treatment+WhiteFungLogis+BlackPathDamW+ThripsDamW+GM_Leaf_Area+ChlorA+(1|gh_bench),data=dat)
summary(fit6) #WhiteFung Dam is not significant


fit7<-lmer(flav_Conc~treatment+BlackPathDamW+ThripsDamW+GM_Leaf_Area+ChlorA+(1|gh_bench),data=dat)
summary(fit7)#Black Path Dam sint significant. 

fit7<-lmer(flav_Conc~treatment+ThripsDamW+GM_Leaf_Area+ChlorA+(1|gh_bench),data=dat)

#is treatment significant? 
fit8<-update(fit7,~.-treatment)
anova(fit7,fit8) #Treatment is highly significant. 

#Therefore, the best model is fit7, with predictors of thrips damage, black pathogen damage, treatment and leaf area all affecting the true plasticity of flavonoid compounds in this experiement. 

#However, the surprising result is that the amount of flavonoids actually decreases with increasing thrips and blackpathogen damage, even after chlorophyll is accounted for. This suggests that when the plant is suffering some sort of fitness consequence, it reduces the amount of glucosinolates and flavonoids before it reduces the amount of chlorophyll, which makes sense. 

plot(fit7)
qqnorm(residuals(fit7))
plot(residuals(fit7))


```



#Is there plasticity to pathogens? 
```{r}

fit<-lm(gluc_Conc~treatment+GM_Leaf_Area+I(GM_Leaf_Area^2)+ThripsDamW*Family,data=dat)
fit2<-lm(gluc_Conc~treatment+GM_Leaf_Area+I(GM_Leaf_Area^2)+ThripsDamW+Family,data=dat)
summary(fit)

anova(fit,fit2)#there is not genetic variation for plasticity...but this isnt a very good framework. 
```








