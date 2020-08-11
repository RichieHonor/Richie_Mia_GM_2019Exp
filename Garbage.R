
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









#What predicts performance in the intraspecific treatment?
```{r}

#Modelling random effects 
fit<-lmer(GM_TotalLeaf_Area~Fern+ThripsDamW*BlackPathDamW*WhiteFungLogis+(1|Family)+(1|gh_bench/gh_col)+(1|competitor),data=dat2)

fit2<-lmer(GM_TotalLeaf_Area~Fern+ThripsDamW*BlackPathDamW*WhiteFungLogis+(1|Family)+(1|gh_bench)+(1|competitor),data=dat2)

fit3<-lmer(GM_TotalLeaf_Area~Fern+ThripsDamW*BlackPathDamW*WhiteFungLogis+(1|gh_bench/gh_col)+(1|competitor),data=dat2)

fit4<-lmer(GM_TotalLeaf_Area~Fern+ThripsDamW*BlackPathDamW*WhiteFungLogis+(1|Family)+(1|competitor),data=dat2)

fit5<-lmer(GM_TotalLeaf_Area~Fern+ThripsDamW*BlackPathDamW*WhiteFungLogis+(1|Family),data=dat2)

#Is collumn important?
anova(fit,fit2) #no it is not.

#Is family important? 
anova(fit,fit3) #Yes it is. 

#Is gh_bench important? 
anova(fit2,fit4) #No it is not. 

#Is competitor important? 
anova(fit4,fit5) #No it is not. So family is imporant, but the competitor is not important. Intersting. 

#Modelling fixed effects 
fit<-lmer(GM_TotalLeaf_Area~Fern+ThripsDamW*BlackPathDamW*WhiteFungLogis+(1|Family),data=dat2)

fit2<-update(fit,~.-ThripsDamW:BlackPathDamW:WhiteFungLogis)
anova(fit,fit2) #There is a significant three way interaction.
summary(fit)

#Removing fern as it is not significant. 
fit2<-lmer(GM_TotalLeaf_Area~ThripsDamW*BlackPathDamW*WhiteFungLogis+(1|Family),data=dat2)

#Rechecking that the three way interaction is significant. 
fit2.5<-update(fit2,~.-ThripsDamW:BlackPathDamW:WhiteFungLogis)
anova(fit2,fit2.5) #it is not still significant... but is close enough i think . 0.055

anova(fit2.5,fit)

summary(fit2) 
```






#What predicts garlic mustard performance in the alone treatment? 
```{r}
#Weighted pathogen damage data are no on too small of a scale. These will be standardized. 


#Modelling random effects 
fit<-lmer(GM_TotalLeaf_Area~Fern+ThripsDamW*BlackPathDamW*WhiteFungLogis+(1|Family)+(1|gh_bench/gh_col),data=dat2)

fit2<-lmer(GM_TotalLeaf_Area~Fern+ThripsDamW*BlackPathDamW*WhiteFungLogis+(1|Family)+(1|gh_bench),data=dat2)

fit3<-lmer(GM_TotalLeaf_Area~Fern+ThripsDamW*BlackPathDamW*WhiteFungLogis+(1|gh_bench/gh_col),data=dat2)

fit4<-lmer(GM_TotalLeaf_Area~Fern+ThripsDamW*BlackPathDamW*WhiteFungLogis+(1|Family),data=dat2)

fit5<-lmer(GM_TotalLeaf_Area~Fern+ThripsDamW*BlackPathDamW*WhiteFungLogis+(1|gh_bench),data=dat2)

#Is collumn important?
anova(fit,fit2) #no it is not.

#Is family important? 
anova(fit,fit5) #nNot at all. 

#Is gh_bench important? 
anova(fit2,fit4) #yes it is



#Modeling fixed effects: 

fit<-lmer(GM_TotalLeaf_Area~Fern+ThripsDamW*BlackPathDamW*WhiteFungLogis+(1|gh_bench),data=dat2)

fit2<-update(fit,~.-ThripsDamW:BlackPathDamW:WhiteFungLogis)
anova(fit,fit2) #Not a significant three way interaction.

fit3<-update(fit2,~.-BlackPathDamW:WhiteFungLogis)
anova(fit3,fit2) #There is not a significant two way interaction, 

fit4<-update(fit3,~.-BlackPathDamW:ThripsDamW)
anova(fit4,fit3)
#There is not a significant interaction between black path dam and thrips dam. 

fit5<-update(fit4,~.-WhiteFungLogis:ThripsDamW)
anova(fit5,fit4) #There is not a significant whitefung dam and thrips dam interaction. 

summary(fit5)#Fern is not significant
fit6<-lmer(GM_TotalLeaf_Area~ThripsDamW+BlackPathDamW+WhiteFungLogis+(1|gh_bench),data=dat2)
summary(fit6)#thrips damage is not significant. 

fit7<-lmer(GM_TotalLeaf_Area~BlackPathDamW+WhiteFungLogis+(1|gh_bench),data=dat2)
summary(fit7)#White fungal damage and black pathogen damage are significant. 



#These coeficients make much more sense, all in the negatives. Weighing was the appropriate thing to do. This should be done in models that predict pathogen response to glucosinolates and flavonoids. - Chlorophyll concentration should be used in those models as well. 
summary(fit7)
hist(residuals(fit7))
plot(residuals(fit7))
plot(fit7)
```


#Effect of maple on GM performance. 
I chose to seperate this model from the last due to the causality inherent in this data. Glucosinolate concentration should affect garlic mustard fitness only through effects on maple and pathogen abundance. Including glucosinolate concentration in a model with both of these variables to determien the effect of glucosinolates on fitness is actually determining how much glucosinolates affect fitness when these terms are controlled for, which should be not at all. Therefore multiple models are needed to test if glucosinolates influence abundance of pathogens and maple performance, after controlling for garlic mustard performance (Using chlorophyll A, because those with high performance might simply be making more glucosinolates). Then i can see how much maple and pathogens influence garlic mustard performance. 
```{r}
#ThripsDam and BlackPathDam represent standardized, logged data which are weighed by leaf size. 


#Modelling random effects 
fit<-lmer(GM_TotalLeaf_Area~Maple_TotalLeafArea_End+log(Fern+1)+ThripsDamW*BlackPathDamW*WhiteFungLogis+(1|Family)+(1|gh_bench/gh_col),data=MapleModel)

fit2<-lmer(GM_TotalLeaf_Area~Maple_TotalLeafArea_End+log(Fern+1)+ThripsDamW*BlackPathDamW*WhiteFungLogis+(1|Family)+(1|gh_bench),data=MapleModel)

fit3<-lmer(GM_TotalLeaf_Area~Maple_TotalLeafArea_End+log(Fern+1)+ThripsDamW*BlackPathDamW*WhiteFungLogis+(1|gh_bench/gh_col),data=MapleModel)


#Is collumn important?
anova(fit,fit2) #Yes it is. 

#Is family important? 
anova(fit,fit3) #Yes family is very important, this is a different result from the one I received prior to weighing the pathogen data. 


#Modeling fixed effects: 

fit<-lmer(GM_TotalLeaf_Area~Maple_TotalLeafArea_End+log(Fern+1)+ThripsDamW*BlackPathDamW*WhiteFungLogis+(1|Family)+(1|gh_bench/gh_col),data=MapleModel)

fit2<-update(fit,~.-ThripsDamW:BlackPathDamW:WhiteFungLogis)
anova(fit,fit2) #Not a significant three way interaction.

fit3<-update(fit2,~.-BlackPathDamW:WhiteFungLogis)
anova(fit3,fit2) #There is not a significant two way interaction.

fit4<-update(fit3,~.-BlackPathDamW:ThripsDamW)
anova(fit4,fit3)
#There is a significant interaction between black path dam and thrips dam interaction p=0.08

fit5<-update(fit3,~.-WhiteFungLogis:ThripsDamW)
anova(fit5,fit3) #There is not significant white path and thripsdamage interaction.

summary(fit5) #White fungal damage not significant. 


fit6<-lmer(GM_TotalLeaf_Area ~ Maple_TotalLeafArea_End + log(Fern + 1) +  
             ThripsDamW + BlackPathDamW  + (1 | Family) +  
             (1 | gh_bench/gh_col) + ThripsDamW:BlackPathDamW,data= MapleModel)

summary(fit6)
plot(fit2)
plot(residuals(fit2))
qqnorm(residuals(fit2)) #looks good. 






#Modelling with full data set. 
```{r}
dat2<-dat2[dat2$treatment=="a",]


#Modelling random effects 
fit<-lmer(GM_TotalLeaf_Area~Fern+ThripsDamW*BlackPathDamW*WhiteFungLogis+(1|Family)+(1|gh_bench/gh_col),data=dat2)

fit2<-lmer(GM_TotalLeaf_Area~Fern+ThripsDamW*BlackPathDamW*WhiteFungLogis+(1|Family)+(1|gh_bench),data=dat2)

fit3<-lmer(GM_TotalLeaf_Area~Fern+ThripsDamW*BlackPathDamW*WhiteFungLogis+(1|gh_bench/gh_col),data=dat2)

fit4<-lmer(GM_TotalLeaf_Area~Fern+ThripsDamW*BlackPathDamW*WhiteFungLogis+(1|Family),data=dat2)


#Is collumn important?
anova(fit,fit2) #no it is not.

#Is family important? 
anova(fit,fit3) #almost, but no it is not

#Is gh_bench important? 
anova(fit2,fit4) #yes it is



#Modeling fixed effects: 

fit<-lmer(GM_TotalLeaf_Area~Fern+ThripsDamW*BlackPathDamW*WhiteFungDam+(1|gh_bench),data=dat2)

fit2<-update(fit,~.-ThripsDamW:BlackPathDamW:WhiteFungDam)
anova(fit,fit2) #Not a significant three way interaction.

fit3<-update(fit2,~.-BlackPathDamW:WhiteFungDam)
anova(fit3,fit2) #There is not a significant two way interaction, 

fit4<-update(fit3,~.-BlackPathDamW:ThripsDamW)
anova(fit4,fit3)
#There is not a significant interaction between black path dam and thrips dam. 

fit5<-update(fit4,~.-WhiteFungDam:ThripsDamW)
anova(fit5,fit4) #There is not a significant whitefung dam and thrips dam interaction. 

summary(fit5)#White Fung Dam is not significant.

fit6<-lmer(GM_TotalLeaf_Area~ThripsDamW+BlackPathDamW+Fern+(1|gh_bench),data=dat2)
summary(fit6)
```








#Modelling: What influences performance (total leaf area)
```{r}
#This is the full model.I expect that the influence of gluc conc and flav conc could vary between treatments because of allelopathy, therefore, i include interactions with these variables. However, i have no reason to think that pathogens or fern abundance could influence fitness differently in different treatments, therefore, no interactions are included.

#First lets ensure we are using the correct random effects.


fitfull<-lmer(GM_TotalLeaf_Area~treatment*gluc_Conc*flav_Conc+BlackPathDam+WhiteFungDam+ThripsDam+Fern+(1|Family)+(1|gh_bench/gh_row), data=dat2)

fitfull2<-lmer(GM_TotalLeaf_Area~treatment*gluc_Conc*flav_Conc+BlackPathDam+WhiteFungDam+ThripsDam+Fern+(1|Family)+(1|gh_bench/gh_col), data=dat2)

fitfull3<-lmer(GM_TotalLeaf_Area~treatment*gluc_Conc*flav_Conc+BlackPathDam+WhiteFungDam+ThripsDam+Fern+(1|Family)+(1|gh_bench), data=dat2)

fitfull4<-lmer(GM_TotalLeaf_Area~treatment*gluc_Conc*flav_Conc+BlackPathDam+WhiteFungDam+ThripsDam+Fern+(1|Family), data=dat2)

fitfull5<-lmer(GM_TotalLeaf_Area~treatment*gluc_Conc*flav_Conc+BlackPathDam+WhiteFungDam+ThripsDam+Fern+(1|gh_bench), data=dat2)

fitfull6<-lmer(GM_TotalLeaf_Area~treatment*gluc_Conc*flav_Conc+BlackPathDam+WhiteFungDam+ThripsDam+Fern+(1|gh_bench/gh_col), data=dat2)

anova(fitfull0,fitfull2,fitfull3,fitfull4,fitfull5,fitfull6) #Best model is fitfull6 and fitfull2. I will compare these two directly
anova(fitfull6,fitfull2) #Family is not a significant predictor and will be excluded, but greenhouse location is and will be included. 

#Modelling fixed effects. -----------

#Full Model
fitfull3<-lmer(GM_TotalLeaf_Area~treatment*gluc_Conc*flav_Conc+BlackPathDam+WhiteFungDam+ThripsDam+Fern+(1|gh_bench/gh_col), data=dat2)

#Removing three way interaction
fit.1<-update(fitfull3, ~.-treatment:gluc_Conc:flav_Conc)
anova(fitfull3,fit.1) #Good to remove

#Removing two way interaction gluc:flav
fit.2<-update(fit.1,~.-gluc_Conc:flav_Conc)
anova(fit.2,fit.1) #Good to remove. 

#Removing treatment:flavonoid interactions
fit.3<-update(fit.2,~.-treatment:flav_Conc)
anova(fit.2,fit.3) #Good to remove. 


#Flavonoid Concentration has missing samples. Therefore, I fit two models, one including flav_Conc and one that does not; however, both use the same limited dataset.
#No flav_Conc
fit.4<-lmer(GM_TotalLeaf_Area ~ treatment + gluc_Conc + BlackPathDam + WhiteFungDam +  
              ThripsDam + Fern +   treatment:gluc_Conc + (1 | gh_bench/gh_col)  ,data=dat2[!is.na(dat2$flav_Conc),])
#Yes flav_Conc
fit.3<-lmer(GM_TotalLeaf_Area ~ treatment + gluc_Conc + BlackPathDam + WhiteFungDam +  
              ThripsDam + Fern +   treatment:gluc_Conc+flav_Conc + (1 | gh_bench/gh_col)  ,data=dat2[!is.na(dat2$flav_Conc),])

anova(fit.3,fit.4) #Flav_Conc is a significant predictor. 

#Now that I know flav_conc is significant, I will conduct the rest of the model simplification using the whole data set (i.e. flav_Conc excluded --fit4) 



fit.nf<-lmer(GM_TotalLeaf_Area ~ treatment + gluc_Conc + BlackPathDam + WhiteFungDam +  
               ThripsDam + Fern +   treatment:gluc_Conc+ (1 | gh_bench/gh_col)  ,data=dat2)

fit.0.nf<-update(fit.nf,~.-treatment:gluc_Conc)
anova(fit.0.nf,fit.nf)  #The interaction is right on the verge of significance, I will keep it included. 

#Creating a subsetted dataframe for the black pathogen damage. 
fit.nf.bp<-lmer(GM_TotalLeaf_Area ~ treatment + gluc_Conc + BlackPathDam + WhiteFungDam +  
                  ThripsDam + Fern +   treatment:gluc_Conc+ (1 | gh_bench/gh_col)  ,data=dat2[!is.na(dat2$BlackPathDam),])

fit.1.nf<-update(fit.nf.bp,~.-BlackPathDam)
anova(fit.nf.bp,fit.1.nf) #Black pathogen damage is a very significant predictor. Did not remove.

fit.2.nf<-update(fit.nf,~.-WhiteFungDam)
anova(fit.nf,fit.2.nf)  #White Fungal Damage was not significant, however, both the AIC and BIC was lower with it included, suggesting it may be an important variable. 


#Refitting a model for thrips dam 
fit.nf<-lmer(GM_TotalLeaf_Area ~ treatment + gluc_Conc + BlackPathDam +  
               ThripsDam + Fern +   treatment:gluc_Conc+ (1 | gh_bench/gh_col)  ,data=dat2[!is.na(dat2$ThripsDam),])

fit.3.nf<-update(fit.nf,~.-ThripsDam)
anova(fit.nf,fit.3.nf)  #Thrips Damage was also not a significant predictor, however, both the AIC and BIC was lower with it included, suggesting it may be an important predictor. 


fit.nf<-lmer(GM_TotalLeaf_Area ~ treatment + gluc_Conc + BlackPathDam +  
               ThripsDam +  Fern+  treatment:gluc_Conc+ (1 | gh_bench/gh_col)  ,data=dat2)

fit.4.nf<-update(fit.nf,~.-Fern)
anova(fit.nf,fit.4.nf)  #Fern abundance is significnat predictor performance.


library(lmerTest)
#The best model is therefore: 

#With limited dataset containing flavonoid data
summary(lmer(GM_TotalLeaf_Area ~ treatment*gluc_Conc+Fern+BlackPathDam +flav_Conc + (1 | gh_bench/gh_col) ,data=dat2))

#with whole data set lacking flavonoid data
summary(lmer(GM_TotalLeaf_Area ~ treatment*gluc_Conc+Fern+BlackPathDam  + (1 | gh_bench/gh_col) ,data=dat2))

#The negative slope associated with the alone treamtment is no longer significant, suggesting that there is only a postive slope of glucosinolates in the maple treatment. However, the issue with this analysis is that is does not account for the competitive strength of the maple competitor, which may be correlated with glucosinolate concentration. 

#Models with Thrips damage and white pathogen damage exhibited lower AIC and BIC values than those without, however, the likelyhood ratio test was not significant. Can I use poisson data as a predictor in this model? Should i log transform it? I am not sure. 



```


Is there an effect of genotype on growth rate? Variance Component Analysis
```{r}

summary(fit1.1)

100*c(6.411059e-08,3.963266e-07)/sum(c(6.411059e-08,3.963266e-07))



#Genotype only accounts for 13% of the variation, wheras there is still 86% of the variation in the residuals. I do not think that genotype is even a worthwhile thing to account for in these models. 

#Interestingly, there appears to be a GXE interaction when testing it out in a linear model framework. 
fit.lm.GXE<-lm(RGR1~treatment*Family,data=datLong)
fit.lm.GandE<-lm(RGR1~treatment+Family,data=datLong)

fit.lm.E<-lm(RGR1~treatment,data=datLong) 
fit.lm.G<-lm(RGR1~Family,data=datLong)

summary(fit.lm.E)
summary(fit.lm.GXE)
summary(fit.lm.GandE)

AIC(fit.lm.E)#This is the best model
AIC(fit.lm.GandE)
AIC(fit.lm.GXE)
AIC(fit.lm.G) #This is the worst model. 



summary(fit.lm.G) 
anova(fit.lm.GXE,fit.lm.E)
anova(fit.lm.GXE,fit.lm.GandE)
anova(fit.lm.GXE,fit.lm.G)



```





#Simply analysis of relative growthrate to determine if there is evidence for treatment effects or genetic variation for growth rate. 
```{r}


first<-dat2C %>% group_by(Genotype) %>% arrange(Seconds) %>% summarize(Seconds = first(Seconds),area=first(area),treatment=first(treatment),Family=first(Family)) %>% rename(Seconds1="Seconds",area1="area")#This data frame has retained the treatment and family to be used in the analysis. 

second<-dat2C %>% group_by(Genotype) %>% arrange(Seconds) %>% summarize(Seconds = nth(Seconds,2),area=nth(area,2))  %>% rename(Seconds2="Seconds",area2="area")

third<-dat2C %>% group_by(Genotype) %>% arrange(Seconds) %>% summarize(Seconds = nth(Seconds,3),area=nth(area,3))  %>% rename(Seconds3="Seconds",area3="area")

fourth<-dat2C %>% group_by(Genotype) %>% arrange(Seconds) %>% summarize(Seconds = nth(Seconds,4),area=nth(area,4)) %>% rename(Seconds4="Seconds",area4="area")

fifth<-dat2C %>% group_by(Genotype) %>% arrange(Seconds) %>% summarize(Seconds = nth(Seconds,5),area=nth(area,5))  %>% rename(Seconds5="Seconds",area5="area")

#Merging these together into one dataframe (i love dplyr)
datLong<-first %>% left_join(second,by="Genotype") %>% left_join(third,by="Genotype")%>% left_join(fourth,by="Genotype")%>% left_join(fifth,by="Genotype")

#generating relative growth rate values (RGR)
attach(datLong)
datLong$RGR1<-(log(area2)-log(area1))/(Seconds2-Seconds1)
datLong$RGR2<-(log(area3)-log(area2))/(Seconds3-Seconds2)
datLong$RGR3<-(log(area4)-log(area3))/(Seconds4-Seconds3)
datLong$RGR4<-(log(area5)-log(area4))/(Seconds5-Seconds4)
detach(datLong)

datLong


```





#Influence of gluc and flav on defence. 
```{r}

hist(dat$Fern)
hist(dat$ThripsDam)
hist(dat$WhiteFungDam)
hist(dat$BlackPathDam)

dat$BlackPathDam

#This data is very zero inflated and a poisson model will be too overdispersed. A Negative binomial distribution will be used. 
```



#Modelling : Negative binomial on Thrips Damage. 
```{r}
library("glmmTMB")

dat$ThripsDam<-ceiling(dat$ThripsDam)

fit_full<-glmmTMB(ThripsDam~treatment+gluc_Conc+flav_Conc+(1|Family/Tag)+(1|gh_bench),family=nbinom2,data=dat)

fit_1<-update(fit_full,~.-treatment)
anova(fit_full,fit_1) #treatment is unimportant. 

fit_2<-update(fit_1,~.+ gluc_Conc:flav_Conc)

anova(fit_1,fit_2)#interaction unimportant. 


fit_1<-glmmTMB(ThripsDam~gluc_Conc+flav_Conc+(1|Family/Tag)+(1|gh_bench),family=nbinom2,data=dat[!is.na(dat$flav_Conc),])

fit_2<-update(fit_1,~.- flav_Conc)
anova(fit_1,fit_2)#
#Flav_Conc is highly important. 



fit_2<-update(fit_1,~.- gluc_Conc)
anova(fit_1,fit_2)# gluc conc does not appear to be significant. I will try using the whole gluc conc data set. 

fit_1<-glmmTMB(ThripsDam~gluc_Conc+(1|Family/Tag)+(1|gh_bench),family=nbinom2,data=dat[!is.na(dat$gluc_Conc),])

fit_2<-update(fit_1,~.- gluc_Conc)
anova(fit_1,fit_2)# gluc conc is highly significant however when flavonoids are excluded. Therefore secondary compounds, both glucosinolates and flavonoids predicts reduction in thrips damage. 

#Best model includes both, however, putting both in the same model will distribute the effect amoungst both of the terms. This might not be appropriate, but i will check for this. 

fit_g<-glmmTMB(ThripsDam~gluc_Conc+(1|Family/Tag)+(1|gh_bench),family=nbinom2,data=dat[!is.na(dat$gluc_Conc),])
fit_f<-glmmTMB(ThripsDam~flav_Conc+(1|Family/Tag)+(1|gh_bench),family=nbinom2,data=dat[!is.na(dat$flav_Conc),])
fit_full<-glmmTMB(ThripsDam~gluc_Conc+flav_Conc+(1|Family/Tag)+(1|gh_bench),family=nbinom2,data=dat)

summary(fit_f)$coef
summary(fit_g)$coef
summary(fit_full)$coef

#As expected, the coefficients are distributed evenly amoungst gluc and flav concentrations, the conclusion is that both reduce thrips abundance and it is not possible to tell which does the most, however, the effect is more significant for flavonoids. 

plot(resid(fit_f))
plot(resid(fit_f)^2)
#The model appears to be a very good fit. 


#Permutation test....
fit_f<-glmmTMB(ThripsDam~flav_Conc+(1|Family/Tag)+(1|gh_bench),family=nbinom2,data=dat)
summary(fit_f)

datPerTest<-dat
zStore<-c()
for(i in 1:500){
  #Randomize flavonoid concentration.
  datPerTest$flav_Conc<-sample(dat$flav_Conc,length(dat$flav_Conc),replace = F)
  
  #Run Model with randomized flavonoid concentration and extract test statistic
  flavZval<-summary(update(fit_f,data=datPerTest))$coef[[1]][2,3]
  
  #Store z value.
  zStore[i]<-flavZval
}

sum(zStore<=-5.356)/length(zStore)
#Wth the permutation test, the p value is still 0, which indicates that this model is a very good fit for the data. 
hist(zStore)

```


#Visualize: Thrips Damage predicted by glucosinolates and flavonoids

```{r}

source("GGPlot_Themes.R")
#Reversing link function, to estimate the data on the response scale. 
flavSlope=function(x){
  y=exp(-2.4621 *x+2.9746)
  return(y)
}
glucSlope=function(x){
  y=exp(-2.7869*x+3.5755)
  return(y)
}

#Determining x range to fit the line to 
flavplot=seq(min(dat$flav_Conc,na.rm = T),max(dat$flav_Conc,na.rm=T),length.out = 688)
glucplot=seq(min(dat$gluc_Conc,na.rm = T),max(dat$gluc_Conc,na.rm=T),length.out = 708)

#Calculating slope values
flavy<-flavSlope(flavplot)
glucy<-glucSlope(glucplot)


#tiff("Defence_Figures/FlavonoidThrips.tiff", units="in", width=10, height=6, res=300)
ggplot(dat[!is.na(dat$flav_Conc),])+
  geom_point(aes(y=ThripsDam,x=flav_Conc))+theme_simple()+
  geom_path(x=flavplot,y=flavy,size=1,colour="#999999")+
  
  scale_y_continuous(breaks=seq(0,70,5))+
  ylab("Thrips Damage")+
  xlab(bquote(bold("[Total Flavonoid] " (mg/ml))))
#dev.off()

ggplot(dat[!is.na(dat$gluc_Conc),])+
  geom_point(aes(y=ThripsDam,x=gluc_Conc))+theme_simple()+
  geom_path(x=glucplot,y=glucy,size=1,colour="#999999")+
  
  scale_y_continuous(breaks=seq(0,70,5))+
  ylab("Thrips Damage")+
  xlab(bquote(bold("[Total Glucosinolate] " (mg/ml))))

```






I think that a logistic regression is more appropriate for fungal abundance, because the count of fungal infection could be arbuitrary, especially when fungal patches were large or the leaf was completely covered. 

#WhitePathDam -- logistic regression 
```{r}

#This is the biggest model that could converge. ... interaction with flavonoid could not.
fit_full_g<-glmer(WhiteFungLogis~treatment*gluc_Conc+flav_Conc+(1|Family),family=binomial,data=dat2[!is.na(dat2$flav_Conc),])

fit.1<-update(fit_full_g,~.-flav_Conc)
anova(fit_full_g,fit.1) #Flavonoid Concentration is not significant. 


#Testing glucosinolate treatment interaction. 
fit.2<-glmer(WhiteFungLogis~treatment*gluc_Conc+(1|Family),family=binomial,data=dat2)
fit.3<-glmer(WhiteFungLogis~treatment+gluc_Conc+(1|Family),family=binomial,data=dat2)
anova(fit.2,fit.3) #Glucosinolate:Treatment interaction is not significant.

#Testing glucosinolate involvment at all.
fit.4<-glmer(WhiteFungLogis~treatment+(1|Family),family=binomial,data=dat2[!is.na(dat2$gluc_Conc),])
fit.3<-glmer(WhiteFungLogis~treatment+gluc_Conc+(1|Family),family=binomial,data=dat2[!is.na(dat2$gluc_Conc),])
anova(fit.4,fit.3) #Glucosinolate Concentration is not a significant predictor

#Testing effect of treatment
fit.4<-glmer(WhiteFungLogis~treatment+(1|Family),family=binomial,data=dat2)
fit.5<-glmer(WhiteFungLogis~1+(1|Family),family=binomial,data=dat2)
anova(fit.4,fit.5)
#treatment is significant......

summary(fit.4) #Garlic mustard in the maple treatment have less occurence of fungal colonization. 

plot(fit.4)

#Permutation test
datPerTest<-dat2
treatzStore<-c()
for(i in 1:2000){
  #Randomize flavonoid concentration.
  datPerTest$treatment<-sample(dat2$treatment,length(dat2$treatment),replace = F)
  #Run Model with randomized flavonoid concentration and extract test statistic
  treatZval<-summary(update(fit.4,data=datPerTest))$coef[3,3]
  
  #Store z value.
  treatzStore[i]<-treatZval
}
sum(treatzStore<=-2.916)/length(treatzStore)
#Simulated p value is about 0.0025, which is very similar to the 0.003 I observed, suggesting this is a good model. 
```



#Visualizing -- effect of treatment on proportion of fungal abundance. 
```{r}
#Summarizing for Display: generating frequency of white funal infection by treatment
plot<-dat2 %>% drop_na(WhiteFungLogis) %>% group_by(treatment) %>% summarize(PercWhitFung=sum(WhiteFungLogis)/length(WhiteFungLogis)*100)


#tiff("Defence_Figures/TreatMeanWhiteFung.tiff", units="in", width=8, height=5, res=300)
ggplot(plot)+
  geom_col(aes(x=treatment,y=PercWhitFung,fill=treatment))+theme_simple()+ylab("Fungal Abundance\n (% Infected)")+
  scale_y_continuous(breaks = seq(0,30,5))+
  scale_x_discrete(name="",labels=c("Alone","Garlic Mustard","Maple"))+
  scale_fill_manual(values=c("#009E73","#56B4E9","#E69F00"),labels=c("Alone","Garlic Mustard","Maple"))+
  theme_simple_multiCol()+theme(axis.title.y =  element_text(color = "black", size = 16, face = "bold",margin=margin(3,20,3,0)))
#dev.off()

```



#Modelling: Negative Binomial -- BlackPathDam 
```{r}
dat$BlackPathDam<-ceiling(dat$BlackPathDam)
fit_full<-glmmTMB(BlackPathDam~treatment+gluc_Conc+flav_Conc+(1|Family/Tag)+(1|gh_bench),family=nbinom2,data=dat)
fit_full_0.1<-glmmTMB(BlackPathDam~treatment+gluc_Conc+flav_Conc+(1|Family/Tag),family=nbinom2,data=dat)

anova(fit_full,fit_full_0.1)#gh_Bench was not an important random effect in this model. 

fit_full_0.1<-glmmTMB(BlackPathDam~treatment+gluc_Conc+flav_Conc+(1|Family/Tag),family=nbinom2,data=dat[!is.na(dat$flav_Conc),])
fit_1<-glmmTMB(BlackPathDam~treatment+gluc_Conc+(1|Family/Tag),family=nbinom2,data=dat[!is.na(dat$flav_Conc),])
anova(fit_full_0.1,fit_1) #Flav conc is a very important predictor of black pathogen damage. 

fit_1<-glmmTMB(BlackPathDam~treatment+gluc_Conc+(1|Family/Tag),family=nbinom2,data=dat[!is.na(dat$gluc_Conc),])
fit_2<-glmmTMB(BlackPathDam~treatment+(1|Family/Tag),family=nbinom2,data=dat[!is.na(dat$gluc_Conc),])
anova(fit_1,fit_2) #Glucosinolates are not a significant predictor at all. 


fit_2<-glmmTMB(BlackPathDam~treatment+flav_Conc+(1|Family/Tag),family=nbinom2,data=dat)
fit_3<-glmmTMB(BlackPathDam~flav_Conc+(1|Family/Tag),family=nbinom2,data=dat)
anova(fit_3,fit_2) #Treatment is a significant predictor of black pathogen damage. 


#Therefore, the best model is one with flavonoids and treatment. 

summary(fit_2)

plot(resid(fit_2)) #Model fits well.

#Permutation test. 

datPerTest<-dat
zStoreflav<-c()
zStoretreat<-c()
for(i in 1:500){
  #Randomize flavonoid concentration.
  datPerTest$flav_Conc<-sample(dat$flav_Conc,length(dat$flav_Conc),replace = F)
  datPerTest$treatment<-sample(dat$treatment,length(dat$treatment),replace = F)
  
  #New Model with randomized treatment and flavonoids
  newMod<-update(fit_2,data=datPerTest)
  
  # Extract test statistic
  flavZval<-summary(newMod)$coef[[1]][4,3]
  treatZval<-summary(newMod)$coef[[1]][3,3]
  
  #Store z value.
  zStoreflav[i]<-flavZval
  zStoretreat[i]<-treatZval
}

#treatment p value 
sum(zStoretreat<=-2.082)/length(zStoretreat) 
#The p value for treatment is 0.022, which is very close to the observed 0.03 p value.
sum(zStoreflav<=-2.8863)/length(zStoreflav)
#The p value of flavonoids is 0.002, which is very close to the observed 0.003 p value.
#Conclusion: This is a good model. 
```


#Visualizing: The effect of treatment and flavonoid abundance on black pathogen abundance. 
```{r}
source("GGPlot_Themes.R")
#Reversing link function, to estimate the data on the response scale. 
flavSlope=function(x,int){
  y=exp(-1.5060*x+int)
  return(y)
}



#Determining x range to fit the line to 
flavplot=seq(min(dat$flav_Conc,na.rm = T),max(dat$flav_Conc,na.rm=T),length.out = 680)


#Calculating slope values
flavyA<-flavSlope(flavplot, 0.7826)
flavyGM<-flavSlope(flavplot, 0.7826+0.1049)
flavyM<-flavSlope(flavplot, 0.7826-0.5204 )


#tiff("Defence_Figures/FlavonoidThrips.tiff", units="in", width=10, height=6, res=300)
ggplot(dat[!is.na(dat$flav_Conc),])+
  geom_point(aes(y=BlackPathDam,x=flav_Conc,colour=treatment))+
  theme_simple()+
  geom_path(x=flavplot,y=flavyA,size=1,colour="#009E73")+
  geom_path(x=flavplot,y=flavyGM,size=1,colour="#56B4E9")+
  geom_path(x=flavplot,y=flavyM,size=1,colour="#E69F00")+
  
  scale_colour_manual(values=c("#009E73","#56B4E9","#E69F00"),labels=c("Alone","Garlic Mustard","Maple"))+
  
  scale_y_continuous(breaks=seq(0,70,5))+
  ylab("Black Pathogen Damage")+
  xlab(bquote(bold("[Total Flavonoid] " (mg/ml))))
```


#Visualization -- Black pathogen damage by treatment
```{r}

plot2<-dat2 %>% drop_na(BlackPathDam) %>% group_by(treatment) %>% summarize(BlackPathAve=mean(BlackPathDam,na.rm=T))

ggplot(plot2)+
  geom_col(aes(x=treatment,y=BlackPathAve,fill=treatment))+theme_simple()+ylab("Black Pathogen Infection\n(spots/leaf)")+xlab("Treatment")+theme(legend.position = "none")

```

#Visualizing genetic variation and greenhouse variation, which will be controlled for. 
```{r}
#GH Bench
ggplot(dat2)+
  geom_point(aes(y=gluc_Conc,x=gh_bench,colour=as.factor(gh_bench)))

#GH Col
ggplot(dat2)+
  geom_point(aes(y=gluc_Conc,x=gh_col,colour=as.factor(gh_bench)))

#Investigating genetic differences by treatment
#gluc_Conc
boxplot(gluc_Conc~Family,data=dat2[dat2$treatment=="a",])
boxplot(gluc_Conc~Family,data=dat2[dat2$treatment=="m",])
boxplot(gluc_Conc~Family,data=dat2[dat2$treatment=="gm",])
#bodymass
boxplot(GM_TotalLeaf_Area~Family,data=dat2[dat2$treatment=="a",])
boxplot(GM_TotalLeaf_Area~Family,data=dat2[dat2$treatment=="m",])
boxplot(GM_TotalLeaf_Area~Family,data=dat2[dat2$treatment=="gm",])
```



#Visualization-- The Detriment of pathogens and ferns
```{r}
summary(lmer(GM_TotalLeaf_Area ~ Fern+
               (1 | Family) + (1 | gh_bench/gh_col) ,data=dat2))

a<-ggplot(dat2)+
  geom_point(aes(y=GM_TotalLeaf_Area,x=BlackPathDam))+theme_simple_multiCol()+
  geom_abline(intercept=8281.76,slope = -105.28,size=1.5)+
  xlab("Black Pathogen Damage")

b<-ggplot(dat2[dat2$WhiteFungDam<30,])+
  geom_point(aes(y=GM_TotalLeaf_Area,x=WhiteFungDam),colour="#999999")+theme_simple_multiCol()+
  theme(axis.title.x = element_text(color = "#999999", size = 16, face = "bold",margin=margin(3,0,3,0)),
  )+xlab("Powdery Mildew Damage")

c<-ggplot(dat2)+
  geom_point(aes(y=GM_TotalLeaf_Area,x=ThripsDam),colour="#E69F00")+theme_simple_multiCol()+xlab("Thrips Damage")+
  theme(axis.title.x = element_text(color = "#E69F00", size = 16, face = "bold",margin=margin(3,0,3,0)))

d<-ggplot(dat2)+
  geom_point(aes(y=GM_TotalLeaf_Area,x=Fern),colour="#009E73")+theme_simple_multiCol()+xlab("Fern Abundance")+
  geom_abline(intercept=8003.44,slope = -179.47,size=1.5,color="#009E73")+
  theme(axis.title.x = element_text(color = "#009E73", size = 16, face = "bold",margin=margin(3,0,3,0)))



plot<-plot_grid(a, b,ncol=2,rel_widths = c(1,1))

plot2<-plot_grid(d, c,ncol=2,rel_widths = c(1,1))

plot3<-plot_grid(a,b,c,d,ncol=2,rel_widths = c(1,1))

plot4<-plot_grid(a,b,c,ncol=1,rel_widths = c(1,1,1))
plot5<-plot_grid(a,b,c,ncol=3,rel_widths = c(1,1,1))


y.grob<-textGrob(bquote(bold("Shoot Area "(mm^2))),gp=gpar(fontface="bold",fontsize=20),rot=90)


#tiff("Selection_Figures/PathogenEffect.tiff", units="in", width=14, height=6, res=300)
grid.arrange(plot,left=y.grob)
#dev.off()

#tiff("Selection_Figures/PathogenEffect2.tiff", units="in", width=14, height=6, res=300)
grid.arrange(plot2,left=y.grob)
#dev.off()

#tiff("Selection_Figures/PathogenEffect3.tiff", units="in", width=14, height=10, res=300)
grid.arrange(plot3,left=y.grob)
#dev.off

grid.arrange(plot4,left=y.grob)


#tiff("Selection_Figures/PathogenEffect5.tiff", units="in", width=16, height=5, res=300)
grid.arrange(plot5,left=y.grob)
#dev.off

#tiff("Selection_Figures/FernEffect.tiff", units="in", width=8, height=6, res=300)
grid.arrange(d,left=y.grob)
#dev.off

```

```{r}

#Percent variance explained: 
337263/ ( 395515 + 5196413+337263 )
580.7/ ( 395515 + 5196413+337263 )


#Does fern affect performance without accounting for leaf area (it is not affected by leaf area)
fit7<-glmmTMB(GM_TotalLeaf_Area~treatment+Fern+(1|Family)+(1|gh_bench), data=dat2)
summary(fit7) #Without account for leaf area, fern is a very significant predictor of performance and is correlated with reduced garlic mustard size. Whether it is doing the reduction or it can simply appear when garlic mustard is a worst performer is unclear. These estimates should be used to visualize the effect of treatment on performance. 

#To what degree does family influence performance?

fitfull<-glmmTMB(GM_TotalLeaf_Area~treatment+Fern+(1|Family)+(1|gh_bench), data=dat2)
fitfull2<-glmmTMB(GM_TotalLeaf_Area~treatment+Fern+(1|gh_bench), data=dat2)

anova(fitfull,fitfull2) #Family only predicts performance when pathogens are accounted for. 
```


#What influences performance? 
```{r}
library(glmmTMB)

#Modelling random effects.
fitfull<-glmmTMB(GM_TotalLeaf_Area~treatment+BlackPathDam+WhiteFungLogis+ThripsDam+Fern+GM_Leaf_Area+(1|Family)+(1|gh_bench), data=dat2)

fitfull1<-glmmTMB(GM_TotalLeaf_Area~treatment+BlackPathDam+WhiteFungLogis+ThripsDam+Fern+GM_Leaf_Area+(1|Family), data=dat2)

fitfull2<-glmmTMB(GM_TotalLeaf_Area~treatment+BlackPathDam+WhiteFungLogis+ThripsDam+Fern+GM_Leaf_Area+(1|gh_bench), data=dat2)

fitfull3<-glmmTMB(GM_TotalLeaf_Area~treatment+BlackPathDam+WhiteFungLogis+ThripsDam+Fern+GM_Leaf_Area+(1|gh_bench)+(1|Family:treatment), data=dat2)

fitfull4<-glmmTMB(GM_TotalLeaf_Area~treatment+BlackPathDam+WhiteFungLogis+ThripsDam+Fern+GM_Leaf_Area+(1|gh_bench)+(1|Family/treatment), data=dat2)

#Is family important? 
anova(fitfull,fitfull2) #Yes, family predicts performance. 
#Is gh bench important? 
anova(fitfull,fitfull1) #Yes, gh_bench predicts performance. 
#Is there a GxE interaction? 
anova(fitfull,fitfull3) 
anova(fitfull,fitfull4) #No there does not appear to be. 

#Modelling fixed effects. 
fit<-glmmTMB(GM_TotalLeaf_Area~treatment+BlackPathDam*WhiteFungLogis*ThripsDam+Fern+GM_Leaf_Area+(1|Family)+(1|gh_bench), data=dat2)

fit2<-update(fit,~.-ThripsDam:BlackPathDam:WhiteFungLogis)
anova(fit,fit2) #Not a significant three way interaction.

fit3<-update(fit2,~.-BlackPathDam:WhiteFungLogis)
anova(fit3,fit2) #There is not a significant two way interaction, 

fit4<-update(fit3,~.-BlackPathDam:ThripsDam)
anova(fit4,fit3)
#There is not a significant interaction between black path dam and thrips dam. 

fit5<-update(fit4,~.-WhiteFungLogis:ThripsDam)
anova(fit5,fit4) #There is not a significant whitefung dam and thrips dam interaction. 

summary(fit5)

#Fern is not significant. 
fit6<-glmmTMB(GM_TotalLeaf_Area~treatment+BlackPathDam+WhiteFungLogis+ThripsDam+GM_Leaf_Area+(1|Family)+(1|gh_bench), data=dat2)

fit6<-glmmTMB(GM_TotalLeaf_Area~treatment+Standardize(BlackPathDam)+Standardize(WhiteFungLogis)+Standardize(ThripsDam)+Standardize(GM_Leaf_Area)+(1|gh_bench)+(1|Family), data=dat2[!is.na(dat2$ThripsDam),])

#Is black path dam significant
fit7<-glmmTMB(GM_TotalLeaf_Area~treatment+Standardize(WhiteFungLogis)+Standardize(ThripsDam)+Standardize(GM_Leaf_Area)+(1|gh_bench)+(1|Family), data=dat2[!is.na(dat2$BlackPathDam),])
anova(fit6,fit7)

#Is  White Path dam significant
fit8<-glmmTMB(GM_TotalLeaf_Area~treatment+Standardize(BlackPathDam)+Standardize(ThripsDam)+Standardize(GM_Leaf_Area)+(1|gh_bench)+(1|Family), data=dat2[!is.na(dat2$WhiteFungLogis),])
anova(fit6,fit8)

#Is thrips  dam significant
fit9<-glmmTMB(GM_TotalLeaf_Area~treatment+Standardize(BlackPathDam)+Standardize(WhiteFungLogis)+Standardize(GM_Leaf_Area)+(1|gh_bench)+(1|Family), data=dat2[!is.na(dat2$ThripsDam),])
anova(fit6,fit9)

summary(fit6) #What happens when i remove the effect of leaf size?
fit7<-glmmTMB(GM_TotalLeaf_Area~treatment+Standardize(BlackPathDam)+Standardize(WhiteFungLogis)+Standardize(ThripsDam)+(1|gh_bench)+(1|Family), data=dat2[!is.na(dat2$ThripsDam),])
summary(fit7) #White fung dam and thrips dam lose significance



#Percent variance explained: 
337263/ ( 395515 + 5196413+337263 )
580.7/ ( 395515 + 5196413+337263 )


#Does fern affect performance without accounting for leaf area (it is not affected by leaf area)
fit7<-glmmTMB(GM_TotalLeaf_Area~treatment+Fern+(1|Family)+(1|gh_bench), data=dat2)
summary(fit7) #Without account for leaf area, fern is a very significant predictor of performance and is correlated with reduced garlic mustard size. Whether it is doing the reduction or it can simply appear when garlic mustard is a worst performer is unclear. These estimates should be used to visualize the effect of treatment on performance. 

#To what degree does family influence performance?

fitfull<-glmmTMB(GM_TotalLeaf_Area~treatment+Fern+(1|Family)+(1|gh_bench), data=dat2)
fitfull2<-glmmTMB(GM_TotalLeaf_Area~treatment+Fern+(1|gh_bench), data=dat2)

anova(fitfull,fitfull2) #Family only predicts performance when pathogens are accounted for. 
```





#Is there a cost/benefit to glucosinolate production? 
```{r}
fit6<-glmmTMB(GM_TotalLeaf_Area~treatment+BlackPathDam+WhiteFungLogis+ThripsDam+GM_Leaf_Area+gluc_Conc+(1|Family)+(1|gh_bench), data=dat2)
summary(fit6) #No, infact glucosinolates are positively correlated with performance. 


#Is the cost dependent on the treatment
fit7<-glmmTMB(GM_TotalLeaf_Area~treatment*gluc_Conc+BlackPathDam+WhiteFungLogis+ThripsDam+GM_Leaf_Area+(1|Family)+(1|gh_bench), data=dat2)
summary(fit7) 
anova(fit6,fit7)#There is an interaction between treatment. However it is postive in all treatments and maple performance was not accounted for so this model is nullified as it is uninformative. 

#Looking at standardized coefficients 
fit<-glmmTMB(GM_TotalLeaf_Area~treatment+Standardize(BlackPathDam)+WhiteFungLogis+Standardize(ThripsDam)+GM_Leaf_Area+(1|Family)+(1|gh_bench), data=dat2)
summary(fit)

fit<-glmmTMB(GM_TotalLeaf_Area~treatment+Standardize(Fern)+(1|Family)+(1|gh_bench), data=dat2)
summary(fit)


fit6<-glmmTMB(GM_TotalLeaf_Area~treatment+BlackPathDam+WhiteFungLogis+ThripsDam+GM_Leaf_Area+flav_Conc+Fern+(1|Family)+(1|gh_bench), data=dat2)
summary(fit6) #As are flavonoids. 

#Is the benefit/cost dependent on the treatment? 
fit7<-glmmTMB(GM_TotalLeaf_Area~treatment*flav_Conc+BlackPathDam+WhiteFungLogis+ThripsDam+GM_Leaf_Area+Fern+(1|Family)+(1|gh_bench), data=dat2)
summary(fit6) 

anova(fit6,fit7) #There is not a significant interaction with flavonoids. 

#Is the benefit/cost dependent on glucosinolates not being in the model?

fit7<-glmmTMB(GM_TotalLeaf_Area~treatment*gluc_Conc+flav_Conc+BlackPathDam+WhiteFungLogis+ThripsDam+GM_Leaf_Area+Fern+(1|Family)+(1|gh_bench), data=dat2)
summary(fit7) 

fit8<-glmmTMB(GM_TotalLeaf_Area~treatment+gluc_Conc+flav_Conc+BlackPathDam+WhiteFungLogis+ThripsDam+GM_Leaf_Area+Fern+(1|Family)+(1|gh_bench), data=dat2)
summary(fit7) 
anova(fit7,fit8)

```



#Permutation test
datPerTest<-datFern
zStoretreatGM<-c()
zStoretreatM<-c()

for(i in 1:500){
  #Randomize flavonoid concentration.
  datPerTest$treatment<-sample(datPerTest$treatment,length(datPerTest$treatment),replace = F)
  
  #New Model with randomized treatment and flavonoids
  newMod<-update(fit_1,data=datPerTest)
  
  # Extract test statistics for maple (M) and garlic mustard (gm) treatments
  MtreatZval<-summary(newMod)$coef[[1]][3,3]
  GMtreatZval<-summary(newMod)$coef[[1]][3,2]
  
  #Store z value.
  zStoretreatGM[i]<-GMtreatZval
  zStoretreatM[i]<-MtreatZval
}

sum(zStoretreatM>=2.316180)/length(zStoretreatM)
#Estimated p value for the maple treatment is 0.01. This is very close to the actual p value of 0.02
sum(zStoretreatGM>=2.231574)/length(zStoretreatGM)
#Estimated p value for the garlic mustard treatment is 0. This is close to the actual p value of 0.02

summary(fit_1)



#Visualizing --- the distribution of pathogens by treatment. 
```{r}


plot4<-dat2 %>% drop_na(Fern) %>% group_by(treatment) %>% summarize(Fern=mean(Fern,na.rm=T))


table(dat$comp_number)
ggplot(plot4)+
  geom_col(aes(x=treatment,y=Fern,fill=treatment))+theme_simple()+ylab("Average Fern Abundance\n(ferns/pot)")+xlab("Treatment")+theme(legend.position = "none")+
  scale_x_discrete(name="",labels=c("Alone","Garlic Mustard","Maple"))+
  scale_fill_manual(values=c("#009E73","#56B4E9","#E69F00"),labels=c("Alone","Garlic Mustard","Maple"))+
  theme_simple_multiCol()+theme(axis.title.y =  element_text(color = "black", size = 16, face = "bold",margin=margin(3,20,3,0)))
#dev.off()

```



# Modelling: Effect of GM on Maple Performance
```{r}
#The influence of family, glucosinolates and flavonoids on competitive ability
hist(Maple$Maple_TotalLeafArea_End,breaks=50)
#There is a limit at zero, so this is not a true normal distribution. It resembles a reverse gaussian distribution except this distribution has values at zero, infact it is slighly zer inflated for the distibution. Therefore i think a tweedie distibution is best. 

#Modelling Random effects 

fitTwFull<-glmmTMB(Maple_TotalLeafArea_End~PC1Tot+gluc_Conc+MapleDamage+(1|Family)+(1|gh_bench),family = tweedie, data=MapleModel)

fitTwFull0.1<-glmmTMB(Maple_TotalLeafArea_End~PC1Tot+gluc_Conc+MapleDamage+(1|gh_bench),family = tweedie, data=MapleModel)

fitTwFull0.2<-glmmTMB(Maple_TotalLeafArea_End~PC1Tot+gluc_Conc+MapleDamage+(1|Family),family = tweedie, data=MapleModel)

fitTwFull0.3<-glmmTMB(Maple_TotalLeafArea_End~PC1Tot+gluc_Conc+MapleDamage,family = tweedie, data=MapleModel)

#Is there an effect of family? 
anova(fitTwFull,fitTwFull0.1) #No

#Is there an effect of bench? 
anova(fitTwFull,fitTwFull0.2) #No

#It is not important to include random effects. 


#Modelling: Fixed effects?

fitTw1<-glmmTMB(Maple_TotalLeafArea_End~PC1Tot*gluc_Conc+MapleDamage,family = tweedie, data=MapleModel)
summary(fitTw1) #The interaction is not important, but glucosinolate concentration and maple damage are. 

fitTw2<-glmmTMB(Maple_TotalLeafArea_End~PC1Tot+gluc_Conc+MapleDamage,family = tweedie, data=MapleModel[!is.na(MapleModel$gluc_Conc),])
summary(fitTw2)
#Glucosinolate concentration and MapleDamage are both highly significant.


#The question is, however, do glucosinolates significant reduce maple performance because they are inhibiting maple (allelopathy) OR because those with high glucosinolate expression are better competitors and aquire more resources to then allocate more resources to defenes? (Chicken or the egg)

#Is the negative effect of glucosinolates on maple performance removed once Chlorophyll A expression is added to the model? This would control for garlic mustard performance because those that produce more chlorophyll can also produce more glucosinolates. 
fitTw3<-glmmTMB(Maple_TotalLeafArea_End~PC1Tot+MapleDamage+gluc_Conc+ChlorA,family = tweedie, data=MapleModel)
summary(fitTw3)
# Including Chlorophyll A completely removed the effect of glucosinolates. This suggests that glucosinolates may be correlated with maple performance due to plants that reduce maple performance having higher fitness and more resources for glucosinolates without glucosinolates actually being responsible. On the other hand, those with more glucosinolates may have higher fitness and be able to produce more Chlorophyll A because they are inhibiting maple.


#Checking Model Assumptions. 
summary(fitTw2)
plot(resid(fitTw2))
hist(resid(fitTw2),breaks=20)#Looks reasonably normal to me.

fitV<-lm(Maple_TotalLeafArea_End~PC1Tot+MapleDamage+gluc_Conc+ChlorA,data=MapleModel)
summary(fitV)

vif(fitV) #Vif for ChlorA and gluc concentration are both below 2. This suggest that gluc_Conc and Chlor concentration are not high enough that the effect of the two are indistinguishable. This suggests that the entire effect of glucosinolates is due to chlorophyll concentration. 


#Permutation test to check the reliability of the p values in the best model.
fitTw3<-glmmTMB(Maple_TotalLeafArea_End~PC1Tot+gluc_Conc+MapleDamage,family = tweedie, data=MapleModel)
summary(fitTw3)
# 
# summary(fitTw3)$coef
# datPerTest<-MapleModel
# zStore<-c()
# # #for(i in 1:500){
#   #Randomize flavonoid concentration.
#   datPerTest$gluc_Conc<-sample(datPerTest$gluc_Conc,length(datPerTest$gluc_Conc),replace = F)
#   
#   #Run Model with randomized flavonoid concentration and extract test statistic
#   glucZval<-summary(update(fitTw3,data=datPerTest))$coef[[1]][3,3]
#   
#   #Store z value.
#   zStore[i]<-glucZval
#}

#sum(zStore<=-2.883)/length(zStore)

#The simulated p value was 0.002, which is close to the observed value of 0.004. This is a reliable model. 

```




#Is there evidence that those in the maple and garlic mustard treatment have the same relationship between glucosinolate and chlorA? If glucosinolates are truly beneficial in one treatment but not in another, we might expect for the relationship between glucosinolates and ChlorA to be stronger in the treatment in which glucosinolates are beneficial. This is because if those that create more chlorophyll a passively allocate more resources to glucosinolates, this relationship may be a weaker relationship than if glucosinolates are increasing fitness and allowing an increased allocation to chlrophyll A. 
```{r}
fit<-lm(gluc_Conc~ChlorA*treatment,data=dat)
fit2<-lm(gluc_Conc~ChlorA+treatment,data=dat)
anova(fit,fit2)
#There is not a significant interaction when logged data are used...... This figure is now obsolete. 
summary(fit)
summary(fit2)


datplot<-dat %>% filter(treatment!="mcnt")
datplot$treatment<-droplevels.factor(datplot$treatment)

source("GGPlot_Themes.R")
tiff("Other_Figures/Gluc_Chlor*trement_Cor.tiff", units="in", width=10, height=6, res=300)
ggplot(datplot)+
  geom_point(aes(x=ChlorA,y=gluc_Conc,colour=treatment))+
  geom_abline(slope = 0.37585,intercept = 0.85057,colour="#009E73",size=1)+#alone
  geom_abline(slope = 0.37585+0.20463,intercept = 0.85057-0.09106,colour="#56B4E9",size=1)+#gm
  geom_abline(slope = 0.37585+0.18856,intercept = 0.85057-0.09837,colour="#E69F00",size=1)+#maple
  scale_colour_manual(values=c("#009E73","#56B4E9","#E69F00","black"),labels=c("Alone","Garlic Mustard","Maple"))+
  ylab(bquote(bold("[Glucosinolate] "(mg/ml))))+
  xlab(bquote(bold("[Chlorophyll A] "(μg/ml))))+
  theme_simple()
dev.off()

#This demonstrates that the relationship between glucosinolate and chlorophyll is steeper in the treatments with competition. Therefore, there may be a benefit to producing glucosinolates in the competition treatments? Those with high glucosinolates in the maple treatment have higher chlorophyll levels than those with high glucosinolates in the alone treatment. In other words, those in the competition treatment have lower levels of chlorophyll at low levels of glucosinolates, and have higher level of chlorphyll at higher levels of glucosinolates. In other words (again) glucosinolates may be increasing performance in the competition treatments. 
```






