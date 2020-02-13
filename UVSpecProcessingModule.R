
#Neccessairy functions for use in main function

#This function will reformat the output from the UVspec to subtrac the control from each sample
#Importantly, the control samples are divided by two for glucosinolates and flavonoids, 
#because the concentration of the controls is twice that of the stained samples (Controls are 100ul while stains are 200ul).Using C1V1-C2V2. Dividing the absorbance by 2 will estimate the control if there were 200ul. 
ControlSubtract<-function(x,compound){
  
  if(compound=='flav'){
    plate='FlavStain'
  }
  if(compound=='gluc'){
    plate='GlucStain'
  }
  
  Cont=x[x$Compound==compound & x$Plate=="Contrl",]
  Stain=x[x$Compound==compound & x$Plate==plate,]
  
  if(any(Cont$Well!=Stain$Well)){
    print("Something Went Wrong")
    stop()
  }
  
  else{
    Stain$AdjData<-Stain$Data-Cont$Data/2
    colnames(Stain)[3]<-"StainData"
    Stain$ControlData<-Cont$Data/2
    return(Stain)
  }
}







#Flavonoid STD from Feb 7th with the final protocol
flavCONC<-function(x){
  y=2.28450507*x+0.09187857
  return(y)
}


#Glucosinolate STD from Feb 7th with the final protocol
glucCONC<-function(x){
  y=1.433029*x+0.226200
  return(y)
}

#CHLOROPHYLL

Cha<-function(A652,A665,A750){
  a652<-A652-A750
  a665<-A665-A750
  
  return( -8.0962*a652 + 16.5169*a665 )
}

Chb<-function(A652,A665,A750){
  a652<-A652-A750
  a665<-A665-A750
  return(27.4405*a652 - 12.1688*a665)
}




#Main function which will Calculate the chlorophyll concentration, and for the glucosinolate and flavonoid samples, subtract the controls for each run we have done. 
CleanData<-function(AllData){
  
  
  #Data before isolating run
  dat0<-read.csv(AllData,stringsAsFactors = F)
  
  #Prep Global data frames
  GlobalBindedGF<-data.frame()
  GlobalBindedCh<-data.frame()
  
  #Isolate data from a run (focal run based on the datafile)
  for (i in unique(dat0$DataFile)){
    
    dat<-dat0[dat0$DataFile==i,]
    
    #Remove blank values
    dat2<-dat[!is.na(dat$Sample),]
    
    
    #Removing errors generated from the UV spec machine. These are outputted simply as "#ERR", this will need to be changed. 
    dat2$Data[dat2$Data=="#ERR"]<-NA
    dat2$Data<-as.numeric(dat2$Data)
    
    
    
    #Create Subtract the control values from the glucsinolate and flavonoid value
    gluc<-ControlSubtract(dat2,'gluc')
    flav<-ControlSubtract(dat2,'flav')
    
    #Add concentration into each dataframe
    gluc$Conc<-glucCONC(gluc$AdjData)
    flav$Conc<-flavCONC(flav$AdjData)
    
    
    #STORE AND APPEND DATA FRAMES
    
    #Binding gluc and flav together
    LocalBinded<-rbind(gluc,flav)
    
    GlobalBindedGF<-rbind(GlobalBindedGF,LocalBinded)
    
    
    
    ###########
    #Chlorophyll------------------------------------------------------
    ###########
    ChlorDat<-dat2[(grepl("Chlor",dat2$Plate)),]
    
    
    Chlor750<-ChlorDat[ChlorDat$Plate=="Chlor-750",]
    Chlor652<-ChlorDat[ChlorDat$Plate=="Chlor-652",]
    Chlor665<-ChlorDat[ChlorDat$Plate=="Chlor-665",]
    
    #Check
    for(i in 1:length(Chlor750)){
      if(any(Chlor750$Well!=Chlor652$Well) | any(Chlor652$Well!=Chlor665$Well)| any(Chlor750$Well!=Chlor665$Well)
      ){
        print("Something went wrong")
        break
      }
    }
    
    #Chlorophyll concentrations
    Chlor750$ChlorA<-Cha(Chlor652$Data,Chlor665$Data,Chlor750$Data)
    Chlor750$ChlorB<-Chb(Chlor652$Data,Chlor665$Data,Chlor750$Data)
    #Adding in all of the data
    Chlor750$Data652<-Chlor652$Data
    Chlor750$Data665<-Chlor665$Data
    #Changing column name from "Data" to "Data750" as it represents this.
    colnames(Chlor750)[3]<-"Data750"
    
    #STORE DATA GLOBAL
    GlobalBindedCh<-rbind(GlobalBindedCh,Chlor750)
    
  }#End loop of focal data
  
  return(list(GlobalBindedGF,GlobalBindedCh))
}#End function


#})



#GFData<-CleanData("All_UV_Data.txt")[[1]]
#ChData<-CleanData("All_UV_Data.txt")[[2]]
