
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

theme_simple<-function(){
  #base theme
  theme_classic() %+replace%
    #theme additions
    theme(
      
      #Legend
      legend.position="top",
      legend.title = element_blank(),
      legend.text = element_text(color = "black", size = 14, face = "bold"),
      legend.background=element_rect(size=0.5,colour="black"),
      #Theme of axis text and ticks
      axis.text.x = element_text(color = "black", size = 14, face = "bold"),
      axis.text.y = element_text(color = "black", size = 14, face = "bold"),
      axis.ticks = element_blank(),
      
      #Theme of Axis titles
      axis.title.x =  element_text(color = "black", size = 18, face = "bold",margin=margin(10,0,0,0)),
      axis.title.y =  element_text(color = "black", size = 18, face = "bold",angle=0,vjust= 0.5,margin=margin(0,10,0,0)),
      
      #Plot Margin
      plot.margin = unit(c(1,1,1,1), "cm"),
      
      #Theme of Axis lines
      axis.line=element_line(colour = "black", size = 1)
    )
}

theme_simple_Legend<-function(){
  #base theme
  theme_classic() %+replace%
    #theme additions
    theme(
      
      #Legend
      legend.position="right",
      legend.title = element_blank(),
      legend.text = element_text(color = "black", size = 14, face = "bold"),
      legend.background=element_rect(size=0.5,colour="black"),
      #Theme of axis text and ticks
      axis.text.x = element_text(color = "black", size = 14, face = "bold"),
      axis.text.y = element_text(color = "black", size = 14, face = "bold"),
      axis.ticks = element_blank(),
      
      #Theme of Axis titles
      axis.title.x =  element_text(color = "black", size = 18, face = "bold",margin=margin(10,0,0,0)),
      axis.title.y =  element_text(color = "black", size = 18, face = "bold",angle=0,vjust= 0.5,margin=margin(0,10,0,0)),
      
      #Plot Margin
      plot.margin = unit(c(1,1,1,1), "cm"),
      
      #Theme of Axis lines
      axis.line=element_line(colour = "black", size = 1)
    )
}

theme_simple_vert<-function(){
  #base theme
  theme_classic() %+replace%
    #theme additions
    theme(
      
      #Legend
      legend.position="top",
      legend.title = element_blank(),
      legend.text = element_text(color = "black", size = 14, face = "bold"),
      legend.background=element_rect(size=0.5,colour="black"),
      #Theme of axis text and ticks
      axis.text.x = element_text(color = "black", size = 14, face = "bold"),
      axis.text.y = element_text(color = "black", size = 14, face = "bold"),
      axis.ticks = element_blank(),
      
      #Theme of Axis titles
      axis.title.x =  element_text(color = "black", size = 16, face = "bold",margin=margin(10,0,0,0)),
      axis.title.y =  element_text(color = "black", size = 16,angle = 90, face = "bold",vjust= 0.5,margin=margin(0,10,0,0)),
      
      #Plot Margin
      plot.margin = unit(c(1,1,1,1), "cm"),
      
      #Theme of Axis lines
      axis.line=element_line(colour = "black", size = 1)
    )
}


theme_simple_multiCol<-function(){
  #base theme
  theme_classic() %+replace%
    #theme additions
    theme(
      
      #Remove legend
      legend.position="none",
      
      #Theme of axis text and ticks
      axis.text.x = element_text(color = "black", size = 16, face = "bold"),
      axis.text.y = element_text(color = "black", size = 12, face = "bold"),
      axis.ticks = element_blank(),
      
      #Theme of Axis titles
      axis.title.x =  element_text(color = "black", size = 16, face = "bold",margin=margin(3,0,3,0)),
      
      plot.title =  element_text(color = "black", size = 16, face = "bold",margin=margin(3,0,3,0)),
      
      axis.title.y =  element_blank(),
      
      #Plot Margin
      plot.margin = unit(c(0,1,0,1), "cm"),
      
      #Theme of Axis lines
      axis.line=element_line(colour = "black", size = 1)
    )
}

theme_simple_multiCol_First<-function(){
  #base theme
  theme_classic() %+replace%
    #theme additions
    theme(
      
      #Remove legend
      legend.position="none",
      
      #Theme of axis text and ticks
      axis.text.x = element_text(color = "black", size = 16, face = "bold"),
      axis.text.y = element_text(color = "black", size = 12, face = "bold"),
      axis.ticks = element_blank(),
      
      #Theme of Axis titles
      axis.title.x =  element_text(color = "black", size = 16, face = "bold",margin=margin(3,0,3,0)),
      
      plot.title =  element_text(color = "black", size = 16, face = "bold",margin=margin(3,0,3,0)),
      
      axis.title.y =  element_text(color = "black", size = 16, face = "bold",margin=margin(0,3,0,0),angle =90),
      
      #Plot Margin
      plot.margin = unit(c(0,1,0,1), "cm"),
      
      #Theme of Axis lines
      axis.line=element_line(colour = "black", size = 1)
    )
}
