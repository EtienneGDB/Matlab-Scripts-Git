clear()
rm()
gc()

# Libraries
library(tidyverse)
library("R.matlab")
library(ggplot2)
library(dplyr)

#impVarNames<-readMat("J:\\Piano_Fatigue\\Data_Exported\\Mean_VIPSCORE_Ha.mat")
impVIP_Ha<-readMat("J:\\Piano_Fatigue\\Data_Exported\\VIPMembres_Ha.mat")
impVIP_Li<-readMat("J:\\Piano_Fatigue\\Data_Exported\\VIPMembres_Li.mat")
impVarNames<-readMat("J:\\Piano_Fatigue\\Data_Exported\\VarNames_dec_save_Ha.mat")

Membres<-data.frame('Pelvis','Trunk','Head','Shoulder','Arm','Forearm','Hand')

VarNames=matrix('Na',nrow=180,ncol=3)
for(i in 1:length(impVarNames[[1]])){
  VarNames[[i]]<-impVarNames[[1]][[i]][[1]]
}
VarNames=matrix(VarNames,nrow=180,ncol=3)
Task = 'Li'
if(Task=='Ha'){Color_bar='blue'}else{Color_bar='red'}

for(j in 1:length(impVIP_Ha[[1]])){
  VIPSCORE<-data.frame(id=c(1:180),X=VarNames,Y=impVIP_Li$VIPMembres[[j]][[1]],Task=rep('Ha',180))
  #VIPSCORE<-data.frame(id=c(1:180),X=VarNames,Y=impVIP_Li$VIPMembres[[j]][[1]],Task=rep('Li',180))
  #VIPSCORE=rbind(VIPSCORE_Ha,VIPSCORE_Li)
  #VIPSCORE=VIPSCORE[order(VIPSCORE[,1]),]
  #VIPSCORE$id<-c(1:nrow(VIPSCORE))

# ----- This section prepare a dataframe for labels ---- #
# Get the name and the y position of each label
label_data <- VIPSCORE

# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)

# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
label_data$hjust<-ifelse( angle < -90, 1, 0)

# flip angle BY to make them readable
label_data$angle<-ifelse(angle < -90, angle+180, angle)
# ----- ------------------------------------------- ---- #
# prepare a data frame for base lines
empty_bar<-0.1
base_data <- VIPSCORE %>% 
  group_by(X.1) %>% 
  summarize(start=min(id)+0, end=max(id)-0) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))
base_data <- base_data[order(base_data[,2]),]

# calculate the ANGLE of the base_data
number_of_bar <- nrow(base_data)
angle <-  0 - 360 * (c(1:15)-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)

# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
base_data$hjust<-ifelse( angle < -90, 1, 0)

# flip angle BY to make them readable
base_data$angle<-ifelse(angle < -90, angle+180, angle)
base_data$angle<-ifelse(base_data$angle < -90, base_data$angle+180, base_data$angle)
#---
s<-' '
O<-data.frame()
for(i in 1:nrow(base_data)){
  temp <- label_data[which(label_data[,2]==as.character(base_data[i,1])),]
  temp[,2]<-matrix(rep(s,nchar(as.character(base_data[i,1]))),nrow=nrow(temp),ncol=1)
  temp[nrow(temp)/2,2]<-as.character(base_data[i,1])
  O<-data.frame(rbind(O,temp))
}
for(i in 1:nrow(O)){
  O$Name[i]<-paste(O$X.2[i],O$X.3[i],sep=' ')
}


# Start the plot
assign(paste('p',j,sep=''),
       ggplot(VIPSCORE, aes(x=as.factor(id), y=Y)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar

  # This add the bars with a blue color
  geom_bar(stat="identity", fill=Color_bar,alpha=0.7) +
    #geom_bar(aes(x=as.factor(id), y=Y, fill=Task), stat="identity", alpha=0.5) +
    #scale_fill_manual(values=c("blue","red")) +
  
  geom_segment(aes(x = c(1:nrow(VIPSCORE)), y = 0, xend = c(2:(nrow(VIPSCORE)+1)), yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(aes(x = c(1:nrow(VIPSCORE)), y = 1, xend = c(2:(nrow(VIPSCORE)+1)), yend = 1), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(aes(x = c(1:nrow(VIPSCORE)), y = 2, xend = c(2:(nrow(VIPSCORE)+1)), yend = 2), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(aes(x = c(1:nrow(VIPSCORE)), y = 3, xend = c(2:(nrow(VIPSCORE)+1)), yend = 3), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-1,6) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    #plot.margin = unit(rep(-1,4), "cm")      # Adjust the margin to make in sort labels are not truncated!
  ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  ggplot2::annotate("text", x = rep(max(VIPSCORE$id),4), y = c(0.2, 1.2, 2.2, 3.2), label = c("0", "1", "2", "3") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0) +
  
  # Add the labels, using the label_data dataframe that we have created before
  geom_text(data=O, aes(x=id, y=3.7, label=Name, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  geom_text(data=base_data, aes(x=c(start[1]+2,start[2],start[3:6]+3,start[7]+1,start[8]+2,start[9:11]+4,end[12:15]-2), y=3.5, label=X.1, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle=base_data$angle, inherit.aes = FALSE ) +
  geom_text(data=Membres, aes(x=1, y=6, label=eval(paste(Membres[j]))), color="black", fontface="bold",alpha=1, size=3.5, inherit.aes = FALSE ) +
  geom_segment(data=base_data, aes(x = start, y = 3.6, xend = end, yend = 3.6), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE )

)
ggsave(eval(parse(text=paste("p", j,sep=""))), file=paste("output", j, Task, ".svg",sep=""), width=10, height=10)
}





