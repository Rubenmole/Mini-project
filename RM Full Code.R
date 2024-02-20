install.packages("usdm")
install.packages("psych")
install.packages("vegan")
install.packages("factoextra")
install.packages("ggeffects")
install.packages("ggpubr")
install.packages("FactoMineR")
install.packages("lsmeans")
install.packages("lmerTest")
install.packages("dplyr")
install.packages("lme4")
install.packages("Matrix")
install.packages("lmerTest")
install.packages("ggplot2")
install.packages("ggdist")
install.packages("ggrepel")
library(lsmeans)
library(emmeans)
library(FactoMineR)
library(ggpubr)
library(usdm)
library(psych)
library(vegan)
library(factoextra)
library(ggeffects)
library(dplyr)
library(lme4)
library(Matrix)
library(lmerTest)
library(ggplot2)
library(ggdist)
library(ggrepel)

#Install all relevant packages


rm(list = ls()) #Clear environment

#Read in both data sets
statedata <- read.csv("RMStateData.csv" , header = T, sep = ",") 
d <- read.csv("RMFullData.csv", header = T, sep = ",")
head(d) # View the first few rows to make sure everything is correct



#Cleaning the data and running a PCA
d <- na.omit(d)
d$SVL..cm. <- as.numeric(d$SVL..cm.)#Making sure R knows this SVL is numeric
d$TL..cm <- as.numeric(d$TL..cm)#Making sure R knows total length is numeric
d <- d[is.finite(d$SVL..cm.) & is.finite(d$TL..cm) & is.finite(d$Prop.of.dorsum.red) & is.finite(d$Prop.of.dorsum.black),] #All data is finite
dpss <- princomp(cbind(d$SVL..cm., d$TL..cm, d$Prop.of.dorsum.red,d$Prop.of.dorsum.black ), cor = TRUE)# combining new data and its renamed
summary(dpss) # looks like 2 PCA are the correct amount
loadings(dpss)

#Lets Dduble check how many componants should be used.


eigenvalues <- dpss$sdev^2 #Double check and looks like 2 Components should still be selected 
eigenvalues # Comp 1 is above 2 and Comp 2 is very close to 2
plot(dpss, type="lines", ylim=c(0,2)) # Scree plot confirms 2 Components should be used
mnmn <- eigenvalues / sum(eigenvalues)
biplot(dpss) 
#Test PCA plot
fviz_pca_biplot(dpss, repel = TRUE, geom.ind = "point", ellipse.level=0.95, col.var = "black", labelsize=4,  ) #See what this looks like (basic plot)


#Isolating the species and Comps for the PCA plot, we know we are using Comp 1 and Comp 2. setting up for the PCA plot
biplot(dpss)
str(dpss)
dpss$scores
dpssc <- cbind(d,dpss$scores[,1:2])
biplot_data <- as.data.frame(dpss$scores[, 1:2])
biplot_data$Species <- d$Species




fviz_pca_biplot( dpss,
                 habillage = d$Species,  # Used to colour points by species
                 repel = TRUE, geom.ind = "point", ellipse.level = 0.95, # 95% ellipse cover 
                 col.var = "black",
                 labelsize = 4,
                 addEllipses = TRUE,
                 xlab = "Body Feature Lengths (PCA Component 1)", #Label 
                 ylab = "Colouration (PCA Component 2)",
                 geom = "ellipse" # Shape of geom
)






#PCA complete 


###############################################################################################################################
#Checking for significant differences using LMMs

#Isolating each species from the data set for LMM

d963 <- d %>% 
  filter(Species == "C. coccinea (East)") 
d304 <- d %>% 
  filter(Species == "M.tener")
d119 <- d %>% 
  filter(Species == "C. coccinea (West)")
d130 <- d %>% 
  filter(Species == "L. gentilis")

#Filtering for every species and creating new column where each of them has the priciple component of just there species

#M.tener and 4 of its mimics are all independently labelled 
#CWW
d5 <- na.omit(d119) #Removing NA from data set
d5$SVL..cm. <- as.numeric(d5$SVL..cm.)#Ensuring SVL is read as numeric.
d5$TL..cm <- as.numeric(d5$TL..cm)#Ensuring TL is read as numeric.
d5 <- d5[is.finite(d5$SVL..cm.) & is.finite(d5$TL..cm) & is.finite(d5$Prop.of.dorsum.red) & is.finite(d5$Prop.of.dorsum.black),]#All data is being read as finite.
dpss5 <- princomp(cbind(d5$SVL..cm., d5$TL..cm, d5$Prop.of.dorsum.red,d5$Prop.of.dorsum.black ), cor = TRUE) #Creating PCA variables for these species specifically.


#LG
d4 <- na.omit(d130)#Removing NA from data set
d4$SVL..cm. <- as.numeric(d4$SVL..cm.)#Ensuring SVL is read as numeric.
d4$TL..cm <- as.numeric(d4$TL..cm)#Ensuring TL is read as numeric.
d4 <- d4[is.finite(d4$SVL..cm.) & is.finite(d4$TL..cm) & is.finite(d4$Prop.of.dorsum.red) & is.finite(d4$Prop.of.dorsum.black),]# All data is being read as finite.
dpss4 <- princomp(cbind(d4$SVL..cm., d4$TL..cm, d4$Prop.of.dorsum.red,d4$Prop.of.dorsum.black ), cor = TRUE)#Creating PCA variables for these species specifically.


#CCE
d3 <- na.omit(d963) #Removing NA from data set
d3$SVL..cm. <- as.numeric(d3$SVL..cm.)#Ensuring SVL is read as numeric.
d3$TL..cm <- as.numeric(d3$TL..cm) #Ensuring TL is read as numeric.
d3 <- d3[is.finite(d3$SVL..cm.) & is.finite(d3$TL..cm) & is.finite(d3$Prop.of.dorsum.red) & is.finite(d3$Prop.of.dorsum.black),] #All data is being read as finite.
dpss3 <- princomp(cbind(d3$SVL..cm., d3$TL..cm, d3$Prop.of.dorsum.red,d3$Prop.of.dorsum.black ), cor = TRUE) #Creating PCA variables for these species specifically.


#MT 
d1 <- na.omit(d304) #Removing NA from data set
d1$SVL..cm. <- as.numeric(d1$SVL..cm.) #Ensuring SVL is read as numeric.
d1$TL..cm <- as.numeric(d1$TL..cm) #Ensuring TL is read as numeric.
d1 <- d1[is.finite(d1$SVL..cm.) & is.finite(d1$TL..cm) & is.finite(d1$Prop.of.dorsum.red) & is.finite(d1$Prop.of.dorsum.black),] #All data is being read as finite.
dpss1 <- princomp(cbind(d1$SVL..cm., d1$TL..cm, d1$Prop.of.dorsum.red,d1$Prop.of.dorsum.black ), cor = TRUE) #Creating PCA variables for these species specifically.



# LMM Comp 1  
# Comments below show the number of data points in each column to remind me which species is which.
#Isolating Comp one from each "dpss score".
n211 <- dpss1$scores[, 1] #304 
n231 <- dpss3$scores[, 1] #963
n241<- dpss4$scores[, 1] #130
n251<- dpss5$scores[, 1] #119

#Converting vectors to data frames

n211 <- data.frame(n211)
n231<- data.frame(n231)
n241 <- data.frame(n241)
n251 <- data.frame(n251)

# Creating column names with just species names

dpss1MT523 <- data.frame(species = rep("MT", 304))
dpss3CCE963 <- data.frame(species = rep("CCE", 963))
dpss4LG130<- data.frame(species = rep("LG", 130))
dpss5CWW119 <- data.frame(species = rep("CWW", 119))

colnames(n211) <- "common_column"
colnames(n231) <- "common_column"
colnames(n241) <- "common_column"
colnames(n251) <- "common_column"
#R required me to name these all the same.
dpss1data <- rbind(n211, n231, n241, n251) #Combining  the 4 species Comp 1 as individual columns within a data frame


colnames(dpss1MT523) <- "common_column"
colnames(dpss3CCE963) <- "common_column"
colnames(dpss4LG130) <- "common_column"
colnames(dpss5CWW119) <- "common_column"

names <- rbind(dpss1MT523, dpss3CCE963, dpss4LG130, dpss5CWW119 ) #combining all the coloums with names.


#Combining species names with Comp 1 PCA results.
dpss1full <- cbind(names,dpss1data )
#Now we add states 
statedata <- read.csv("thestatedata.csv")
newstate <- cbind(statedata$State, statedata$County)
colnames(newstate) <- c("State", "County")
finald1 <-cbind(dpss1full, newstate )
colnames(finald1) <- c("Sp", "dpss1", "County", "State") #In the order MT CCE LG CWW


#Running a LMM with the nested effects of state and county
model21 <- lmer(dpss1 ~ Sp + (1|State / County), data = finald1) #County nested in state
model211 <- lmer(dpss1 ~ Sp + (1|State), data = finald1) # A simpler model without county
lr_test_result <- anova(model211, model21, test = "Chisq") #Testing which model is better
lr_test_result # model21 is significantly better so county should remain in the model 

lsmeans_model <- lsmeans(model21, specs = "Sp")
summary(model21)


# Perform pairwise comparisons to double-check that there are no sig interactions 
pairwise_comparisons <- pairs(lsmeans_model, adjust = "tukey")
pairwise_comparisons


# LMM Comp 2  
# Comments below show the number of data points in each column to remind me which species is which.
#Isolating Comp two from each "dpss score".

#Isolating Comp from the data sets (Comp 2)
n210 <- dpss1$scores[, 2] #523
n230 <- dpss3$scores[, 2] #961
n240<- dpss4$scores[, 2] #130
n250<- dpss5$scores[, 2] #119

#Converitng them into data frames.
n21 <- data.frame(n210)
n23<- data.frame(n230)
n24 <- data.frame(n240)
n25 <- data.frame(n250)
#Creating a data set with species names
dpss1MT523 <- data.frame(species = rep("MT", 304))
dpss3CCE963 <- data.frame(species = rep("CCE", 963))
dpss4LG130<- data.frame(species = rep("LG", 130))
dpss5CWW119 <- data.frame(species = rep("CWW", 119))


colnames(n21) <- "common_column"
colnames(n23) <- "common_column"
colnames(n24) <- "common_column"
colnames(n25) <- "common_column"
#Renaming the column names or they can not be combined
dpss2data <- rbind(n21, n23, n24, n25)#Combining column names 


colnames(dpss1MT523) <- "common_column"
colnames(dpss3CCE963) <- "common_column"
colnames(dpss4LG130) <- "common_column"
colnames(dpss5CWW119) <- "common_column"
#Renaming the column names or they can not be combined

names <- rbind(dpss1MT523, dpss3CCE963, dpss4LG130, dpss5CWW119 ) #Combining column names

#Combining specific names with Comp 2 data
dpss2full <- cbind(names,dpss2data )
#Now adding states
statedata <- read.csv("thestatedata.csv")
newstate <- cbind(statedata$State, statedata$County)
colnames(newstate) <- c("State", "County")
finald <-cbind(dpss2full, newstate )
colnames(finald) <- c("Sp", "dpss2", "County", "State")#Ensuring they are names for further analysis.

#Create 2 LMMs
model211<- lmer(dpss2 ~ Sp + (1|State / County), data = finald) #First model with state and county nested
model222 <- lmer(dpss2 ~ Sp + (1|State), data = finald)#Simpler model that removes state

lr_test_result <- anova(model211, model222, test = "Chisq") #Comparing models (this method was used as Sp is categorical)
print(lr_test_result) #show results, it shows that there is no significant diffrence so we should use the simpler model


summary(model222)

lsmeans_model <- lsmeans(model222, specs = "Sp")

# Perform pairwise comparisons to confirm what the summary() showed
pairwise_comparisons2 <- pairs(lsmeans_model, adjust = "tukey")
pairwise_comparisons2

par(mfrow = c(1, 2)) #Setting up for a boxplot



#Final plots

#Comp 1 Boxplot
combined_data <- c(n251,n231,n241,n211)
# Create a box plot
boxplot(combined_data, col = c("white", "white", "white", "grey"),     names = c((italic("C. coccinea"))~"(West)", (italic("C. coccinea"))~"(East)", expression(italic("L.gentilis")), expression(italic("M.tener"))),    main = "Boxplot of PCA Component 1 Scores",        ylab = "Component 1 scores (Body Feature Lengths)",   xlab = expression(italic("Species")))



#Comp 2 Boxplot
combined_data <- c(n25,n23,n24,n21)
# Create a box plot
boxplot(combined_data, col = c("white", "white", "white", "grey"),   names = c((italic("C. coccinea"))~"(West)", (italic("C. coccinea"))~"(East)", expression(italic("L.gentilis")), expression(italic("M.tener"))),  main = "Boxplot of PCA Component 2 scores",  ylab = "Component 2 scores (Body colouration)",  xlab = expression(italic("Species")))



