rm(list=ls())

library(ggpubr)
library(tidyverse)
library(readxl)
library(ggrepel)

set.seed(429466026)

volcano_plots <- function(sample_type, text_size, data_type) {

# read in chemical data 
df_chems <- read.csv(paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/data/", data_type, "chemical_abundances_", sample_type,".csv"))
rownames(df_chems) <- df_chems$X 
df_chems$X <- NULL

chems <- colnames(df_chems %>% dplyr::select(-study_id))

# read in demographic data 

df_c <- read.csv(paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/data/demographics_", sample_type,".csv"))

# relevel factors to make sure reference categories are correct 
df_c$mat_race_eth <- factor(df_c$mat_race_eth, levels=c("White","Latina","Black","Asian/PI","Other or multi-race"))
df_c$mat_edu <- factor(df_c$mat_edu, levels=c("Less than high school","High school diploma or GED","Some college course work for credit or AA degree","Bachelor's degree (4 years)"))
df_c$employ_stat <- factor(df_c$employ_stat, levels=c("Not in labor force","Employed/Student","Unemployed/Other"))
df_c$occupation <- factor(df_c$occupation, levels=c("Not in labor force or unknown","Health care", "Research", 
                                                    "Business or legal", "Education, arts, or caregiving", "Cleaning, food, or other service"))

# want to see if chemical abundances differ based on maternal characteristics 
# main SES variables of interest: race/ethnicity, educational attainment, employment status 
# see if you can get insurance type? 
# other control variables: age? could also be on the pathway 
# considered controlling for gestational age at sample collection, butI think these are all delivery samples, so may be mediator


# after log10 transform, chemical appears approximately normal 

results_ls <- list()
pvals_ls <- list() 


for (i in 1:length(chems)) {
  
  #print(paste0(i,"/", length(chems)))
  # create dataframe with demographics and chemical of interest 
  dat <-  df_c 
  dat$study_id <- paste0(sample_type, dat$ppt_id)
  dat <- left_join(dat, df_chems %>% dplyr::select(study_id, chems[[i]]))
  
  # remove 1 person missing chemical data 
  dim(dat)
  dat <- dat %>% filter(!is.na(get(chems[[i]])))
  dim(dat)

  

m1 <- glm(log10(get(chems[[i]])) ~ mat_race_eth, data=dat)
m2 <- glm(log10(get(chems[[i]])) ~ mat_edu, data=dat)
m3 <- glm(log10(get(chems[[i]])) ~ employ_stat, data=dat)
m4 <- glm(log10(get(chems[[i]])) ~ hospital, data=dat)
m5 <- glm(log10(get(chems[[i]])) ~ age_cat, data=dat)
m6 <- glm(log10(get(chems[[i]])) ~ us_born, data=dat)
m7 <- glm(log10(get(chems[[i]])) ~ occupation, data=dat)
m8 <- glm(log10(get(chems[[i]])) ~ food_insecure, data=dat)
m9 <- glm(log10(get(chems[[i]])) ~ financial_strain, data=dat)


# make note of reference categories 
if (data_type!="aolin30_") {
results_ls[[i]] <- data.frame(cbind(chem_id = chems[[i]], 
                            Latina = m1$coefficients["mat_race_ethLatina"], 
                            Black = m1$coefficients["mat_race_ethBlack"], 
                            Asian = m1$coefficients["mat_race_ethAsian/PI"], 
                            Other = m1$coefficients["mat_race_ethOther or multi-race"], 
                            High_School = m2$coefficients["mat_eduHigh school diploma or GED"], 
                            Some_College = m2$coefficients["mat_eduSome college course work for credit or AA degree"], 
                            College_Grad = m2$coefficients["mat_eduBachelor's degree (4 years)"], 
                            Employed = m3$coefficients["employ_statEmployed/Student"], 
                            Unemployed = m3$coefficients["employ_statUnemployed/Other"], 
                            SFGH = m4$coefficients["hospitalSFGH"], 
                            #Unknown_Hospital = m4$coefficients["hospitalUnknown"], 
                            Age25_30 = m5$coefficients["age_cat(25,30]"], 
                            Age30_35 = m5$coefficients["age_cat(30,35]"], 
                            Age35_47 = m5$coefficients["age_cat(35,48]"], 
                            Not_US_born = m6$coefficients["us_bornNot born in US"], 
                            #US_born_missing = m6$coefficients["us_bornMissing"], 
                            Health_care = m7$coefficients["occupationHealth care"], 
                            Research = m7$coefficients["occupationResearch"], 
                            Service = m7$coefficients["occupationCleaning, food, or other service"], 
                            Education = m7$coefficients["occupationEducation, arts, or caregiving"], 
                            Business = m7$coefficients["occupationBusiness or legal"], 
                            Food_insecure = m8$coefficients["food_insecure"], 
                            Financial_strain = m9$coefficients["financial_strain"]))


pvals_ls[[i]] <- data.frame(cbind(chem_id = chems[[i]], 
                            Latina = summary(m1)$coefficients["mat_race_ethLatina",4], 
                            Black = summary(m1)$coefficients["mat_race_ethBlack",4], 
                            Asian = summary(m1)$coefficients["mat_race_ethAsian/PI",4], 
                            Other = summary(m1)$coefficients["mat_race_ethOther or multi-race",4], 
                            High_School = summary(m2)$coefficients["mat_eduHigh school diploma or GED",4], 
                            Some_College = summary(m2)$coefficients["mat_eduSome college course work for credit or AA degree",4], 
                            College_Grad = summary(m2)$coefficients["mat_eduBachelor's degree (4 years)",4], 
                            Employed = summary(m3)$coefficients["employ_statEmployed/Student",4], 
                            Unemployed = summary(m3)$coefficients["employ_statUnemployed/Other",4], 
                            SFGH = summary(m4)$coefficients["hospitalSFGH",4], 
                            #Unknown_Hospital = summary(m4)$coefficients["hospitalUnknown",4], 
                            Age25_30 = summary(m5)$coefficients["age_cat(25,30]", 4], 
                            Age30_35 = summary(m5)$coefficients["age_cat(30,35]",4], 
                            Age35_47 = summary(m5)$coefficients["age_cat(35,48]", 4], 
                            Not_US_born = summary(m6)$coefficients["us_bornNot born in US",4], 
                            #US_born_missing = summary(m6)$coefficients["us_bornMissing", 4], 
                            Health_care = summary(m7)$coefficients["occupationHealth care", 4], 
                            Research = summary(m7)$coefficients["occupationResearch", 4], 
                            Service = summary(m7)$coefficients["occupationCleaning, food, or other service", 4], 
                            Education = summary(m7)$coefficients["occupationEducation, arts, or caregiving", 4], 
                            Business = summary(m7)$coefficients["occupationBusiness or legal",4], 
                            Food_insecure = summary(m8)$coefficients["food_insecure", 4], 
                            Financial_strain = summary(m9)$coefficients["financial_strain", 4]))
}

else {
  results_ls[[i]] <- data.frame(cbind(chem_id = chems[[i]], 
                                      Latina = m1$coefficients["mat_race_ethLatina"], 
                                      Black = m1$coefficients["mat_race_ethBlack"], 
                                      Asian = m1$coefficients["mat_race_ethAsian/PI"], 
                                      #Other = m1$coefficients["mat_race_ethOther or multi-race"], 
                                      High_School = m2$coefficients["mat_eduHigh school diploma or GED"], 
                                      Some_College = m2$coefficients["mat_eduSome college course work for credit or AA degree"], 
                                      College_Grad = m2$coefficients["mat_eduBachelor's degree (4 years)"], 
                                      Employed = m3$coefficients["employ_statEmployed/Student"], 
                                      Unemployed = m3$coefficients["employ_statUnemployed/Other"], 
                                      SFGH = m4$coefficients["hospitalSFGH"], 
                                      #Unknown_Hospital = m4$coefficients["hospitalUnknown"], 
                                      Age25_30 = m5$coefficients["age_cat(25,30]"], 
                                      Age30_35 = m5$coefficients["age_cat(30,35]"], 
                                      Age35_47 = m5$coefficients["age_cat(35,48]"], 
                                      Not_US_born = m6$coefficients["us_bornNot born in US"], 
                                      #US_born_missing = m6$coefficients["us_bornMissing"], 
                                      Health_care = m7$coefficients["occupationHealth care"], 
                                      Research = m7$coefficients["occupationResearch"], 
                                      Service = m7$coefficients["occupationCleaning, food, or other service"], 
                                      Education = m7$coefficients["occupationEducation, arts, or caregiving"], 
                                      Business = m7$coefficients["occupationBusiness or legal"], 
                                      Food_insecure = m8$coefficients["food_insecure"], 
                                      Financial_strain = m9$coefficients["financial_strain"]))
  
  
  pvals_ls[[i]] <- data.frame(cbind(chem_id = chems[[i]], 
                                    Latina = summary(m1)$coefficients["mat_race_ethLatina",4], 
                                    Black = summary(m1)$coefficients["mat_race_ethBlack",4], 
                                    Asian = summary(m1)$coefficients["mat_race_ethAsian/PI",4], 
                                    #Other = summary(m1)$coefficients["mat_race_ethOther or multi-race",4], 
                                    High_School = summary(m2)$coefficients["mat_eduHigh school diploma or GED",4], 
                                    Some_College = summary(m2)$coefficients["mat_eduSome college course work for credit or AA degree",4], 
                                    College_Grad = summary(m2)$coefficients["mat_eduBachelor's degree (4 years)",4], 
                                    Employed = summary(m3)$coefficients["employ_statEmployed/Student",4], 
                                    Unemployed = summary(m3)$coefficients["employ_statUnemployed/Other",4], 
                                    SFGH = summary(m4)$coefficients["hospitalSFGH",4], 
                                    #Unknown_Hospital = summary(m4)$coefficients["hospitalUnknown",4], 
                                    Age25_30 = summary(m5)$coefficients["age_cat(25,30]", 4], 
                                    Age30_35 = summary(m5)$coefficients["age_cat(30,35]",4], 
                                    Age35_47 = summary(m5)$coefficients["age_cat(35,48]", 4], 
                                    Not_US_born = summary(m6)$coefficients["us_bornNot born in US",4], 
                                    #US_born_missing = summary(m6)$coefficients["us_bornMissing", 4], 
                                    Health_care = summary(m7)$coefficients["occupationHealth care", 4], 
                                    Research = summary(m7)$coefficients["occupationResearch", 4], 
                                    Service = summary(m7)$coefficients["occupationCleaning, food, or other service", 4], 
                                    Education = summary(m7)$coefficients["occupationEducation, arts, or caregiving", 4], 
                                    Business = summary(m7)$coefficients["occupationBusiness or legal",4], 
                                    Food_insecure = summary(m8)$coefficients["food_insecure", 4], 
                                    Financial_strain = summary(m9)$coefficients["financial_strain", 4]))  
  
}
}

results <- data.frame(do.call(rbind, results_ls))
rownames(results) <- results$chem_id
pvals <- data.frame(do.call(rbind, pvals_ls))
rownames(pvals) <- pvals$chem_id


annotations <- read.csv(paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/data/", data_type, "chemical_annotations_", sample_type,".csv"))

results <- left_join(results, annotations, by="chem_id")

results <- results %>% rowwise() %>% mutate(Sum = sum(c(Pharma_Presence, Pest_Presence, PFAS_Presence, Plast_Presence, Cosm_Presence)))
#Flame_Presence, -- dont' have this in Aolin's data 
table(results$Sum)

# plots by race ------------------------------------------------------------------------------------------------ 

results$Latina <- as.numeric(as.character(results$Latina))
results$Black <- as.numeric(as.character(results$Black))
results$Asian <- as.numeric(as.character(results$Asian))
if (data_type!="aolin30_"){
results$Other <- as.numeric(as.character(results$Other))
r_results <- results %>% dplyr::select(chem_id, Latina, Black, Asian, Other)
r_results <- pivot_longer(r_results, cols=Latina:Other, names_to = "Race", values_to = "Diff_from_ref")

}
if (data_type=="aolin30_"){
r_results <- results %>% dplyr::select(chem_id, Latina, Black, Asian)
r_results <- pivot_longer(r_results, cols=Latina:Asian, names_to = "Race", values_to = "Diff_from_ref")
}
r_results <- left_join(r_results, annotations, by="chem_id")

# pvalues 
if (data_type!="aolin30_"){
  r_pvals <- pvals %>% dplyr::select(chem_id, Latina, Black, Asian, Other)
  r_pvals<- pivot_longer(r_pvals, cols=Latina:Other, names_to = "Race", values_to = "pvalue")
}
if (data_type=="aolin30_"){
r_pvals <- pvals %>% dplyr::select(chem_id, Latina, Black, Asian)
r_pvals<- pivot_longer(r_pvals, cols=Latina:Asian, names_to = "Race", values_to = "pvalue")
}

r_results <- left_join(r_results, r_pvals, by=c("chem_id","Race"))

# plot only those that are significant 
r_results$pvalue <- as.numeric(as.character(r_results$pvalue))
r_results$FDR_pvalue <- p.adjust(r_results$pvalue, method="BH")

# plots by education ------------------------------------------------------------------------------------------------ 

results$High_School <- as.numeric(as.character(results$High_School))
results$Some_College <- as.numeric(as.character(results$Some_College))
results$College_Grad <- as.numeric(as.character(results$College_Grad))

ed_results <- results %>% dplyr::select(chem_id, High_School, Some_College, College_Grad)


ed_results <- pivot_longer(ed_results, cols=High_School:College_Grad, names_to = "Education", values_to = "Diff_from_ref")

ed_results <- left_join(ed_results, annotations, by="chem_id")

# pvalues 
ed_pvals <- pvals %>% dplyr::select(chem_id, High_School, Some_College, College_Grad)
ed_pvals<- pivot_longer(ed_pvals, cols=High_School:College_Grad, names_to = "Education", values_to = "pvalue")

ed_results <- left_join(ed_results, ed_pvals, by=c("chem_id","Education"))

# plot only those that are significant 
ed_results$pvalue <- as.numeric(as.character(ed_results$pvalue))
ed_results$FDR_pvalue <- p.adjust(ed_results$pvalue, method="BH")

# plots by employment status  ------------------------------------------------------------------------------------------------ 

results$Employed <- as.numeric(as.character(results$Employed))
results$Unemployed <- as.numeric(as.character(results$Unemployed))

emp_results <- results %>% dplyr::select(chem_id, Employed, Unemployed)


emp_results <- pivot_longer(emp_results, cols=Employed:Unemployed, names_to = "Employment", values_to = "Diff_from_ref")

emp_results <- left_join(emp_results, annotations, by="chem_id")

# pvalues 
emp_pvals <- pvals %>% dplyr::select(chem_id, Employed, Unemployed)


emp_pvals<- pivot_longer(emp_pvals, cols=Employed:Unemployed, names_to = "Employment", values_to = "pvalue")

emp_results <- left_join(emp_results, emp_pvals, by=c("chem_id","Employment"))

# plot only those that are significant 
emp_results$pvalue <- as.numeric(as.character(emp_results$pvalue))
emp_results$FDR_pvalue <- p.adjust(emp_results$pvalue, method="BH")

# plots by hospital ------------------------------------------------------------------------------------------------ 

results$SFGH <- as.numeric(as.character(results$SFGH))
#results$Unknown_Hospital <- as.numeric(as.character(results$Unknown_Hospital))

#h_results <- results %>% select(chem_id, SFGH, Unknown_Hospital)
h_results <- results %>% dplyr::select(chem_id, SFGH)


#h_results <- pivot_longer(h_results, cols=SFGH:Unknown_Hospital, names_to = "Hospital", values_to = "Diff_from_ref")
h_results <- pivot_longer(h_results, cols=SFGH, names_to = "Hospital", values_to = "Diff_from_ref")

h_results <- left_join(h_results, annotations, by="chem_id")

# pvalues 
#h_pvals <- pvals %>% select(chem_id, SFGH, Unknown_Hospital)
h_pvals <- pvals %>% dplyr::select(chem_id, SFGH)


#h_pvals<- pivot_longer(h_pvals, cols=SFGH:Unknown_Hospital, names_to = "Hospital", values_to = "pvalue")
h_pvals<- pivot_longer(h_pvals, cols=SFGH, names_to = "Hospital", values_to = "pvalue")

h_results <- left_join(h_results, h_pvals, by=c("chem_id","Hospital"))

# plot only those that are significant 
h_results$pvalue <- as.numeric(as.character(h_results$pvalue))
h_results$FDR_pvalue <- p.adjust(h_results$pvalue, method="BH")

# plots by age ------------------------------------------------------------------------------------------------ 

results$Age25_30 <- as.numeric(as.character(results$Age25_30))
results$Age30_35 <- as.numeric(as.character(results$Age30_35))
results$Age35_47 <- as.numeric(as.character(results$Age35_47))

a_results <- results %>% dplyr::select(chem_id, Age25_30, Age30_35, Age35_47)


a_results <- pivot_longer(a_results, cols=Age25_30:Age35_47, names_to = "Age", values_to = "Diff_from_ref")

a_results <- left_join(a_results, annotations, by="chem_id")

# pvalues 
a_pvals <- pvals %>% dplyr::select(chem_id, Age25_30, Age30_35, Age35_47)
a_pvals<- pivot_longer(a_pvals, cols=Age25_30:Age35_47, names_to = "Age", values_to = "pvalue")

a_results <- left_join(a_results, a_pvals, by=c("chem_id","Age"))

# plot only those that are significant 
a_results$pvalue <- as.numeric(as.character(a_results$pvalue))
a_results$FDR_pvalue <- p.adjust(a_results$pvalue, method="BH")

# plots by nativity ------------------------------------------------------------------------------------------------ 

results$Not_US_born <- as.numeric(as.character(results$Not_US_born))
#results$US_born_missing <- as.numeric(as.character(results$US_born_missing))

#n_results <- results %>% select(chem_id, Not_US_born, US_born_missing)
n_results <- results %>% dplyr::select(chem_id, Not_US_born)


#n_results <- pivot_longer(n_results, cols=Not_US_born:US_born_missing, names_to = "Nativity", values_to = "Diff_from_ref")
n_results <- pivot_longer(n_results, cols=Not_US_born, names_to = "Nativity", values_to = "Diff_from_ref")

n_results <- left_join(n_results, annotations, by="chem_id")

# pvalues 
#n_pvals <- pvals %>% select(chem_id, Not_US_born, US_born_missing)
n_pvals <- pvals %>% dplyr::select(chem_id, Not_US_born)


#n_pvals<- pivot_longer(n_pvals, cols=Not_US_born:US_born_missing, names_to = "Nativity", values_to = "pvalue")
n_pvals<- pivot_longer(n_pvals, cols=Not_US_born, names_to = "Nativity", values_to = "pvalue")

n_results <- left_join(n_results, n_pvals, by=c("chem_id","Nativity"))

# plot only those that are significant 
n_results$pvalue <- as.numeric(as.character(n_results$pvalue))
n_results$FDR_pvalue <- p.adjust(n_results$pvalue, method="BH")


# plots by occupation ------------------------------------------------------------------------------------------------ 

results$Health_care <- as.numeric(as.character(results$Health_care))
results$Research <- as.numeric(as.character(results$Research))
results$Service <- as.numeric(as.character(results$Service))
results$Education <- as.numeric(as.character(results$Education))
results$Business <- as.numeric(as.character(results$Business))

o_results <- results %>% dplyr::select(chem_id, Health_care, Research, Service, Education, Business)


o_results <- pivot_longer(o_results, cols=Health_care:Business, names_to = "Occupation", values_to = "Diff_from_ref")

o_results <- left_join(o_results, annotations, by="chem_id")

# pvalues 
o_pvals <- pvals %>% dplyr::select(chem_id, Health_care, Research, Service, Education, Business)


o_pvals<- pivot_longer(o_pvals, cols=Health_care:Business, names_to = "Occupation", values_to = "pvalue")

o_results <- left_join(o_results, o_pvals, by=c("chem_id","Occupation"))

# plot only those that are significant 
o_results$pvalue <- as.numeric(as.character(o_results$pvalue))
o_results$FDR_pvalue <- p.adjust(o_results$pvalue, method="BH")



# plots by food insecurity ------------------------------------------------------------------------------------------------ 

results$Food_insecure <- as.numeric(as.character(results$Food_insecure))

fi_results <- results %>% dplyr::select(chem_id, Food_insecure)


fi_results <- pivot_longer(fi_results, cols=Food_insecure, names_to = "Food_Security", values_to = "Diff_from_ref")

fi_results <- left_join(fi_results, annotations, by="chem_id")

# pvalues 
fi_pvals <- pvals %>% dplyr::select(chem_id, Food_insecure)


fi_pvals<- pivot_longer(fi_pvals, cols=Food_insecure, names_to = "Food_Security", values_to = "pvalue")

fi_results <- left_join(fi_results, fi_pvals, by=c("chem_id","Food_Security"))

# plot only those that are significant 
fi_results$pvalue <- as.numeric(as.character(fi_results$pvalue))
fi_results$FDR_pvalue <- p.adjust(fi_results$pvalue, method="BH")



# plots by financial strain ------------------------------------------------------------------------------------------------ 

results$Financial_strain <- as.numeric(as.character(results$Financial_strain))

fs_results <- results %>% dplyr::select(chem_id, Financial_strain)


fs_results <- pivot_longer(fs_results, cols=Financial_strain, names_to = "Financial_Strain", values_to = "Diff_from_ref")

fs_results <- left_join(fs_results, annotations, by="chem_id")

# pvalues 
fs_pvals <- pvals %>% dplyr::select(chem_id, Financial_strain)


fs_pvals<- pivot_longer(fs_pvals, cols=Financial_strain, names_to = "Financial_Strain", values_to = "pvalue")

fs_results <- left_join(fs_results, fs_pvals, by=c("chem_id","Financial_Strain"))

# plot only those that are significant 
fs_results$pvalue <- as.numeric(as.character(fs_results$pvalue))
fs_results$FDR_pvalue <- p.adjust(fs_results$pvalue, method="BH")




# plots 

p1 <- ggplot(a_results) + geom_point(aes(x=Diff_from_ref, y=-log(FDR_pvalue), color=factor(Age), shape=factor(Age)), size=4) + 
  scale_color_manual(name="", labels=c("26-30 vs. 15-25", "31-35 vs. 15-25", "36-47 vs. 15-26"), 
                     values=c("#7fcdbb","#1d91c0", "#0c2c84")) + 
  scale_shape_manual(name="", labels=c("26-30 vs. 15-25", "31-35 vs. 15-25", "36-47 vs. 15-26"), 
                     values=c(0, 1, 2)) +
  theme_bw() + lims(y=c(0,10), x=c(-1, 1)) + theme(legend.position = "top", legend.text = element_text(size=7)) + 
  geom_hline(yintercept = -log(0.1), linetype=3) +
  labs(y = "-ln(FDR corrected p-value)", x=expression(paste("difference in ", log[10]," abundance from reference group",sep=""))) + 
  geom_text_repel(data=a_results %>% filter(-log(FDR_pvalue)> -log(0.1)), nudge_y = 0.75,
                  aes(x=Diff_from_ref, y=-log(FDR_pvalue)), label=a_results$PREFERRED_NAME[a_results$FDR_pvalue<0.1], size=text_size)
  
 #p1
#+ scale_y_continuous(trans="reverse")

p2 <- ggplot(r_results) + geom_point(aes(x=Diff_from_ref, y=-log(FDR_pvalue), color=factor(Race), shape=factor(Race)), size=4) + 
  scale_color_manual(name="", labels=c("Asian vs. White","Black vs. White","Latina vs. White", "Other vs. White"), 
                     values=c("#d73027","#ff7f00","#66bd63","#006837")) + 
  scale_shape_manual(name="", labels=c("Asian vs. White","Black vs. White","Latina vs. White", "Other vs. White"), 
                     values=c(3,4,5,6)) + 
  theme_bw() + lims(y=c(0,10), x=c(-1, 1)) + theme(legend.position = "top", legend.text = element_text(size=7)) + 
  geom_hline(yintercept = -log(0.1), linetype=3) +
  labs(y = "-ln(FDR corrected p-value)", x=expression(paste("difference in ", log[10]," abundance from reference group",sep=""))) + 
  geom_text_repel(data=r_results %>% filter(-log(FDR_pvalue)> -log(0.1)), nudge_y = 0.75,
                  aes(x=Diff_from_ref, y=-log(FDR_pvalue)), label=r_results$PREFERRED_NAME[r_results$FDR_pvalue<0.1], size=text_size)


#p2
p3 <- ggplot(ed_results) + geom_point(aes(x=Diff_from_ref, y=-log(FDR_pvalue), color=factor(Education), shape=factor(Education)), size=4) + 
  scale_color_manual(name="", labels=c("College graduate vs. \nLess than high school", 
                                                   "High school graduate vs. \nLess than high school", 
                                                   "Some college vs. \nLess than high school"), 
                     values=c("#fecc5c","#f03b20", "#bd0026")) + 
  scale_shape_manual(name="", labels=c("College graduate vs. \nLess than high school", 
                                       "High school graduate vs. \nLess than high school", 
                                       "Some college vs. \nLess than high school"), 
                     values=c(7,8,9)) + 
  theme_bw() + lims(y=c(0,10), x=c(-1, 1)) + theme(legend.position = "top", legend.text = element_text(size=7)) + 
  geom_hline(yintercept = -log(0.1), linetype=3) +
  labs(y = "-ln(FDR corrected p-value)", x=expression(paste("difference in ", log[10]," abundance from reference group",sep=""))) +
  geom_text_repel(data=ed_results %>% filter(-log(FDR_pvalue)> -log(0.1)), nudge_y = 0.75,
                  aes(x=Diff_from_ref, y=-log(FDR_pvalue)), label=ed_results$PREFERRED_NAME[ed_results$FDR_pvalue<0.1], size=text_size)



#p3
p4 <- ggplot(emp_results) + geom_point(aes(x=Diff_from_ref, y=-log(FDR_pvalue), color=factor(Employment), shape=factor(Employment)), size=4) + 
  scale_color_manual(name="", 
                     labels=c("Not in labor force vs.\n Employed or student", "Not in labor force vs.\n Unemployed or other"), 
                     values=c("#8c6bb1", "#6e016b")) + 
  scale_shape_manual(name="", 
                     labels=c("Not in labor force vs.\n Employed or student", "Not in labor force vs.\n Unemployed or other"), 
                     values=c(10,11)) + 
  theme_bw() + lims(y=c(0,10), x=c(-1, 1)) + theme(legend.position = "top", legend.text = element_text(size=7)) + 
  geom_hline(yintercept = -log(0.1), linetype=3) +
  labs(y = "-ln(FDR corrected p-value)", x=expression(paste("difference in ", log[10]," abundance from reference group",sep=""))) + 
  geom_text_repel(data=emp_results %>% filter(-log(FDR_pvalue)> -log(0.1)), nudge_y = 0.75,
                  aes(x=Diff_from_ref, y=-log(FDR_pvalue)), label=emp_results$PREFERRED_NAME[emp_results$FDR_pvalue<0.1], size=text_size)

#p4
p5 <- ggplot(n_results) + geom_point(aes(x=Diff_from_ref, y=-log(FDR_pvalue), color=factor(Nativity), shape=factor(Nativity)), size=4) + 
  scale_color_manual(name="", 
                     labels=c("Not born in US vs. Born in US", "Missing vs. Born in US"), 
                     values=c('#fa9fb5','#c51b8a')) + 
  scale_shape_manual(name="", 
                     labels=c("Not born in US vs. Born in US", "Missing vs. Born in US"), 
                     values=c(12,13)) + 
  theme_bw() + lims(y=c(0,10), x=c(-1, 1)) + theme(legend.position = "top", legend.text = element_text(size=7)) + 
  geom_hline(yintercept = -log(0.1), linetype=3) +
  labs(y = "-ln(FDR corrected p-value)", x=expression(paste("difference in ", log[10]," abundance from reference group",sep=""))) + 
  geom_text_repel(data=n_results %>% filter(-log(FDR_pvalue)> -log(0.1)), nudge_y = 0.75,
                  aes(x=Diff_from_ref, y=-log(FDR_pvalue)), label=n_results$PREFERRED_NAME[n_results$FDR_pvalue<0.1], size=text_size)


#p5


p6 <- ggplot(o_results) + geom_point(aes(x=Diff_from_ref, y=-log(FDR_pvalue), color=factor(Occupation), shape=factor(Occupation)), size=3) + 
  scale_color_manual(name="", 
                     labels=c("Health care vs. \nNot in labor force", "Research vs. \nNot in labor force", "Food, cleaning, \nor other service vs. \nNot in labor force", "Education, arts, \nor childcare vs. \nNot in labor force", "Business or legal vs. \nNot in labor force"), 
                     values=c('#b2e2e2','#66c2a4','#2ca25f','#006d2c','#00441b')) + 
  scale_shape_manual(name="", 
                     labels=c("Health care vs. \nNot in labor force", "Research vs. \nNot in labor force", "Food, cleaning, \nor other service vs. \nNot in labor force", "Education, arts, \nor childcare vs. \nNot in labor force", "Business or legal vs. \nNot in labor force"), 
                     values=c(15,16,17,18, 19)) + 
  theme_bw() + lims(y=c(0,10), x=c(-1, 1)) + theme(legend.position = "top", legend.text = element_text(size=5)) + 
  geom_hline(yintercept = -log(0.1), linetype=3) +
  labs(y = "-ln(FDR corrected p-value)", x=expression(paste("difference in ", log[10]," abundance from reference group",sep=""))) + 
  geom_text_repel(data=o_results %>% filter(-log(FDR_pvalue)> -log(0.1)), nudge_y = 0.75,
                  aes(x=Diff_from_ref, y=-log(FDR_pvalue)), label=o_results$PREFERRED_NAME[o_results$FDR_pvalue<0.1], size=text_size)


#p6

p7 <- ggplot(h_results %>% filter(Hospital=="SFGH")) + geom_point(aes(x=Diff_from_ref, y=-log(FDR_pvalue), color=factor(Hospital), shape=factor(Hospital)), size=4) + 
  scale_color_manual(name="", 
                     labels=c("SFGH vs. Mission Bay"), 
                     values=c('#cb181d')) + 
  scale_shape_manual(name="", 
                     labels=c("SFGH vs. Mission Bay"), 
                     values=c(14)) + 
  theme_bw() + lims(y=c(0,10), x=c(-1, 1)) + theme(legend.position = "top", legend.text = element_text(size=7)) + 
  geom_hline(yintercept = -log(0.1), linetype=3) +
  labs(y = "-ln(FDR corrected p-value)", x=expression(paste("difference in ", log[10]," abundance from reference group",sep=""))) + 
  geom_text_repel(data=h_results %>% filter(-log(FDR_pvalue)> -log(0.1)), nudge_y = 0.75,
                  aes(x=Diff_from_ref, y=-log(FDR_pvalue)), label=h_results$PREFERRED_NAME[h_results$FDR_pvalue<0.1], size=text_size)



p8 <- ggplot(fi_results) + geom_point(aes(x=Diff_from_ref, y=-log(FDR_pvalue), color=factor(Food_Security), shape=factor(Food_Security)), size=4) + 
  scale_color_manual(name="", 
                     labels=c("Food insecure vs. Food secure"), 
                     values=c('#238443')) + 
  scale_shape_manual(name="", 
                     labels=c("Food insecure vs. Food secure"), 
                     values=c(24)) + 
  theme_bw() + lims(y=c(0,10), x=c(-1, 1)) + theme(legend.position = "top", legend.text = element_text(size=7)) + 
  geom_hline(yintercept = -log(0.1), linetype=3) +
  labs(y = "-ln(FDR corrected p-value)", x=expression(paste("difference in ", log[10]," abundance from reference group",sep=""))) + 
  geom_text_repel(data=fi_results %>% filter(-log(FDR_pvalue)> -log(0.1)), nudge_y = 0.75,
                  aes(x=Diff_from_ref, y=-log(FDR_pvalue)), label=fi_results$PREFERRED_NAME[fi_results$FDR_pvalue<0.1], size=text_size)


p9 <- ggplot(fs_results) + geom_point(aes(x=Diff_from_ref, y=-log(FDR_pvalue), color=factor(Financial_Strain), shape=factor(Financial_Strain)), size=4) + 
  scale_color_manual(name="", 
                     labels=c("Financial strain vs. No financial strain"), 
                     values=c('#66c2a4')) + 
  scale_shape_manual(name="", 
                     labels=c("Financial strain vs. No financial strain"), 
                     values=c(23)) + 
  theme_bw() + lims(y=c(0,10), x=c(-1, 1)) + theme(legend.position = "top", legend.text = element_text(size=7)) + 
  geom_hline(yintercept = -log(0.1), linetype=3) +
  labs(y = "-ln(FDR corrected p-value)", x=expression(paste("difference in ", log[10]," abundance from reference group",sep=""))) + 
  geom_text_repel(data=fs_results %>% filter(-log(FDR_pvalue)> -log(0.1)), nudge_y = 0.75,
                  aes(x=Diff_from_ref, y=-log(FDR_pvalue)), label=fs_results$PREFERRED_NAME[fs_results$FDR_pvalue<0.1], size=text_size)

#p7
mplot <- ggarrange(p1,p3,p4,p5,p6,p8,p9, ncol=2, nrow=5)
#mplot

ggsave(mplot, file=paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/results/plots/",data_type,"volcano_plots_annotated_color_",sample_type,"_7_names.pdf"), width=10, height=14)

# save results for supplement 

# age
a_results$sample_type <- sample_type 

df_a <- a_results %>% filter(FDR_pvalue<0.1)
write.csv(df_a, file=paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/results/age_diff_volcano_chems_",sample_type,"_FDR01.csv"))

write.csv(a_results, file=paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/results/age_diff_volcano_chems_",sample_type,".csv"))

#education 
ed_results$sample_type <- sample_type 

df_ed <- ed_results %>% filter(FDR_pvalue<0.1)
write.csv(df_ed, file=paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/results/edu_diff_volcano_chems_",sample_type,"_FDR01.csv"))

write.csv(ed_results, file=paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/results/edu_diff_volcano_chems_",sample_type,".csv"))


#employment 
emp_results$sample_type <- sample_type 

df_emp <- emp_results %>% filter(FDR_pvalue<0.1)
write.csv(df_emp, file=paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/results/emp_diff_volcano_chems_",sample_type,"_FDR01.csv"))

write.csv(emp_results, file=paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/results/emp_diff_volcano_chems_",sample_type,".csv"))


#nativity 
n_results$sample_type <- sample_type 

df_n <- n_results %>% filter(FDR_pvalue<0.1)
write.csv(df_n, file=paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/results/nativity_diff_volcano_chems_",sample_type,"_FDR01.csv"))

write.csv(n_results, file=paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/results/nativity_diff_volcano_chems_",sample_type,".csv"))


#occupation 
o_results$sample_type <- sample_type 

df_o <- o_results %>% filter(FDR_pvalue<0.1)
write.csv(df_o, file=paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/results/occupation_diff_volcano_chems_",sample_type,"_FDR01.csv"))

write.csv(o_results, file=paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/results/occupation_diff_volcano_chems_",sample_type,".csv"))


#food insecurity 
fi_results$sample_type <- sample_type 

df_fi <- fi_results %>% filter(FDR_pvalue<0.1)
write.csv(df_fi, file=paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/results/food_diff_volcano_chems_",sample_type,"_FDR01.csv"))

write.csv(fi_results, file=paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/results/food_diff_volcano_chems_",sample_type,".csv"))

#financial strain
fs_results$sample_type <- sample_type 

df_fs <- fs_results %>% filter(FDR_pvalue<0.1)
write.csv(df_fs, file=paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/results/finance_diff_volcano_chems_",sample_type,"_FDR01.csv"))

write.csv(fs_results, file=paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/results/finance_diff_volcano_chems_",sample_type,".csv"))

}



volcano_plots("M", 3, "")
volcano_plots("C", 3, "")


