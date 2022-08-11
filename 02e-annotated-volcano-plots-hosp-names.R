rm(list=ls())

library(tidyverse)
library(readxl)
library(ggrepel)
library(ggpubr)

set.seed(429466026)

hosp_volcano_plots <- function(sample_type, text_size) {
  
  # read in chemical data 
  df_chems <- read.csv(paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/data/chemical_abundances_", sample_type,".csv"))
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
  
  results <- data.frame(do.call(rbind, results_ls))
  rownames(results) <- results$chem_id
  pvals <- data.frame(do.call(rbind, pvals_ls))
  rownames(pvals) <- pvals$chem_id
  
  
  annotations <- read.csv(paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/data/chemical_annotations_", sample_type,".csv"))
  
  results <- left_join(results, annotations, by="chem_id")
  
  results <- results %>% rowwise() %>% mutate(Sum = sum(c(Pharma_Presence, Pest_Presence, Flame_Presence, PFAS_Presence, Plast_Presence, Cosm_Presence)))
  table(results$Sum)
  
 
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
  
 
  # plots 

  
  p1 <- ggplot(h_results %>% filter(Hospital=="SFGH")) + geom_point(aes(x=Diff_from_ref, y=-log(FDR_pvalue), color=factor(Hospital), shape=factor(Hospital)), size=5) + 
    scale_color_manual(name="", 
                       labels=c("SFGH vs. Mission Bay"), 
                       values=c('#cb181d')) + 
    scale_shape_manual(name="", 
                       labels=c("SFGH vs. Mission Bay"), 
                       values=c(14)) + 
    theme_bw() + lims(y=c(0,10), x=c(-1, 1)) + theme(legend.position = "none", axis.text = element_text(size=16), legend.text = element_text(size=16), axis.title = element_text(size=16)) +
    geom_hline(yintercept = -log(0.1), linetype=3) +
    labs(y = "-ln(FDR corrected p-value)", x=expression(paste("difference in ", log[10]," abundance from reference group",sep=""))) + 
    geom_text_repel(data=h_results %>% filter(Hospital=="SFGH" & -log(FDR_pvalue)> -log(0.1)), nudge_y = 0.75,
                    aes(x=Diff_from_ref, y=-log(FDR_pvalue)), label=h_results$PREFERRED_NAME[h_results$FDR_pvalue<0.1 & h_results$Hospital=="SFGH"], size=text_size, max.overlaps = 20)
  
  
  ggsave(p1, file=paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/results/plots/hosp_volcano_plots_annotated_color_",sample_type,"_names.pdf"), width=10)

  h_results$sample_type <- sample_type 
  
  df_h <- h_results %>% filter(FDR_pvalue<0.1)
  write.csv(df_h, file=paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/results/hosp_diff_volcano_chems_",sample_type,"_FDR01.csv"))
  
  write.csv(h_results, file=paste0("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/results/hosp_diff_volcano_chems_",sample_type,".csv"))
  
  }



hosp_volcano_plots("M", 4)
hosp_volcano_plots("C", 4)



df_h1 <- read.csv("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/results/hosp_diff_volcano_chems_M.csv")
df_h2 <- read.csv("/Users/danagoin/Documents/Research projects/R01 GSS New Methods/results/hosp_diff_volcano_chems_C.csv")

df_h <- data.frame(rbind(df_h1, df_h2))
df_h <- df_h %>% group_by(chem_id) %>% mutate(n_overlap = n() -1)
