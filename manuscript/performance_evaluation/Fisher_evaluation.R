
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(grid)
#install.packages("devtools")
#library(devtools)
#devtools::install_github("kassambara/easyGgplot2")
library(easyGgplot2)

rm(list=ls())

########################################################
#Section 1: collect the data from Centrifuge and MicroFisher
#setwd("F:/PostDoc_dataset/bioinfomatic/performance evaluation")
setwd("E:/PostDoc_dataset/bioinfomatic/performance evaluation/classification_result_collection_weight_length")
list.files(getwd())


communitys=c("simulating_50species_1", "simulating_50species_2", "simulating_50species_3", "simulating_50species_4", "simulating_50species_5",
            "simulating_100species_1", "simulating_100species_2", "simulating_100species_3", "simulating_100species_4", "simulating_100species_5",
            "simulating_200species_1", "simulating_200species_2", "simulating_200species_3", "simulating_200species_4", "simulating_200species_5")

#calculate the real community abundance
dataFrame = data.frame()
for (community in communitys){
  #real_community = read_tsv(paste(community, ".taxonomy.txt", sep = ""))
  real_community = read.table(paste(community, ".taxonomy.txt", sep = ""), fill = T)
  real_community_1=data.frame(seqID=real_community$V1, TaxID=real_community$V2, taxonomy=paste(real_community$V3, real_community$V4, sep = " "))
  real_community_2=real_community_1$taxonomy
  real_community_3=strsplit(real_community_2, ";")
  real_community_4=data.frame(t(data.frame(real_community_3)))
  real_community_4[,1]=str_replace(real_community_4[,1], "k__", "")
  real_community_4[,2]=str_replace(real_community_4[,2], "p__", "")
  real_community_4[,3]=str_replace(real_community_4[,3], "c__", "")
  real_community_4[,4]=str_replace(real_community_4[,4], "o__", "")
  real_community_4[,5]=str_replace(real_community_4[,5], "f__", "")
  real_community_4[,6]=str_replace(real_community_4[,6], "g__", "")
  real_community_4[,7]=str_replace(real_community_4[,7], "s__", "")
  real_community_5=data.frame(Kingdom=real_community_4[,1], Phylum=real_community_4[,2], Class=real_community_4[,3], Order=real_community_4[,4],
                              Family=real_community_4[,5], Genus=real_community_4[,6], Species=real_community_4[,7]) %>% dplyr::distinct()
  
  real_abundance=read.table(paste(community, ".short_read_abundance.txt", sep = ""))
  real_abundance[,1]=str_replace(real_abundance[,1], paste(paste("/home/microbiome/data_storage/SATA2/Fisher_test/Test_Fungi_RefSeq/mock_community/", community, sep = ""),
                                                           "/splite_seq/", sep = ""), "")
  real_abundance[,1]=str_replace(real_abundance[,1], ".fasta", "")
  real_abundance[,1]=str_replace(real_abundance[,1], "_", " ")
  colnames(real_abundance) = c("Species", "Abundance")
  
  real_community_abundance = real_community_5 %>% dplyr::left_join(real_abundance, by = "Species")
  real_community_abundance = real_community_abundance %>% dplyr::mutate(Mock_community = community)  %>% dplyr::distinct()

    
  ter_K<-aggregate(real_community_abundance[,c("Abundance")] ~ Kingdom, data = real_community_abundance, sum) %>% mutate(Rank = "kingdom")
  colnames(ter_K) <- c("Rank_name", "Abundance", "Rank")
  ter_P<-aggregate(real_community_abundance[,c("Abundance")] ~ Phylum, data = real_community_abundance, sum) %>% mutate(Rank = "phylum")
  colnames(ter_P) <- c("Rank_name", "Abundance", "Rank")
  ter_C<-aggregate(real_community_abundance[,c("Abundance")] ~ Class, data = real_community_abundance, sum) %>% mutate(Rank = "class")
  colnames(ter_C) <- c("Rank_name", "Abundance", "Rank")
  ter_O<-aggregate(real_community_abundance[,c("Abundance")] ~ Order, data = real_community_abundance, sum) %>% mutate(Rank = "order")
  colnames(ter_O) <- c("Rank_name", "Abundance", "Rank")
  ter_F<-aggregate(real_community_abundance[,c("Abundance")] ~ Family, data = real_community_abundance, sum) %>% mutate(Rank = "family")
  colnames(ter_F) <- c("Rank_name", "Abundance", "Rank")
  ter_G<-aggregate(real_community_abundance[,c("Abundance")] ~ Genus, data = real_community_abundance, sum) %>% mutate(Rank = "genus")
  colnames(ter_G) <- c("Rank_name", "Abundance", "Rank")
  ter_S<-aggregate(real_community_abundance[,c("Abundance")] ~ Species, data = real_community_abundance, sum) %>% mutate(Rank = "species")
  colnames(ter_S) <- c("Rank_name", "Abundance", "Rank")
  
  real_community_abundance_1 <- rbind(ter_K, ter_P, ter_C, ter_O, ter_F, ter_G, ter_S) %>% mutate(Mock_community=community)
  

  dataFrame = dataFrame %>% bind_rows(real_community_abundance_1)
  #write.csv(real_community_abundance, paste(paste("F:/PostDoc_dataset/bioinfomatic/performance evaluation/Processed_data/", community,sep = ""), ".csv", sep = ""))
}

write.csv(dataFrame,"F:/PostDoc_dataset/bioinfomatic/performance evaluation/Processed_data/Real_community_abundance.csv")


#calculate the predicted combined community abundance
dateFrame_community_length = data.frame()
for (length in c(70, 80, 90, 100, 110, 120, 130, 140, 150)) {
  print(length)
  dateFrame_community=data.frame()
  for (community in communitys){
    dateFrame = data.frame()
    for (db in c("ITS", "LSU", "ITS_LSU")) {
      for (rank in c("class", "order", "family", "genus", "species")) {
        for (filter in c("_filtered", "")) {
          file_name=paste(paste(paste(paste(paste(paste(paste(paste(paste("result_", community, sep = ""), ".short_read_min", sep = ""), length, sep=""), db, sep = "_")
                      , "_merged_output", sep = ""), filter, sep = ""), "_taxa_", sep = ""), rank, sep=""), ".tsv", sep = "")
          df = read_tsv(file_name)
          colnames(df) <- c("TaxonomyID", "Rank_name", 
                            #"Read_number", 
                            "Abundance")
          df_1=df %>% mutate(DB = db) %>% mutate(Rank = rank) %>% mutate(Filter = filter)
          
          dateFrame = dateFrame %>% dplyr::bind_rows(df_1)
          
        }
      }
    }
    dateFrame = dateFrame %>% dplyr::mutate(Mock_community = community)
    
    dateFrame_community = dateFrame_community %>% dplyr::bind_rows(dateFrame)
    
  }
  dateFrame_community=dateFrame_community %>% dplyr::mutate(MiniLength = length)
  
  dateFrame_community_length = dateFrame_community_length %>% dplyr::bind_rows(dateFrame_community)
}  
dateFrame_community_length$Mode = "raw"
write.csv(dateFrame_community_length,"F:/PostDoc_dataset/bioinfomatic/performance evaluation/Processed_data/Predicted_community_abundance_combine_mode_raw.csv")
dateFrame_community_length$Mode = "weighted"
write.csv(dateFrame_community_length,"F:/PostDoc_dataset/bioinfomatic/performance evaluation/Processed_data/Predicted_community_abundance_combine_mode_weighted.csv")
dateFrame_community_length$Mode = "weighted_length"
write.csv(dateFrame_community_length,"F:/PostDoc_dataset/bioinfomatic/performance evaluation/Processed_data/Predicted_community_abundance_combine_weighted_length.csv")

#calculate the predicted single database community abundance
dateFrame_community_length_cen = data.frame()
for (length in c(70, 80, 90, 100, 110, 120, 130, 140, 150)) {
  print(length)
  dateFrame_community=data.frame()
  for (community in communitys){
    dateFrame = data.frame()
    for (db in c("dbITS1_fisher", "dbITS2_fisher", "dbLSU_D1_fisher_new", "dbLSU_D2_fisher_new")) {
      filename = paste(paste(paste(paste(paste("result_", community, sep = ""), ".short_read_min", sep = ""), length, sep = ""), db, sep = "_"), "_report_kreport.tsv", sep = "")
      df = read_tsv(filename )
      df = df[,c(1, 4:6)]
      colnames(df) = c("Abundance" ,"rank", "TaxonomyID", "Rank_name")
      df_1 = df[which(df$rank != "-" & df$rank != "D" & df$rank != "K"),]
      df_1$Rank = ""
      df_1$Rank[which(df_1$rank == "P")] = "phylum"
      df_1$Rank[which(df_1$rank == "C")] = "class"
      df_1$Rank[which(df_1$rank == "O")] = "order"
      df_1$Rank[which(df_1$rank == "F")] = "family"
      df_1$Rank[which(df_1$rank == "G")] = "genus"
      df_1$Rank[which(df_1$rank == "S")] = "species"
      
      df_1$DB = db
      df_1$Mock_community = community
      df_1$MiniLength = length
      dateFrame_community_length_cen = dateFrame_community_length_cen %>% bind_rows(df_1)
      
    }
  }
}
write.csv(dateFrame_community_length_cen,"E:/PostDoc_dataset/bioinfomatic/performance evaluation/Processed_data/Predicted_community_abundance_single_DB.csv")




####################################################
#performance evluation 
#Section 2: evaluate the performance of Centrifuge and MicroFisher in taxa classification and abundance estimation
####################################################
rm(list=ls())
#setwd("F:/PostDoc_dataset/bioinfomatic/performance evaluation")
setwd("E:/PostDoc_dataset/bioinfomatic/performance evaluation/Processed_data")
list.files(getwd())

#calculate the evaluation performance of combined prediction
real_community = read.csv("Real_community_abundance.csv", row.names = 1)
#predicted_community_combine_raw = read.csv("Predicted_community_abundance_combine_mode_raw.csv") %>% dplyr::select(-Read_number)
predicted_community_combine_weighted = read.csv("Predicted_community_abundance_combine_mode_weighted.csv")
predicted_community_single_DB = read.csv("Predicted_community_abundance_single_DB.csv")
#predicted_community_combine_weighted_length = read.csv("Predicted_community_abundance_combine_weighted_length.csv")
predicted_community_combine = rbind( predicted_community_combine_weighted)

##identify the ture and false predictions for Microfisher output
##Histogram plotting of abundance of true and false prediction 
performance_evaluation = data.frame()
True_predicted = data.frame()
False_predicted = data.frame()
for (mode in c("weighted")) {
  for (length in c(70, 80, 90, 100, 110, 120, 130, 140, 150)) {
    print(length)
    dateFrame_community=data.frame()
    for (community in communitys){
      dateFrame = data.frame()
      for (db in c("ITS", "LSU", "ITS_LSU")) {
        for (rank in c("class", "order", "family", "genus", "species")) {
          for (filter in c("_filtered")) {
            predicted_df = predicted_community_combine %>% dplyr::filter(DB == db) %>% dplyr::filter(Rank == rank) %>% dplyr::filter(Filter == filter) %>%
              dplyr::filter(Mock_community == community) %>% dplyr::filter(MiniLength == length) %>% dplyr::filter(Mode == mode)
            real_df = real_community %>% dplyr::filter(Mock_community == community) %>% dplyr::filter(Rank == rank)  %>%  dplyr::distinct()
            
            df_evaluate = predicted_df %>% dplyr::inner_join(real_df, by = "Rank_name") %>% 
                          dplyr::rename(Predicted_abundance = Abundance.x, Real_abundance = Abundance.y, Mock_community = Mock_community.x, Rank = Rank.x)
            flase_evaluate = predicted_df %>% dplyr::full_join(real_df, by = "Rank_name") %>% 
                             dplyr::rename(Predicted_abundance = Abundance.x, Real_abundance = Abundance.y, Mock_community = Mock_community.x, Rank = Rank.x) %>%
                             dplyr::filter(is.na(Real_abundance))
            
            evaluate = data.frame(length = length, community = community, db = db, rank = rank, filter = filter, Mode = mode,
                                  Ture_Predicted = length(df_evaluate$Rank_name),
                                  Real_number = length(real_df$Rank_name),
                                  Predicted_Number = length(predicted_df$Rank_name))
            True_predicted = True_predicted %>% bind_rows(df_evaluate)
            False_predicted = False_predicted %>% bind_rows(flase_evaluate)
            performance_evaluation = performance_evaluation %>% bind_rows(evaluate)
          }
        }
      }
    }
  }
}
performance_evaluation = performance_evaluation %>% mutate(Sensitivity = performance_evaluation$Ture_Predicted/performance_evaluation$Real_number,
                                                           Precision = performance_evaluation$Ture_Predicted/performance_evaluation$Predicted_Number,
                                                           Accurate = performance_evaluation$Ture_Predicted/(performance_evaluation$Predicted_Number + performance_evaluation$Real_number - performance_evaluation$Ture_Predicted))
write.csv(True_predicted, "True_predicted_values_combine_mode_raw_weight.csv")
write.csv(False_predicted, "False_predicted_values_combine_mode_raw_weight.csv")
write.csv(performance_evaluation, "performance_evaluation_combine_mode_raw_weight.csv")


#calculate the evaluation performance of single database prediction 
##identify the ture and false predictions for single database output
##Histogram plotting of abundance of true and false prediction 
real_community = read.csv("Real_community_abundance.csv", row.names = 1)
predicted_community = read.csv("Predicted_community_abundance_single_DB.csv")
performance_evaluation = data.frame()
True_predicted = data.frame()
False_predicted = data.frame()
for (length in c(70, 80, 90, 100, 110, 120, 130, 140, 150)) {
  print(length)
  dateFrame_community=data.frame()
  for (community in communitys){
    dateFrame = data.frame()
    for (db in c("dbITS1_fisher", "dbITS2_fisher", "dbLSU_D1_fisher_new", "dbLSU_D2_fisher_new")) {
      for (rank_a in c("class", "order", "family", "genus", "species")) {
        for (filter in c("_filtered")) {
          predicted_df = predicted_community %>% dplyr::filter(DB == db ) %>%  dplyr::filter(Rank == rank_a) %>%
            dplyr::filter(Mock_community == community) %>% dplyr::filter(MiniLength == length)  
          real_df = real_community %>% dplyr::filter(Mock_community == community) %>% dplyr::filter(Rank == rank_a)  %>%  dplyr::distinct()
          
          df_evaluate = predicted_df %>% dplyr::inner_join(real_df, by = "Rank_name")
          
          df_evaluate = predicted_df %>% dplyr::inner_join(real_df, by = "Rank_name") %>% 
            dplyr::rename(Predicted_abundance = Abundance.x, Real_abundance = Abundance.y, Mock_community = Mock_community.x, Rank = Rank.x)
          flase_evaluate = predicted_df %>% dplyr::full_join(real_df, by = "Rank_name") %>% 
            dplyr::rename(Predicted_abundance = Abundance.x, Real_abundance = Abundance.y, Mock_community = Mock_community.x, Rank = Rank.x) %>%
            dplyr::filter(is.na(Real_abundance))
          
          evaluate = data.frame(length = length, community = community, db = db, rank = rank_a, filter = "single_DB_raw", 
                                Ture_Predicted = length(df_evaluate$Rank_name),
                                Real_number = length(real_df$Rank_name),
                                Predicted_Number = length(predicted_df$Rank_name))
          True_predicted = True_predicted %>% bind_rows(df_evaluate)
          False_predicted = False_predicted %>% bind_rows(flase_evaluate)
          performance_evaluation = performance_evaluation %>% bind_rows(evaluate)
        }
      }
    }
  }
}

performance_evaluation = performance_evaluation %>% mutate(Sensitivity = performance_evaluation$Ture_Predicted/performance_evaluation$Real_number,
                                                           Precision = performance_evaluation$Ture_Predicted/performance_evaluation$Predicted_Number,
                                                           Accurate = performance_evaluation$Ture_Predicted/(performance_evaluation$Predicted_Number + performance_evaluation$Real_number - performance_evaluation$Ture_Predicted))
write.csv(True_predicted, "True_predicted_values_singleDB.csv")
write.csv(False_predicted, "False_predicted_values_singleDB.csv")
write.csv(performance_evaluation, "performance_evaluation_single_DB.csv")



###################################################################################
# Section 3: plotting
#performance in taxa detection
performance_evaluation_combine = read.csv("performance_evaluation_combine_mode_raw_weight.csv")
performance_evaluation_single = read.csv("performance_evaluation_single_DB.csv") 
performance_evaluation = rbind(performance_evaluation_combine, mutate(performance_evaluation_single, Mode = "SingleDB"))
performance_evaluation_combine_new = data.frame()
#calculate the performance values
for (value in c("Sensitivity", "Precision", "Accurate")) {
  df = performance_evaluation %>% dplyr::select(1:10) %>% dplyr::mutate(Value = performance_evaluation[, value]) %>% mutate(Performance = value)
  performance_evaluation_combine_new = performance_evaluation_combine_new %>% bind_rows(df)
}
#generate the performance figures: minimum hit length
for (db_a in c("ITS", "LSU", "ITS_LSU", "dbITS1_fisher", "dbITS2_fisher", "dbLSU_D1_fisher_new", "dbLSU_D2_fisher_new")) {
      for (rank_a in c("class", "order", "family", "genus", "species")) {
        for (filter_a in c("_filtered")) {
          for (mode_a in c("raw", "weighted","SingleDB")) {
            df = performance_evaluation_combine_new %>% dplyr::filter( rank == rank_a & filter == filter_a & db == db_a & Mode == mode_a)
            plot_name = paste(paste(paste(paste("HitLength_", mode_a, sep = ""), db_a, sep = "_" ), rank_a, sep = "_"), filter_a, sep = "")
            plot <- ggplot(data = df, aes(x = length, y = Value, colour=Performance)) +
              geom_point(size = 4, shape = 16) +
              theme_bw() +
              expand_limits(y=c(0,1)) +
              stat_smooth(method ="loess", span = 1, se = TRUE, level = 0.95) +
              ggtitle(plot_name)+
              scale_y_continuous(minor_breaks = seq(0, 1, 0.5))+
              theme(legend.position = 'none',
                    axis.title.x = element_blank(),
                    axis.text.x  = element_blank(),# 字体的大小
                    axis.title.y = element_blank(),
                    axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                                vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                                size=10),
                    #axis.ticks.x=element_blank()
              )
              ggsave(paste(paste("figures/", plot_name, sep = ""), ".pdf" , sep = ""),height=3, width=4)
          }
        }
      }
    }
#generate the performance figures:species numbers       
for (db_a in c("ITS", "LSU", "ITS_LSU", "dbITS1_fisher", "dbITS2_fisher", "dbLSU_D1_fisher_new", "dbLSU_D2_fisher_new")) {
  for (value in c("Sensitivity", "Precision", "Accurate"))
    for (filter_a in c("_filtered")) {
      for (mode_a in c("weighted")) {
        df = performance_evaluation_combine_new  %>% mutate(speciesNum = community) %>% dplyr::filter( Performance == value & filter == filter_a & db == db_a, length == 120)
        df[,"speciesNum"] = str_replace( df[,"speciesNum"], "simulating_", "")
        df[,"speciesNum"] = str_replace( df[,"speciesNum"], "species_1", "")
        df[,"speciesNum"] = str_replace( df[,"speciesNum"], "species_2", "")
        df[,"speciesNum"] = str_replace( df[,"speciesNum"], "species_3", "")
        df[,"speciesNum"] = str_replace( df[,"speciesNum"], "species_4", "")
        df[,"speciesNum"] = str_replace( df[,"speciesNum"], "species_5", "")
        plot_name = paste(paste(paste(paste("speciesNum_", mode_a, sep = ""), db_a, sep = "_" ), value, sep = "_"), filter_a, sep = "")
        df$rank <- factor(df$rank, levels = c("species","genus","family","order","class"))
        plot <- ggplot(data = df, aes(x = rank, y = Value, group = speciesNum, color = speciesNum, shape = speciesNum)) +
          geom_point(size = 4) +
          theme_bw() +
          expand_limits(y=c(0,1)) +
          scale_shape_manual(values = c(0,1,2, 4,9,3,5)) +
          #stat_smooth(method ="loess", span = 1, se = TRUE, level = 0.95) +
          ggtitle(plot_name)+
          scale_y_continuous(minor_breaks = seq(0, 1, 0.5))+
          theme(legend.position = "none",
                axis.title.x = element_blank(),
                axis.text.x  = element_blank(),# 字体的大小
                axis.title.y = element_blank(),
                axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                            vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                            size=15),
                #axis.ticks.x=element_blank(),
                #axis.ticks.y=element_blank()
                )
         ggplot(data = df, aes(x = rank, y = Value,fill = speciesNum)) +
           geom_boxplot(alpha=1 ) +
           geom_jitter(width=0.25, alpha=0.5, color="black") +
           expand_limits(y=c(0,1)) +
           theme_bw() +
          #stat_smooth(method ="loess", span = 1, se = TRUE, level = 0.95) +
          ggtitle(plot_name)+
          scale_y_continuous(minor_breaks = seq(0, 1, 0.5))+
          theme(legend.position = "none",
                axis.title.x = element_blank(),
                axis.text.x  = element_blank(),# 字体的大小
                axis.title.y = element_blank(),
                axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                            vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                            size=15),
                #axis.ticks.x=element_blank(),
                #axis.ticks.y=element_blank()
          )
        ggsave(paste(paste("figures/", plot_name, sep = ""), ".pdf" , sep = ""),height=4, width=6)
      }
    }
  }
 
#comparison the relative abundance between ture and false prediction of MicroFisher output
True_predicted = read_csv("True_predicted_values_combine_mode_raw_weight.csv")
False_predicted = read_csv("False_predicted_values_combine_mode_raw_weight.csv")
for (rank_a in c("class", "order", "family", "genus", "species")){
  for (db in c("ITS", "LSU", "ITS_LSU")) { 
    df_true = True_predicted %>% dplyr::filter(Rank == rank_a & Filter == "_filtered" & MiniLength == 120 & Mode == "weighted" & DB == db) %>% dplyr::mutate(true_false = "true")
  df_false = False_predicted %>% dplyr::filter(Rank == rank_a & Filter == "_filtered" & MiniLength == 120 & Mode == "weighted" & DB == db) %>% dplyr::mutate(true_false = "flase")
  df_combine = df_true %>% dplyr::bind_rows(df_false)
  df = data.frame(value = df_combine$Predicted_abundance, true_false = df_combine$true_false)
  
  plot_name = paste(paste("density_plot", rank_a, sep = "_"), db, sep = "_")
  cols = c("#72D8FF", "#F76D5E")
  ggplot(df, aes(x=log10(value*100), fill = true_false)) + 
    geom_density(alpha = 0.7)+
    scale_fill_manual(values = cols)+
    #geom_violin(trim=FALSE)+
   labs(title="",x="Log10(Predict abundance (%))", y = "Density") +
    theme(legend.position = "right",
          axis.title.x = element_text(angle=0,# 设置旋转的角度
                                      vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                      size=10),
          axis.text.x  = element_text(angle=0,# 设置旋转的角度
                                      vjust=0.7,# 设置纵向廉价距离 hjust为横向偏移距离
                                      size=10,
                                      color = "black"),# 字体的大小
          axis.title.y = element_text(angle=90,# 设置旋转的角度
                                      vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                      size=10),
          axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                      vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                      size=10,
                                      color = "black")) 
  
  ggsave(paste(paste("figures/", plot_name, sep = ""), ".pdf" , sep = ""),height=2.5, width=3.5)
  
  
  df$value_2 = log10(df$value*100)
  ggplot2.histogram(data=df, xName='value_2',
                    groupName='true_false', 
                    position = "identity",
                    addMeanLine = F,
                    meanLineColor = "grey", 
                    meanLineType = "dashed", 
                    addDensityCurve = F, 
                    densityFill = "green",
                    densityAlpha = 0.2,
                    scale = c("frequency", "density"),
                    legendPosition="right") +
    expand_limits(y=c(0,180)) +
    labs(title="",x="Log10[Predict abundance (%)]", y = "Number of taxa")
  
  ggsave(paste(paste("figures/Histogram_of_", plot_name, sep = ""), "_2.pdf" , sep = ""),height=2.5, width=3.5)
  
  
  
  
  false = log10((df_false$Predicted_abundance)*100)
  ture = log10((df_true$Predicted_abundance)*100)
  
  #pdf(paste(paste("figures/Histogram_of_", plot_name, sep = ""), ".pdf" , sep = ""),height=2.5, width=3.5)
  # First distribution
  hist(ture, 
       breaks=1e2 ,
       #xlim=c(0,300), 
       col=rgb(1,0,0,0.5), 
       xlab="log10[Predicted abundance (%)]", 
       ylab="Number of taxa", 
       main= paste0(plot_name) )
  hist(false, 
       breaks=1e2, 
       #xlim=c(0,300), 
       col=rgb(0,0,1,0.5), 
       add=T)
  #dev.off()
  }
}

#comparison the relative abundance between ture and false prediction of single database output
True_predicted = read_csv("True_predicted_values_singleDB.csv")
False_predicted = read_csv("False_predicted_values_singleDB.csv")
for (rank_a in c("class", "order", "family", "genus", "species")){
  for (db in c("dbITS1_fisher" , "dbITS2_fisher" ,"dbLSU_D1_fisher_new", "dbLSU_D2_fisher_new")) { 
    df_true = True_predicted %>% dplyr::filter(Rank == rank_a & MiniLength == 120 & DB == db) %>% dplyr::mutate(true_false = "true")
    df_true$Real_abundance = df_true$Real_abundance*100
    df_false = False_predicted %>% dplyr::filter(Rank == rank_a  & MiniLength == 120  & DB == db) %>% dplyr::mutate(true_false = "flase")
    df_combine = df_true %>% dplyr::bind_rows(df_false)
    df = data.frame(value = df_combine$Predicted_abundance, true_false = df_combine$true_false)
    
    plot_name = paste(paste("density_plot", rank_a, sep = "_"), db, sep = "_")
    cols = c("#72D8FF", "#F76D5E")
    ggplot(df, aes(x=log10(value*100), fill = true_false)) + 
      geom_density(alpha = 0.7)+
      scale_fill_manual(values = cols)+
      #geom_violin(trim=FALSE)+
      labs(title="",x="Log10(Predict abundance (%))", y = "Density") +
      theme(legend.position = "right",
            axis.title.x = element_text(angle=0,# 设置旋转的角度
                                        vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                        size=10),
            axis.text.x  = element_text(angle=0,# 设置旋转的角度
                                        vjust=0.7,# 设置纵向廉价距离 hjust为横向偏移距离
                                        size=10,
                                        color = "black"),# 字体的大小
            axis.title.y = element_text(angle=90,# 设置旋转的角度
                                        vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                        size=10),
            axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                        vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                        size=10,
                                        color = "black")) 
    
    ggsave(paste(paste("figures/", plot_name, sep = ""), ".pdf" , sep = ""),height=2.5, width=3.5)
    
    
    df$value_2 = log10(df$value)
    ggplot2.histogram(data=df, xName='value_2',
                      groupName='true_false', 
                      position = "identity",
                      addMeanLine = F,
                      meanLineColor = "grey", 
                      meanLineType = "dashed", 
                      addDensityCurve = F, 
                      densityFill = "green",
                      densityAlpha = 0.2,
                      scale = c("frequency", "density"),
                      legendPosition="right") +
      labs(title="",x="Log10[Predict abundance (%)]", y = "Number of taxa")
    
    ggsave(paste(paste("figures/Histogram_of_", plot_name, sep = ""), ".pdf" , sep = ""),height=2.5, width=3.5)
    
    
    
    
    false = log10((df_false$Predicted_abundance)*100)
    ture = log10((df_true$Predicted_abundance)*100)
    
    #pdf(paste(paste("figures/Histogram_of_", plot_name, sep = ""), ".pdf" , sep = ""),height=2.5, width=3.5)
    # First distribution
    hist(ture, 
         breaks=1e2 ,
         #xlim=c(0,300), 
         col=rgb(1,0,0,0.5), 
         xlab="log10[Predicted abundance (%)]", 
         ylab="Number of taxa", 
         main= paste0(plot_name) )
    hist(false, 
         breaks=1e2, 
         #xlim=c(0,300), 
         col=rgb(0,0,1,0.5), 
         add=T)
    #dev.off()
  }
}

#Generate the histogram plot of actual taxa abundance
actual_abundance = read_csv("Real_community_abundance.csv")
for (rank_a in c("class", "order", "family", "genus", "species")){
  df = actual_abundance %>% dplyr::filter(Rank == rank_a)
  df = data.frame(value = df$Abundance, true_false = "True")
  df$value_2 = log10(df$value*100)
  plot_name = paste0("actual_abundance_", rank_a)
  ggplot2.histogram(data=df, xName='value_2',
                    groupName='true_false', 
                    groupColors = c("grey20"),
                    position = "identity",
                    addMeanLine = F,
                    meanLineColor = "grey", 
                    meanLineType = "dashed", 
                    addDensityCurve = F, 
                    densityFill = "green",
                    densityAlpha = 0.8,
                    scale = c("frequency", "density"),
                    legendPosition="right") +
    expand_limits(y=c(0,180)) +
    labs(title="",x="Log10[Predict abundance (%)]", y = "Number of taxa")
  ggsave(paste(paste("figures/Histogram_of_", plot_name, sep = ""), "_2.pdf" , sep = ""),height=2.5, width=3.5)
}




#Performance in abundance estimation
#########################################

#Error rank
#######################
predict_value=read.csv("True_predicted_values_combine_mode_raw_weight.csv")
True_predicted = read_csv("True_predicted_values_singleDB.csv")
True_predicted = True_predicted %>% mutate(Error = Predicted_abundance-Real_abundance)
#use MAE and RMSE for all the dataset at different taxonomy level
for (community in communitys) {
  for (rank_a in c("class", "order", "family", "genus", "species")) {
    for (filter_a in c("_filtered")) {
      for (mode_a in c("weighted")) {
        for (length_a in c(70,80,90,100,110,120,130,140,150)) {
            df = predict_value %>% dplyr::filter(Rank == rank_a & Filter == filter_a & MiniLength == length_a & Mode == mode_a & Mock_community == community)
            df_errors = data.frame()
            for (db_a in c("ITS", "LSU", "ITS_LSU")) {
              df_1 = df %>% filter(DB == db_a)
                df_error_0 = data.frame()
                for (i in unique(df_1$Rank_name)) {
                   df_tmp = df_1 %>% filter(Rank_name == i) %>% mutate(Error = (Predicted_abundance - Real_abundance)*100)
                   MAE = sum(abs(df_tmp$Error))/length(df_tmp$Error)
                   RMSE = sqrt(sum((df_tmp$Error)^2)/length(df_tmp$Error))
                   df_error = data.frame(Rank_Name = i, MAE = MAE, RMSE = RMSE, DB=db_a)
                   df_error_0 = df_error_0 %>% bind_rows(df_error)
                }
                df_error_0 = df_error_0 %>% dplyr::arrange(desc(RMSE)) %>% dplyr::mutate(RanK = 1:length(df_error_0$RMSE))
                df_errors = df_errors %>% bind_rows(df_error_0)
            }
            
            plot_name = paste(paste(paste(paste("MAE_", mode_a, sep = ""), length_a, sep = "_" ), rank_a, sep = "_"), filter_a, sep = "")
            ggplot(data = df_errors, aes(x = RanK, 
                                         y = RMSE, 
                                         colour=DB, 
                                         group = DB,
                                        # shape = DB
                                         )) +
              geom_line(aes(color = DB)) +
              geom_point(aes(color = DB) ) +
              scale_color_brewer(palette = 'Set1') +
              theme_bw() +
              ggtitle(plot_name)+
              expand_limits(y=c(0,max((df_errors$MAE))), x=c(0,length(unique(df_errors$RanK)))) +
              theme(legend.position = "right",
                    #axis.title.x = element_blank(),
                    #axis.text.x  = element_blank(),# 字体的大小
                    #axis.title.y = element_blank(),
                    axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                                vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                                size=10),
                    #axis.ticks.x=element_blank(),
                    #axis.ticks.y=element_blank()
              )         
            ggsave(paste(paste("figures/", plot_name, sep = ""), ".pdf" , sep = ""),height=3, width=4)
        }
      }
    }
  }
}
                    
#for seprate dataset at different taxonomy level
for (community in communitys) {
  for (rank_a in c("class", "order", "family", "genus", "species")) {
    for (filter_a in c("_filtered")) {
      for (mode_a in c( "weighted")) {
        #for (length_a in c(70,80,90,100,110,120,130,140,150)) {
        for (length_a in c(120)) {
          df = predict_value %>% dplyr::filter(Rank == rank_a & Filter == filter_a & MiniLength == length_a & Mode == mode_a & Mock_community == community)
          df_errors = data.frame()
          for (db_a in c("ITS", "LSU", "ITS_LSU")) {
            df_db = df %>% filter(DB == db_a) %>% mutate(Error =  abs((Predicted_abundance - Real_abundance)*100)) %>% 
              dplyr::select(Rank_name, Predicted_abundance, Filter, DB, Mock_community, Real_abundance, Error) %>%
              arrange(desc(Error)) 
            df_db = df_db %>% mutate(Error_rank = 1:length(df_db$Error))
            df_errors = df_errors %>% bind_rows(df_db)
          }
          
          RMSE_error = data.frame()
          for (db_a in c("ITS", "LSU", "ITS_LSU")) {
            df_tmp = df %>% filter(DB == db_a) %>% mutate(Error =  abs((Predicted_abundance - Real_abundance)*100))
            MAE = sum(abs(df_tmp$Error))/length(df_tmp$Error)
            RMSE = sqrt(sum((df_tmp$Error)^2)/length(df_tmp$Error))
            df_error = data.frame(DB = db_a, MAE = MAE, RMSE = RMSE)
            RMSE_error = RMSE_error %>% bind_rows(df_error)
          }
          ITS_RMSE = RMSE_error[1,3]
          LSU_RMSE = RMSE_error[2,3]
          ITS_LSU_RMSE = RMSE_error[3,3]
          write.table(RMSE_error, paste(paste("figures/", plot_name, sep = ""), ".txt" , sep = ""))
          
          
          plot_name = paste(paste(paste(paste(paste("Error_",community, sep = ""), mode_a, sep = "_"), length_a, sep = "_" ), rank_a, sep = "_"), filter_a, sep = "")
          p1=ggplot(data = df_errors, aes(x = Error_rank, 
                                       y = Error, 
                                       colour=DB, 
                                       group = DB,
                                       #shape = DB
                )) +
            geom_line(aes(color = DB)) +
            geom_point(aes(color = DB) ,shape = 1 ,size = 1) +
            scale_color_brewer(palette = 'Set1') +
            scale_shape_manual(values = 1)+
            scale_color_brewer(palette = 'Set1') +
            theme_bw() +
            ggtitle(paste(paste(paste("Taxonomy Level:", rank_a, sep = " "), community, sep = " ("), ")", sep = ""))+
            expand_limits(y=c(0,max((df_errors$Error))), x=c(0,length(unique(df_errors$RanK)))) +
            scale_color_discrete(name = "Marker gene",
                                 breaks = c("ITS", "LSU", "ITS_LSU"),
                                 labels = c(paste(paste("ITS1 & ITS2 (r.m.s. error:",round(ITS_RMSE, 2) , sep = " "), ")", sep = ""),
                                            paste(paste("LSU D1 & D2 (r.m.s. error:",round(LSU_RMSE, 2), sep = " "), ")", sep = ""),
                                            paste(paste("ITS & LSU (r.m.s. error:",round(ITS_LSU_RMSE, 2), sep = " "), ")", sep = ""))) +
            ylab("Absolute error in abundance estimation")+
            xlab(paste("Error ranked", rank_a , sep = " "))+
            theme(legend.position = c(0.67,  0.85),
                  axis.title.x = element_text(angle=0,# 设置旋转的角度
                                               vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                               size=15),
                  axis.text.x  = element_text(angle=0,# 设置旋转的角度
                                               vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                               size=15,
                                              color = "black"),# 字体的大小
                  axis.title.y = element_text(angle=90,# 设置旋转的角度
                                               vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                               size=15),
                  axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                              vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                              size=15,
                                              color = "black"),
                  axis.line = element_line(1,1:3),
                  #axis.ticks.x=element_blank(),
                  #axis.ticks.y=element_blank()
            )         
         # ggsave(paste(paste("figures/", plot_name, sep = ""), ".pdf" , sep = ""),height=6, width=7)
          
          p2= ggplot(df_errors, aes(x=DB, y=Error, fill=DB)) + 
            geom_violin(trim=FALSE)+
            geom_boxplot(width=0.1, fill="white")+
            labs(title="",x="", y = "Absolute error") +
            scale_fill_brewer(palette="Set1") +
            scale_x_discrete(breaks = c("ITS", "LSU", "ITS_LSU"),
                             labels = c("ITS1 & ITS2","LSU D1 & D2","ITS & LSU"))+
            theme(legend.position = "none",
                  axis.title.x = element_text(angle=30,# 设置旋转的角度
                                                                        vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                                                        size=10),
                  axis.text.x  = element_text(angle=30,# 设置旋转的角度
                                              vjust=0.7,# 设置纵向廉价距离 hjust为横向偏移距离
                                              size=10,
                                              color = "black"),# 字体的大小
                  axis.title.y = element_text(angle=90,# 设置旋转的角度
                                              vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                              size=10),
                  axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                              vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                              size=10,
                                              color = "black")) 
          
          p3 = qplot(DB, Error, data = df_errors, 
                    geom=c("violin"),
                     fill = DB ,colour = I("grey38")) + 
               xlab("")+
               ylab("Absolute error")+
                theme(legend.position = "none")
         # ggsave(paste(paste("figures/", plot_name, sep = ""), "_boxplot.pdf" , sep = ""),height=3, width=3.5)
          ggdraw() +
            draw_plot(p1, 0,0,1,1) 
          # draw_plot(p2, 0.5,0.22,0.4,0.47) 
          ggsave(paste(paste("figures/", plot_name, sep = ""), ".pdf" , sep = ""),height=5, width=5.5)
          

          
        }
      }
    }
  }
}
  

#for all the communities at different taxonomy level
for (rank_a in c("class", "order", "family", "genus", "species")) {
    for (filter_a in c("_filtered")) {
      for (mode_a in c( "weighted")) {
        #for (length_a in c(70,80,90,100,110,120,130,140,150)) {
        for (length_a in c(120)) {
          df = predict_value %>% dplyr::filter(Rank == rank_a & Filter == filter_a & MiniLength == length_a & Mode == mode_a )
          df_errors = data.frame()
          for (db_a in c("ITS", "LSU", "ITS_LSU")) {
            df_db = df %>% filter(DB == db_a) %>% mutate(Error =  abs((Predicted_abundance - Real_abundance)*100)) %>% 
              dplyr::select(Rank_name, Predicted_abundance, Filter, DB, Mock_community, Real_abundance, Error) %>%
              arrange(desc(Error)) 
            df_db = df_db %>% mutate(Error_rank = 1:length(df_db$Error))
            df_errors = df_errors %>% bind_rows(df_db)
          }
          
          plot_name = paste(paste(paste(paste(paste("Error_","all_community", sep = ""), mode_a, sep = "_"), length_a, sep = "_" ), rank_a, sep = "_"), filter_a, sep = "")
          
          RMSE_error = data.frame()
          for (db_a in c("ITS", "LSU", "ITS_LSU")) {
            df_tmp = df %>% filter(DB == db_a) %>% mutate(Error =  abs((Predicted_abundance - Real_abundance)*100))
            MAE = sum(abs(df_tmp$Error))/length(df_tmp$Error)
            RMSE = sqrt(sum((df_tmp$Error)^2)/length(df_tmp$Error))
            df_error = data.frame(DB = db_a, MAE = MAE, RMSE = RMSE)
            RMSE_error = RMSE_error %>% bind_rows(df_error)
          }
          ITS_RMSE = RMSE_error[1,3]
          LSU_RMSE = RMSE_error[2,3]
          ITS_LSU_RMSE = RMSE_error[3,3]
          write.table(RMSE_error, paste(paste("figures/", plot_name, sep = ""), ".txt" , sep = ""))
          
          

          p1=ggplot(data = df_errors, aes(x = Error_rank, 
                                          y = Error, 
                                          colour=DB, 
                                          group = DB,
                                          #shape = DB
                                          )) +
            geom_line(aes(color = DB)) +
            geom_point(aes(color = DB) ,shape = 1 ,size = 1) +
            scale_color_brewer(palette = 'Set1') +
            scale_shape_manual(values = 1)+
            scale_color_brewer(palette = 'Set1') +
            theme_bw() +
            ggtitle(paste(paste(paste("Taxonomy Level:", rank_a, sep = " "), "all_community", sep = " ("), ")", sep = ""))+
            expand_limits(y=c(0,max((df_errors$Error))), x=c(0,length(unique(df_errors$RanK)))) +
            scale_color_discrete(name = "Marker gene",
                                 breaks = c("ITS", "LSU", "ITS_LSU"),
                                 labels = c(paste(paste("ITS1 & ITS2 (r.m.s. error:",round(ITS_RMSE, 2) , sep = " "), ")", sep = ""),
                                            paste(paste("LSU D1 & D2 (r.m.s. error:",round(LSU_RMSE, 2), sep = " "), ")", sep = ""),
                                            paste(paste("ITS & LSU (r.m.s. error:",round(ITS_LSU_RMSE, 2), sep = " "), ")", sep = ""))) +
            ylab("Absolute error in abundance estimation")+
            xlab(paste("Error ranked", rank_a , sep = " "))+
            theme(legend.position = c(0.67,  0.85),
                  axis.title.x = element_text(angle=0,# 设置旋转的角度
                                              vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                              size=15),
                  axis.text.x  = element_text(angle=0,# 设置旋转的角度
                                              vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                              size=15,
                                              color = "black"),# 字体的大小
                  axis.title.y = element_text(angle=90,# 设置旋转的角度
                                              vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                              size=15),
                  axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                              vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                              size=15,
                                              color = "black"),
                  axis.line = element_line(1,1:3),
                  #axis.ticks.x=element_blank(),
                  #axis.ticks.y=element_blank()
            )         
          # ggsave(paste(paste("figures/", plot_name, sep = ""), ".pdf" , sep = ""),height=6, width=7)
          
          p2= ggplot(df_errors, aes(x=DB, y=Error, fill=DB)) + 
            geom_violin(trim=FALSE)+
            geom_boxplot(width=0.2, fill="white")+
            labs(title="",x="", y = "Absolute error") +
            scale_fill_brewer(palette="Set1") +
            scale_x_discrete(breaks = c("ITS", "LSU", "ITS_LSU"),
                             labels = c("ITS1 & ITS2","LSU D1 & D2","ITS & LSU"))+
            theme(legend.position = "none",
                  axis.title.x = element_text(angle=30,# 设置旋转的角度
                                              vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                              size=10),
                  axis.text.x  = element_text(angle=30,# 设置旋转的角度
                                              vjust=0.7,# 设置纵向廉价距离 hjust为横向偏移距离
                                              size=10,
                                              color = "black"),# 字体的大小
                  axis.title.y = element_text(angle=90,# 设置旋转的角度
                                              vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                              size=10),
                  axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                              vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                              size=10,
                                              color = "black")) 
          
          p3 = qplot(DB, Error, data = df_errors, 
                     geom=c("violin"),
                     fill = DB ,colour = I("grey38")) + 
            xlab("")+
            ylab("Absolute error")+
            theme(legend.position = "none")
          # ggsave(paste(paste("figures/", plot_name, sep = ""), "_boxplot.pdf" , sep = ""),height=3, width=3.5)
          ggdraw() +
            draw_plot(p1, 0,0,1,1) 
            #draw_plot(p2, 0.5,0.22,0.4,0.47) 
          ggsave(paste(paste("figures/", plot_name, sep = ""), ".pdf" , sep = ""),height=5, width=5.5)
          
          
          
        }
      }
    }
  }

#for all the communities at different database level
for (db_a in c("ITS", "LSU", "ITS_LSU")) {
  for (filter_a in c("_filtered")) {
    for (mode_a in c( "weighted")) {
      #for (length_a in c(70,80,90,100,110,120,130,140,150)) {
      for (length_a in c(120)) {
        df = predict_value %>% dplyr::filter(DB == db_a & Filter == filter_a & MiniLength == length_a & Mode == mode_a )
        df_errors = data.frame()
        for (rank_a in c("class", "order", "family", "genus", "species"))  {
          df_rank = df %>% filter(Rank == rank_a) %>% mutate(Error =  abs((Predicted_abundance - Real_abundance)*100)) %>% 
            dplyr::select(Rank_name, Rank, Predicted_abundance, Filter, DB, Mock_community, Real_abundance, Error) %>%
            arrange(desc(Error)) 
          df_rank = df_rank %>% mutate(Error_rank = 1:length(df_rank$Error))
          df_rank = df_rank %>% mutate(Error_rank_norm = Error_rank/length(df_rank$Error))
          df_errors = df_errors %>% bind_rows(df_rank)
        }
        
        plot_name = paste(paste(paste(paste(paste("Error_","all_community", sep = ""), mode_a, sep = "_"), length_a, sep = "_" ), db_a, sep = "_"), filter_a, sep = "")
        
        RMSE_error = data.frame()
        for (rank_a in c("class", "order", "family", "genus", "species")) {
          df_tmp = df %>% filter(Rank == rank_a) %>% mutate(Error =  abs((Predicted_abundance - Real_abundance)*100))
          MAE = sum(abs(df_tmp$Error))/length(df_tmp$Error)
          RMSE = sqrt(sum((df_tmp$Error)^2)/length(df_tmp$Error))
          df_error = data.frame(Rank = rank_a, MAE = MAE, RMSE = RMSE)
          RMSE_error = RMSE_error %>% bind_rows(df_error)
        }
        class_RMSE = RMSE_error[1,3]
        len_class = length(which(df_errors$Rank == "class") == "True")
        order_RMSE = RMSE_error[2,3]
        len_order = length(which(df_errors$Rank == "order") == "True")
        family_RMSE = RMSE_error[3,3]
        len_family = length(which(df_errors$Rank == "family") == "True")
        genus_RMSE = RMSE_error[4,3]
        len_genus = length(which(df_errors$Rank == "genus") == "True")
        species_RMSE = RMSE_error[5,3]
        len_species = length(which(df_errors$Rank == "species") == "True")
        write.table(RMSE_error, paste(paste("figures/", plot_name, sep = ""), ".txt" , sep = ""))
        
        
        
        p1=ggplot(data = df_errors, aes(x = Error_rank_norm, 
                                        y = Error, 
                                        colour=Rank, 
                                        group = Rank,
                                        #shape = DB
                                        )) +
          geom_line(aes(color = Rank)) +
          geom_point(aes(color = Rank) ,shape = 1 ,size = 0.5) +
          scale_color_brewer(palette = 'Set1') +
          scale_shape_manual(values = 1)+
          scale_color_brewer(palette = 'Set1') +
          theme_bw() +
          ggtitle(paste(paste(paste("Database:", db_a, sep = " "), "all_community", sep = " ("), ")", sep = ""))+
          expand_limits(y=c(0,max((df_errors$Error))), x=c(0,length(unique(df_errors$RanK)))) +
          scale_color_discrete(name = "Taxonomy levels",
                               breaks = c("class", "order", "family", "genus", "species"),
                               labels = c(paste(paste(paste("Class (r.m.s. error:",round(class_RMSE, 2) , sep = " "), ")", sep = ""), len_class, sep = " ("),
                                          paste(paste(paste("Order (r.m.s. error:",round(order_RMSE, 2), sep = " "), ")", sep = ""), len_order, sep = " ("),
                                          paste(paste(paste("Family (r.m.s. error:",round(family_RMSE, 2), sep = " "), ")", sep = ""), len_family, sep = " ("),
                                          paste(paste(paste("Genus (r.m.s. error:",round(genus_RMSE, 2), sep = " "), ")", sep = ""), len_genus, sep = " ("),
                                          paste(paste(paste("Species (r.m.s. error:",round(species_RMSE, 2), sep = " "), ")", sep = ""), len_species, sep = " (")
                                          )) +
          ylab("Absolute error in abundance estimation")+
          xlab(paste("Error rank", "" , sep = ""))+
          theme(legend.position = c(0.67,  0.55),
                axis.title.x = element_text(angle=0,# 设置旋转的角度
                                            vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                            size=15),
                axis.text.x  = element_text(angle=0,# 设置旋转的角度
                                            vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                            size=15,
                                            color = "black"),# 字体的大小
                axis.title.y = element_text(angle=90,# 设置旋转的角度
                                            vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                            size=15),
                axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                            vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                            size=15,
                                            color = "black"),
                axis.line = element_line(1,1:3),
                #axis.ticks.x=element_blank(),
                #axis.ticks.y=element_blank()
          )         
        # ggsave(paste(paste("figures/", plot_name, sep = ""), ".pdf" , sep = ""),height=6, width=7)
        
        p2= ggplot(df_errors, aes(x=DB, y=Error, fill=DB)) + 
          geom_violin(trim=FALSE)+
          geom_boxplot(width=0.2, fill="white")+
          labs(title="",x="", y = "Absolute error") +
          scale_fill_brewer(palette="Set1") +
          scale_x_discrete(breaks = c("ITS", "LSU", "ITS_LSU"),
                           labels = c("ITS1 & ITS2","LSU D1 & D2","ITS & LSU"))+
          theme(legend.position = "none",
                axis.title.x = element_text(angle=30,# 设置旋转的角度
                                            vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                            size=10),
                axis.text.x  = element_text(angle=30,# 设置旋转的角度
                                            vjust=0.7,# 设置纵向廉价距离 hjust为横向偏移距离
                                            size=10,
                                            color = "black"),# 字体的大小
                axis.title.y = element_text(angle=90,# 设置旋转的角度
                                            vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                            size=10),
                axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                            vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                            size=10,
                                            color = "black")) 
        
        p3 = qplot(DB, Error, data = df_errors, 
                   geom=c("violin"),
                   fill = DB ,colour = I("grey38")) + 
          xlab("")+
          ylab("Absolute error")+
          theme(legend.position = "none")
        # ggsave(paste(paste("figures/", plot_name, sep = ""), "_boxplot.pdf" , sep = ""),height=3, width=3.5)
        ggdraw() +
          draw_plot(p1, 0,0,1,1) 
        #draw_plot(p2, 0.5,0.22,0.4,0.47) 
        ggsave(paste(paste("figures/", plot_name, sep = ""), ".pdf" , sep = ""),height=3, width=5.5)
        
        
        
      }
    }
  }
}

#for all the communities at different classification tools
predict_value = predict_value %>% dplyr::filter(Filter == "_filtered", Mode == "weighted")
singleDB = data.frame(Rank = True_predicted$Rank, DB = True_predicted$DB, Mock_community = True_predicted$Mock_community, 
                      MiniLength = True_predicted$MiniLength, TaxonomyID = True_predicted$TaxonomyID, Rank_name = True_predicted$Rank_name,
                      Real_abundance = True_predicted$Real_abundance, 
                      Predicted_abundance = (True_predicted$Predicted_abundance)/100)

combineDB = data.frame(Rank = predict_value$Rank, DB = predict_value$DB, Mock_community = predict_value$Mock_community, 
                       MiniLength = predict_value$MiniLength, TaxonomyID = predict_value$TaxonomyID, Rank_name = predict_value$Rank_name,
                       Real_abundance = predict_value$Real_abundance, 
                       Predicted_abundance = predict_value$Predicted_abundance)
all_abundance_data = singleDB %>% dplyr::bind_rows(combineDB)


DBs = c("dbITS1_fisher", "dbITS2_fisher", "dbLSU_D1_fisher_new", "dbLSU_D2_fisher_new", "ITS", "LSU", "ITS_LSU")

    for (rank_a in c("class", "order", "family", "genus", "species")){
      #for (length_a in c(70,80,90,100,110,120,130,140,150)) {
      for (length_a in c(120)) {
        df = all_abundance_data %>% dplyr::filter( MiniLength == length_a, Rank == rank_a)
        df_errors = data.frame()
        for (db_a in c("dbITS1_fisher", "dbITS2_fisher", "dbLSU_D1_fisher_new", "dbLSU_D2_fisher_new", "ITS", "LSU", "ITS_LSU"))  {
          df_rank = df %>% filter(DB == db_a) %>% mutate(Error =  abs((Predicted_abundance - Real_abundance)*100)) %>% 
            arrange(desc(Error)) 
          df_rank = df_rank %>% mutate(Error_rank = 1:length(df_rank$Error))
          df_rank = df_rank %>% mutate(Error_rank_norm = Error_rank/length(df_rank$Error))
          df_errors = df_errors %>% bind_rows(df_rank)
        }
        
        plot_name = paste(paste(paste("ErrorRank_","all_community", sep = ""), length_a, sep = "_" ), rank_a, sep = "")
        
        RMSE_error = data.frame()
        for (db_a in c("dbITS1_fisher", "dbITS2_fisher", "dbLSU_D1_fisher_new", "dbLSU_D2_fisher_new", "ITS", "LSU", "ITS_LSU")) {
          df_tmp = df %>% filter(DB == db_a) %>% mutate(Error =  abs((Predicted_abundance - Real_abundance)*100))
          MAE = sum(abs(df_tmp$Error))/length(df_tmp$Error)
          RMSE = sqrt(sum((df_tmp$Error)^2)/length(df_tmp$Error))
          df_error = data.frame(DB = db_a, MAE = MAE, RMSE = RMSE)
          RMSE_error = RMSE_error %>% bind_rows(df_error)
        }
        ITS1_RMSE = RMSE_error[1,3]
        len_ITS1 = length(which(df_errors$DB == "dbITS1_fisher") == "True")
        ITS2_RMSE = RMSE_error[2,3]
        len_ITS2 = length(which(df_errors$DB == "dbITS2_fisher") == "True")
        D1_RMSE = RMSE_error[3,3]
        len_D1 = length(which(df_errors$DB == "dbLSU_D1_fisher_new") == "True")
        D2_RMSE = RMSE_error[4,3]
        len_D2 = length(which(df_errors$DB == "dbLSU_D2_fisher_new") == "True")
        ITS_RMSE = RMSE_error[5,3]
        len_ITS = length(which(df_errors$DB == "ITS") == "True")
        LSU_RMSE = RMSE_error[6,3]
        len_LSU = length(which(df_errors$DB == "LSU") == "True")
        ITS_LSU_RMSE = RMSE_error[7,3]
        len_ITS_LSU = length(which(df_errors$DB == "ITS_LSU") == "True")
        write.table(RMSE_error, paste(paste("figures/", plot_name, sep = ""), ".txt" , sep = ""))
        
        
        
    p1 = ggplot(data = df_errors, aes(x = Error_rank_norm, 
                                        y = Error, 
                                        colour=DB, 
                                        group = DB,
                                        #shape = DB
        )) +
          geom_line(aes(color = DB)) +
          geom_point(aes(color = DB) ,shape = 1 ,size = 0.5) +
          scale_color_brewer(palette = 'Set2') +
          scale_shape_manual(values = 1)+
          xlim(0, 0.1)+
          theme_bw() +
          ggtitle(paste(paste(paste("Taxonomic level:", rank_a, sep = " "), "All_community", sep = " ("), ")", sep = ""))+
          expand_limits(y=c(0,max((df_errors$Error))), x=c(0,length(unique(df_errors$RanK)))) +
          scale_color_discrete(name = "HMDs",
                               breaks = c("dbITS1_fisher", "dbITS2_fisher", "dbLSU_D1_fisher_new", "dbLSU_D2_fisher_new", "ITS", "LSU", "ITS_LSU"),
                               labels = c(paste0("ITS1 (r.m.s. error:",round(ITS1_RMSE, 2), "; ", "Total taxa:", len_ITS1, ")"),
                                          paste0("ITS2 (r.m.s. error:",round(ITS2_RMSE, 2), "; ", "Total taxa:", len_ITS2, ")"),
                                          paste0("LSU D1 (r.m.s. error:",round(D1_RMSE, 2), "; ", "Total taxa:", len_D1, ")"),
                                          paste0("LSU D2 (r.m.s. error:",round(D2_RMSE, 2), "; ", "Total taxa:", len_D2, ")"),
                                          paste0("ITS (r.m.s. error:",round(ITS_RMSE, 2), "; ", "Total taxa:", len_ITS, ")"),
                                          paste0("LSU (r.m.s. error:",round(LSU_RMSE, 2), "; ", "Total taxa:", len_LSU, ")"),
                                          paste0("ITS_LSU (r.m.s. error:",round(ITS_LSU_RMSE, 2), "; ", "Total taxa:", len_ITS_LSU, ")")
                               )) +
          ylab("Absolute error in abundance estimation (%)")+
          xlab(paste("Fraction of estimated taxa", "" , sep = ""))+
          theme(legend.position = c(0.7,  0.65),
                axis.title.x = element_text(angle=0,# 设置旋转的角度
                                            vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                            size=15),
                axis.text.x  = element_text(angle=0,# 设置旋转的角度
                                            vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                            size=15,
                                            color = "black"),# 字体的大小
                axis.title.y = element_text(angle=90,# 设置旋转的角度
                                            vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                            size=15),
                axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                            vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                            size=15,
                                            color = "black"),
                axis.line = element_line(1,1:3),
                #axis.ticks.x=element_blank(),
                #axis.ticks.y=element_blank()
          )         
        # ggsave(paste(paste("figures/", plot_name, sep = ""), ".pdf" , sep = ""),height=6, width=7)
        
        p2= ggplot(df_errors, aes(x=DB, y=Error, fill=DB)) + 
          geom_violin(trim=FALSE)+
          geom_boxplot(width=0.2, fill="white")+
          labs(title="",x="", y = "Absolute error") +
          scale_fill_brewer(palette="Set1") +
          scale_x_discrete(breaks = c("ITS", "LSU", "ITS_LSU"),
                           labels = c("ITS1 & ITS2","LSU D1 & D2","ITS & LSU"))+
          theme(legend.position = "none",
                axis.title.x = element_text(angle=30,# 设置旋转的角度
                                            vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                            size=10),
                axis.text.x  = element_text(angle=30,# 设置旋转的角度
                                            vjust=0.7,# 设置纵向廉价距离 hjust为横向偏移距离
                                            size=10,
                                            color = "black"),# 字体的大小
                axis.title.y = element_text(angle=90,# 设置旋转的角度
                                            vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                            size=10),
                axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                            vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                            size=10,
                                            color = "black")) 
        
        p3 = qplot(DB, Error, data = df_errors, 
                   geom=c("violin"),
                   fill = DB ,colour = I("grey38")) + 
          xlab("")+
          ylab("Absolute error")+
          theme(legend.position = "none")
        # ggsave(paste(paste("figures/", plot_name, sep = ""), "_boxplot.pdf" , sep = ""),height=3, width=3.5)
        ggdraw() +
          draw_plot(p1, 0,0,1,1) 
        #draw_plot(p2, 0.5,0.22,0.4,0.47) 
        ggsave(paste(paste("figures/", plot_name, sep = ""), ".pdf" , sep = ""),height=4, width=7)
        
        
        
      }
    }






#correlations of Predicted- and True- abundance of taxa 
#generate the correlation plot for all communities
for (db in c("ITS", "LSU", "ITS_LSU")) {            
predict_value_selected = predict_value %>% filter(Filter == "_filtered", DB == db, MiniLength == 120, Mode == "weighted") %>% 
                         dplyr::mutate(MiniLength = replace(MiniLength, MiniLength ==80, "80")) %>%
                         dplyr::mutate(MiniLength = replace(MiniLength, MiniLength ==120, "120")) %>%
                         dplyr::filter( Rank == "genus" |Rank == ""| Rank == "class")
predict_value_selected_2 = predict_value_selected %>% 
  filter(Real_abundance <= 0.05, Predicted_abundance <= 0.05) 

plot_name = paste(paste("correlations_", db, sep = ""), "_minLen_rank", sep = "")


corr_class = cor.test(predict_value_selected$Real_abundance[which(predict_value_selected$Rank == "class")],
                    predict_value_selected$Predicted_abundance[which(predict_value_selected$Rank == "class")])
corr_genus = cor.test(predict_value_selected$Real_abundance[which(predict_value_selected$Rank == "genus")],
                    predict_value_selected$Predicted_abundance[which(predict_value_selected$Rank == "genus")])
corr_class_log = cor.test(log10(predict_value_selected_2$Real_abundance[which(predict_value_selected_2$Rank == "class")]),
                          log10(predict_value_selected_2$Predicted_abundance[which(predict_value_selected_2$Rank == "class")]))
corr_genus_log = cor.test(log10(predict_value_selected_2$Real_abundance[which(predict_value_selected_2$Rank == "genus")]),
                        log10(predict_value_selected_2$Predicted_abundance[which(predict_value_selected_2$Rank == "genus")]))

corr_table = data.frame(Database=c("class", "genus","class_log", "genus_log"),
                        info=c(corr_class$data.name, corr_genus$data.name,corr_class_log$data.name, corr_genus_log$data.name),
                        corr_value=c(corr_class$estimate, corr_genus$estimate, corr_class_log$estimate, corr_genus_log$estimate),
                        pvalue=c(corr_class$p.value,corr_genus$p.value,corr_class_log$p.value, corr_genus_log$p.value))


write.table(corr_table,paste(paste("figures/", plot_name, sep = ""), "corr.txt" , sep = "") )

corr_value_class=corr_class$estimate
corr_value_genus=corr_genus$estimate
corr_value_class_log=corr_class_log$estimate
corr_value_genus_log=corr_genus_log$estimate

corr_Pvalue_class=corr_class$p.value
corr_Pvalue_genus=corr_genus$p.value
corr_Pvalue_class_log=corr_class_log$p.value
corr_Pvalue_genus_log=corr_genus_log$p.value


p1=ggplot(data = predict_value_selected, aes(x = Real_abundance*100, y = Predicted_abundance*100, group = Rank,shape = Rank, colour=Rank)) +
  geom_point(aes(color = Rank, fill = Rank), size = 2.5, alpha = 0.7) +
  theme_bw() +
  #stat_smooth(method ="loess", span = 2, se = TRUE, level = 0.95) +
  #scale_color_brewer(palette = 'Set1') +
  scale_colour_manual(values = c("forestgreen","firebrick1","forestgreen",  "darkorchid")) +
  scale_shape_manual(values = c(16,4, 17, 4,9,3,5)) +
  #expand_limits(y=c(0,5), x=c(0,5)) +
  ##expand_limits(y=c(0,max((df_errors$RMSE))), x=c(0,length(unique(df_errors$RanK)))) +
  scale_colour_manual(name = "Rank",values = c("forestgreen","firebrick1", "forestgreen", "darkorchid"),
                      breaks = c("class", "genus"),
                      labels = c(paste(paste(paste(paste("Class (corr=",round(corr_value_class, 2) , sep = ""), "p-value:", sep = ";"),corr_Pvalue_class, sep = ""),")",sep=""),
                                 paste(paste(paste(paste("Genus (corr=",round(corr_value_genus, 2) , sep = ""), "p-value:", sep = ";"),corr_Pvalue_genus, sep = ""),")",sep=""))) +

  ggtitle(plot_name)+
  theme(legend.position = c(0.25, 0.75),
        axis.title.x = element_text(angle=0,# 设置旋转的角度
                                    vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                    size=15),
        axis.text.x  = element_text(angle=0,# 设置旋转的角度
                                    vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                    size=15,
                                    color = "black"),# 字体的大小
        axis.title.y = element_text(angle=90,# 设置旋转的角度
                                    vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                    size=15),
        axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                    vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                    size=15,
                                    color = "black"),
        axis.line = element_line(1,1:3),
        #axis.ticks.x=element_blank(),
        #axis.ticks.y=element_blank()
  )+
  xlab("Expected abundance") +
  ylab("Predicted abundance")


 




 p2 = ggplot(data = predict_value_selected_2, aes(x = log10(Real_abundance*100), y = log10(Predicted_abundance*100), colour=Rank, shape = Rank, group = Rank)) +
  geom_point(aes(color = Rank, fill = Rank), size = 1.5, alpha = 1) +
  theme_bw() +
  #stat_smooth(method ="loess", span = 2, se = TRUE, level = 0.95) +
  #scale_color_brewer(palette = 'Set1') +
  scale_colour_manual(name = "Rank",values = c("forestgreen","firebrick1","forestgreen",  "brown"),
                      breaks = c("class","genus"),
                      labels = c(paste(paste(paste(paste("Class (corr=",round(corr_value_class_log, 2) , sep = ""), "p-value:", sep = ";"),corr_Pvalue_class_log, sep = ""),")",sep=""),
                                 paste(paste(paste(paste("Genus (corr=",round(corr_value_genus_log, 2) , sep = ""), "p-value:", sep = ";"),corr_Pvalue_genus_log, sep = ""),")",sep=""))) +

  scale_shape_manual(values = c(16,4)) +
  #expand_limits(y=c(0,5), x=c(0,5)) +
  ##expand_limits(y=c(0,max((df_errors$RMSE))), x=c(0,length(unique(df_errors$RanK)))) +
  #ggtitle(plot_name)+
  scale_x_continuous(position = "top")+
  scale_y_continuous(position = "left")+
  xlab("Log10(Expected abundance)") +
  ylab("Log10(Predicted abundance)")+
  theme(legend.position = c(0.25,0.75),
        #axis.title.x = element_blank(),
        # axis.text.x  = element_blank(),# 字体的大小
        # axis.title.y = element_blank(),
        axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                    vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                    size=10),
        axis.line = element_line(1,1:3),
        #axis.ticks.x=element_blank(),
        #axis.ticks.y=element_blank()
  )   

p2 = p2 + geom_point(data = predict_value_selected_2[which(predict_value_selected_2$Rank == "class"),], 
                aes(color = Rank, fill = Rank), size = 1.5, alpha = 1)
  

predict_value_selected_2 = predict_value_selected %>% 
  filter(Real_abundance <= 0.05, Predicted_abundance <= 0.05) %>% filter(Rank == "species")
p3 = ggplot(data = predict_value_selected_2, aes(x = log10(Real_abundance*100), y = log10(Predicted_abundance*100), colour=MiniLength, group = MiniLength)) +
  geom_point(aes(color = MiniLength, fill = MiniLength), size = 0.7) +
  theme_bw() +
  #stat_smooth(method ="loess", span = 2, se = TRUE, level = 0.95) +
  #scale_color_brewer(palette = 'Set1') +
  scale_colour_manual(values = c( "darkblue","red2")) +
  scale_shape_manual(values = c(21,25, 4,9,3,5)) +
  #expand_limits(y=c(0,5), x=c(0,5)) +
  ##expand_limits(y=c(0,max((df_errors$RMSE))), x=c(0,length(unique(df_errors$RanK)))) +
  #ggtitle(plot_name)+
  scale_x_continuous(position = "bottom")+
  scale_y_continuous(position = "right")+
  xlab("Log10(Expected abundance)") +
  ylab("Log10(predicted abundance)")+
  theme(legend.position = c(0.25, 0.75),
        #axis.title.x = element_blank(),
        # axis.text.x  = element_blank(),# 字体的大小
        # axis.title.y = element_blank(),
        axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                    vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                    size=10),
        axis.line = element_line(1,1:3),
        #axis.ticks.x=element_blank(),
        #axis.ticks.y=element_blank()
  )   
#ggsave(paste(paste("figures/", plot_name, sep = ""), "species_scale_lim.pdf" , sep = ""),height=3, width=4)

ggdraw() +
  draw_plot(p1, 0,0,1,1) +
  draw_plot(p2, 0.6,0.08,0.4,0.4)
#+  draw_plot(p3, 0.08,0.6,0.35,0.35)
ggsave( paste(paste("figures/", plot_name, sep = ""), "_scale_lim.pdf" , sep = ""),height=7, width=7)
}

#generate the correlation plot for separate community
for (rank_a in c("class", "order", "family", "genus", "species")) {      
  for (mock_community in communitys) {
  predict_value_selected = predict_value %>% filter(Filter == "_filtered", Rank == rank_a, MiniLength == 120, Mode == "weighted", Mock_community == mock_community) %>% 
    dplyr::mutate(MiniLength = replace(MiniLength, MiniLength ==120, "Min_120")) 

  plot_name = paste(paste(paste("correlations_", rank_a, sep = ""), "minLen_rank", sep = "_"), mock_community, sep = "_")
  
  
  
  corr_overall = cor.test(log10(predict_value_selected$Real_abundance), log10(predict_value_selected$Predicted_abundance))
  corr_ITS = cor.test(log10(predict_value_selected$Real_abundance[which(predict_value_selected$DB == "ITS")]),
                      log10(predict_value_selected$Predicted_abundance[which(predict_value_selected$DB == "ITS")]))
  corr_LSU = cor.test(log10(predict_value_selected$Real_abundance[which(predict_value_selected$DB == "LSU")]),
                      log10(predict_value_selected$Predicted_abundance[which(predict_value_selected$DB == "LSU")]))
  corr_ITS_LSU = cor.test(log10(predict_value_selected$Real_abundance[which(predict_value_selected$DB == "ITS_LSU")]),
                          log10(predict_value_selected$Predicted_abundance[which(predict_value_selected$DB == "ITS_LSU")]))
  corr_table = data.frame(Database=c("overall", "ITS", "LSU", "ITS_LSU"),
                          info=c(corr_overall$data.name, corr_ITS$data.name, corr_LSU$data.name, corr_ITS_LSU$data.name),
                          corr_value=c(corr_overall$estimate, corr_value=corr_ITS$estimate, corr_value=corr_LSU$estimate, corr_value=corr_ITS_LSU$estimate),
                          pvalue=c(corr_overall$p.value,corr_ITS$p.value, corr_LSU$p.value, corr_ITS_LSU$p.value))
  
  
  write.table(corr_table,paste(paste("figures/", plot_name, sep = ""), "corr.txt" , sep = "") )
  
  corr_value_ITS=corr_ITS$estimate
  corr_value_LSU=corr_LSU$estimate
  corr_value_ITS_LSU=corr_ITS_LSU$estimate
  
  corr_Pvalue_ITS=corr_ITS$p.value
  corr_Pvalue_LSU=corr_LSU$p.value
  corr_Pvalue_ITS_LSU=corr_ITS_LSU$p.value
  
  
  
  
  ggplot(data = predict_value_selected, aes(x = log10(Real_abundance*100), y = log10(Predicted_abundance*100), group = DB,shape = DB, colour=DB)) +
    geom_point(aes(color = DB, fill = DB), size = 1.5) +
    theme_bw() +
    #stat_smooth(method ="loess", span = 2, se = TRUE, level = 0.95) +
    #scale_color_brewer(palette = 'Set1') +
    scale_colour_manual(values = c("red1", "darkorchid", "forestgreen")) +
    scale_shape_manual(values = c(16,17, 4,9,3,5)) +
    #expand_limits(y=c(0,5), x=c(0,5)) +
    ##expand_limits(y=c(0,max((df_errors$RMSE))), x=c(0,length(unique(df_errors$RanK)))) +
    scale_colour_manual(name = "Marker gene",values = c("red1", "darkorchid", "forestgreen"),
                         breaks = c("ITS", "LSU", "ITS_LSU"),
                         labels = c(paste(paste(paste(paste("ITS1 & ITS2 (corr=",round(corr_value_ITS, 2) , sep = ""), "p-value:", sep = ";"),corr_Pvalue_ITS, sep = ""),")",sep=""),
                                    paste(paste(paste(paste("LSU D1 & D2 (corr=",round(corr_value_LSU, 2) , sep = ""), "p-value:", sep = ";"),corr_Pvalue_LSU, sep = ""),")",sep=""),
                                    paste(paste(paste(paste("ITS & LSU (corr=",round(corr_value_ITS_LSU, 2) , sep = ""), "p-value:", sep = ";"),corr_Pvalue_ITS_LSU, sep = ""),")",sep=""))) +
    scale_shape_manual(name = "Marker gene",values = c(16,17, 4,9,3,5),
                       breaks = c("ITS", "LSU", "ITS_LSU"),
                       labels = c(paste(paste(paste(paste("ITS1 & ITS2 (corr=",round(corr_value_ITS, 2) , sep = ""), "p-value:", sep = ";"),corr_Pvalue_ITS, sep = ""),")",sep=""),
                                  paste(paste(paste(paste("LSU D1 & D2 (corr=",round(corr_value_LSU, 2) , sep = ""), "p-value:", sep = ";"),corr_Pvalue_LSU, sep = ""),")",sep=""),
                                  paste(paste(paste(paste("ITS & LSU (corr=",round(corr_value_ITS_LSU, 2) , sep = ""), "p-value:", sep = ";"),corr_Pvalue_ITS_LSU, sep = ""),")",sep=""))) +
    ggtitle(plot_name)+
    theme(legend.position = c(0.25, 0.75),
          axis.title.x = element_text(angle=0,# 设置旋转的角度
                                      vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                      size=15),
          axis.text.x  = element_text(angle=0,# 设置旋转的角度
                                      vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                      size=15,
                                      color = "black"),# 字体的大小
          axis.title.y = element_text(angle=90,# 设置旋转的角度
                                      vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                      size=15),
          axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                      vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                      size=15,
                                      color = "black"),
          axis.line = element_line(1,1:3),
          #axis.ticks.x=element_blank(),
          #axis.ticks.y=element_blank()
    )+
    xlab("log10(Expected abundance)") +
    ylab("log10(Predicted abundance)")
  
  ggsave(paste(paste("figures/", plot_name, sep = ""), ".pdf" , sep = ""),height=5, width=5)
  

  }
}
  
for (rank_a in c("class", "order", "family", "genus", "species")) {      
    predict_value_selected = predict_value %>% filter(Filter == "_filtered", Rank == rank_a, MiniLength == 120, Mode == "weighted") %>% 
      dplyr::mutate(MiniLength = replace(MiniLength, MiniLength ==120, "Min_120")) 
    
    plot_name = paste(paste("correlations_", rank_a, sep = ""), "minLen_rank", sep = "_")
    
    
    
    corr_overall = cor.test(log10(predict_value_selected$Real_abundance), log10(predict_value_selected$Predicted_abundance))
    corr_ITS = cor.test(log10(predict_value_selected$Real_abundance[which(predict_value_selected$DB == "ITS")]),
                        log10(predict_value_selected$Predicted_abundance[which(predict_value_selected$DB == "ITS")]))
    corr_LSU = cor.test(log10(predict_value_selected$Real_abundance[which(predict_value_selected$DB == "LSU")]),
                        log10(predict_value_selected$Predicted_abundance[which(predict_value_selected$DB == "LSU")]))
    corr_ITS_LSU = cor.test(log10(predict_value_selected$Real_abundance[which(predict_value_selected$DB == "ITS_LSU")]),
                            log10(predict_value_selected$Predicted_abundance[which(predict_value_selected$DB == "ITS_LSU")]))
    corr_table = data.frame(Database=c("overall", "ITS", "LSU", "ITS_LSU"),
                            info=c(corr_overall$data.name, corr_ITS$data.name, corr_LSU$data.name, corr_ITS_LSU$data.name),
                            corr_value=c(corr_overall$estimate, corr_value=corr_ITS$estimate, corr_value=corr_LSU$estimate, corr_value=corr_ITS_LSU$estimate),
                            pvalue=c(corr_overall$p.value,corr_ITS$p.value, corr_LSU$p.value, corr_ITS_LSU$p.value))
    
    
    write.table(corr_table,paste(paste("figures/", plot_name, sep = ""), "corr.txt" , sep = "") )
    
    corr_value_ITS=corr_ITS$estimate
    corr_value_LSU=corr_LSU$estimate
    corr_value_ITS_LSU=corr_ITS_LSU$estimate
    
    corr_Pvalue_ITS=corr_ITS$p.value
    corr_Pvalue_LSU=corr_LSU$p.value
    corr_Pvalue_ITS_LSU=corr_ITS_LSU$p.value
    
    
    
    
    ggplot(data = predict_value_selected, aes(x = log10(Real_abundance*100), y = log10(Predicted_abundance*100), group = DB,shape = DB, colour=DB)) +
      geom_point(aes(color = DB, fill = DB), size = 1.5) +
      theme_bw() +
      #stat_smooth(method ="loess", span = 2, se = TRUE, level = 0.95) +
      #scale_color_brewer(palette = 'Set1') +
      scale_colour_manual(values = c("red1", "darkorchid", "forestgreen")) +
      scale_shape_manual(values = c(16,17, 4,9,3,5)) +
      #expand_limits(y=c(0,5), x=c(0,5)) +
      ##expand_limits(y=c(0,max((df_errors$RMSE))), x=c(0,length(unique(df_errors$RanK)))) +
      scale_colour_manual(name = "Marker gene",values = c("red1", "darkorchid", "forestgreen"),
                         breaks = c("ITS", "LSU", "ITS_LSU"),
                         labels = c(paste(paste(paste(paste("ITS1 & ITS2 (corr=",round(corr_value_ITS, 2) , sep = ""), "p-value:", sep = ";"),corr_Pvalue_ITS, sep = ""),")",sep=""),
                                    paste(paste(paste(paste("LSU D1 & D2 (corr=",round(corr_value_LSU, 2) , sep = ""), "p-value:", sep = ";"),corr_Pvalue_LSU, sep = ""),")",sep=""),
                                    paste(paste(paste(paste("ITS & LSU (corr=",round(corr_value_ITS_LSU, 2) , sep = ""), "p-value:", sep = ";"),corr_Pvalue_ITS_LSU, sep = ""),")",sep=""))) +
      scale_shape_manual(name = "Marker gene",values = c(16,17, 4,9,3,5),
                           breaks = c("ITS", "LSU", "ITS_LSU"),
                           labels = c(paste(paste(paste(paste("ITS1 & ITS2 (corr=",round(corr_value_ITS, 2) , sep = ""), "p-value:", sep = ";"),corr_Pvalue_ITS, sep = ""),")",sep=""),
                                      paste(paste(paste(paste("LSU D1 & D2 (corr=",round(corr_value_LSU, 2) , sep = ""), "p-value:", sep = ";"),corr_Pvalue_LSU, sep = ""),")",sep=""),
                                      paste(paste(paste(paste("ITS & LSU (corr=",round(corr_value_ITS_LSU, 2) , sep = ""), "p-value:", sep = ";"),corr_Pvalue_ITS_LSU, sep = ""),")",sep=""))) +
      ggtitle(plot_name)+
      theme(legend.position = c(0.45, 0.85),
            axis.title.x = element_text(angle=0,# 设置旋转的角度
                                        vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                        size=15),
            axis.text.x  = element_text(angle=0,# 设置旋转的角度
                                        vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                        size=15,
                                        color = "black"),# 字体的大小
            axis.title.y = element_text(angle=90,# 设置旋转的角度
                                        vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                        size=15),
            axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                        vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                        size=15,
                                        color = "black"),
            axis.line = element_line(1,1:3),
            #axis.ticks.x=element_blank(),
            #axis.ticks.y=element_blank()
      )+
      xlab("log10(Expected abundance)") +
      ylab("log10(Predicted abundance)")
    
    ggsave(paste(paste("figures/", plot_name, sep = ""), ".pdf" , sep = ""),height=5, width=5)
    
    
  }























































  
  predict_value_selected_2 = predict_value_selected %>% 
    filter(Real_abundance <= 0.05, Predicted_abundance <= 0.05) %>% filter(Rank == "genus")
  p2 = ggplot(data = predict_value_selected_2, aes(x = log10(Real_abundance*100), y = log10(Predicted_abundance*100), colour=MiniLength, group = MiniLength)) +
    geom_point(aes(color = MiniLength, fill = MiniLength), size = 0.7) +
    theme_bw() +
    #stat_smooth(method ="loess", span = 2, se = TRUE, level = 0.95) +
    #scale_color_brewer(palette = 'Set1') +
    scale_colour_manual(values = c( "steelblue", "limegreen")) +
    scale_shape_manual(values = c(21,25, 4,9,3,5)) +
    #expand_limits(y=c(0,5), x=c(0,5)) +
    ##expand_limits(y=c(0,max((df_errors$RMSE))), x=c(0,length(unique(df_errors$RanK)))) +
    #ggtitle(plot_name)+
    scale_x_continuous(position = "top")+
    scale_y_continuous(position = "left")+
    theme(legend.position = c(0.2,0.8),
          #axis.title.x = element_blank(),
          # axis.text.x  = element_blank(),# 字体的大小
          # axis.title.y = element_blank(),
          axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                      vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                      size=10),
          axis.line = element_line(1,1:3),
          #axis.ticks.x=element_blank(),
          #axis.ticks.y=element_blank()
    )   
  
  #ggsave(paste(paste("figures/", plot_name, sep = ""), "genus_scale_lim.pdf" , sep = ""),height=3, width=4)
  
  predict_value_selected_2 = predict_value_selected %>% 
    filter(Real_abundance <= 0.05, Predicted_abundance <= 0.05) %>% filter(Rank == "species")
  p3 = ggplot(data = predict_value_selected_2, aes(x = log10(Real_abundance*100), y = log10(Predicted_abundance*100), colour=MiniLength, group = MiniLength)) +
    geom_point(aes(color = MiniLength, fill = MiniLength), size = 0.7) +
    theme_bw() +
    #stat_smooth(method ="loess", span = 2, se = TRUE, level = 0.95) +
    #scale_color_brewer(palette = 'Set1') +
    scale_colour_manual(values = c( "steelblue", "limegreen")) +
    scale_shape_manual(values = c(21,25, 4,9,3,5)) +
    #expand_limits(y=c(0,5), x=c(0,5)) +
    ##expand_limits(y=c(0,max((df_errors$RMSE))), x=c(0,length(unique(df_errors$RanK)))) +
    #ggtitle(plot_name)+
    scale_x_continuous(position = "bottom")+
    scale_y_continuous(position = "right")+
    theme(legend.position = c(0.2, 0.8),
          #axis.title.x = element_blank(),
          # axis.text.x  = element_blank(),# 字体的大小
          # axis.title.y = element_blank(),
          axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                      vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                      size=10),
          axis.line = element_line(1,1:3),
          #axis.ticks.x=element_blank(),
          #axis.ticks.y=element_blank()
    )   
  #ggsave(paste(paste("figures/", plot_name, sep = ""), "species_scale_lim.pdf" , sep = ""),height=3, width=4)
  
  ggdraw() +
    draw_plot(p1, 0,0,1,1) +
    draw_plot(p2, 0.5,0.05,0.35,0.4) +
    draw_plot(p3, 0.05,0.5,0.35,0.4)
  ggsave( paste(paste("figures/", plot_name, sep = ""), "_scale_lim.pdf" , sep = ""),height=10, width=12)
}

























db="ITS"
community="simulating_200species_2"
mode_a = "raw"
rank_a="class"
filter_a="_filtered"
db_a = "ITS"
length_a=120
length=120

# setwd("F:/PostDoc_dataset/bioinfomatic/performance evaluation")
# ##########
# #Read the data 1
# 
# combine <- read.table("hitlength_DBs_kreport_combine.txt",header = T)
# separate <- read.table("hitlength_DBs_separated.txt",header = T)
# df <- rbind(separate, combine)
# colnames(df)
# df$Sensitivity <- df$Ture_Positive/(df$Ture_Positive + df$False_negative)
# df$Accuracy <- df$Ture_Positive/(df$Ture_Positive + df$False_negative + df$False_Positive)
# df$Precision <- df$Ture_Positive/(df$Ture_Positive  + df$False_Positive)
# 
# 
# Sensitivity <- df[,c("group_id","speciesNum","Database","length","taxa_level", "Sensitivity")]
# colnames(Sensitivity) <- c("group_id","speciesNum","Database","length","taxa_level","value"  )
# Sensitivity$Performance <- "Sensitivity"
# Accuracy <- df[,c("group_id","speciesNum","Database","length","taxa_level","Accuracy"  )]
# colnames(Accuracy) <- c("group_id","speciesNum","Database","length","taxa_level","value"  )
# Accuracy$Performance <- "Accuracy"
# Precision <- df[,c("group_id","speciesNum","Database","length","taxa_level","Precision"  )]
# colnames(Precision) <- c("group_id","speciesNum","Database","length","taxa_level","value"  )
# Precision$Performance <- "Precision"
# 
# final_table <- rbind(Sensitivity, Accuracy, Precision)


#combine real data classification
rm(list=ls())
setwd("F:/PostDoc_dataset/bioinfomatic/existing_comparison/final_results")
file_list = read.table("eukdetect/file.txt")
#combine classification result from eukdetect

euk_result_combine_count = data.frame(Rank = "", taxaName = "")
euk_result_combine_rank = data.frame(Rank = "", taxaName = "")
for (file in file_list$V1) {
  fileName = paste(file, "_filtered_hits_taxonomy.txt", sep = "")
  print(fileName)
  euk_result = read_tsv(paste("eukdetect/", fileName, sep = ""))
  euk_result_readCount = data.frame(Rank = euk_result$Rank, 
                                    taxaName = euk_result$Name,
                                    file = euk_result$Marker_read_count)
  colnames(euk_result_readCount) = c("Rank", "taxaName", file)
  euk_result_readCount_count = euk_result_readCount[, 2:3]
  euk_result_readCount_rank = euk_result_readCount[, 2:1]
  euk_result_combine_count = euk_result_combine_count %>%
                             dplyr::full_join(euk_result_readCount_count, by = "taxaName")
  euk_result_combine_rank = euk_result_combine_rank %>%
                            dplyr::bind_rows(euk_result_readCount_rank)                
}
euk_result_combine_rank = unique(euk_result_combine_rank)
euk_result_combine = dplyr::full_join(euk_result_combine_count, euk_result_combine_rank,  by = "taxaName")
write.table(euk_result_combine, "eukdetect/combine.txt")
write.csv(euk_result_combine, "eukdetect/combine.csv")
#combine classification result from kaiju
file_list = read.table("file.txt")
kaiju_result_table_final = data.frame()
levels = c("class", "family", "genus", "order", "phylum", "species")
for (taxa_level in levels) {
  taxonomy_file = paste(paste("kaiju/kaiju_summary_", taxa_level, sep = ""), ".table", sep = "")
  kaiju_result_table = data.frame(taxon_name="")
  for (file in file_list$V1) {
    fileName = paste(file, "_kaiju.out", sep = "")
    abund_tab = read_tsv(taxonomy_file)
    
    abundance_table = abund_tab %>%
      dplyr::filter(file == fileName) %>%
      dplyr::select(taxon_name, reads)
    colnames(abundance_table) = c("taxon_name", file)
    
    kaiju_result_table = kaiju_result_table %>%
      dplyr::full_join(abundance_table, by = "taxon_name")
  }
  kaiju_result_table$tax_level = taxa_level
  kaiju_result_table_final = kaiju_result_table_final %>%
                            dplyr::bind_rows(kaiju_result_table)
}

write.table(kaiju_result_table_final, "kaiju/kaiju_summary_all_final.txt")

  
#combine classification result from microfisher
file_list = read.table("microfisher/file.txt")


levels = c("class", "family", "genus", "order", "species")
    result_combine_4 = data.frame()
    for (type in c("DNA", "RNA")) {
       result_combine_3 = data.frame()
       for (length in c(80, 120, 150)) {
         result_combine_2 = data.frame()
         for (taxon_level in levels) {
            result_combine = data.frame(name = "", proportion = "")
            for (file in file_list$V1) {
              folder = paste(paste("microfisher/merged_results_meta", paste(type, length, sep = "_"), sep = ""), file, sep = "_")
              fileName = paste(paste("merged_output_filtered_taxa_", taxon_level, sep = ""), ".tsv", sep = "")
              
              result = read_tsv(paste(folder, fileName, sep = "/"))
              result_2 = result %>% dplyr::select(name, proportion)
              colnames(result_2) <- c("name", file)
              result_combine = result_combine %>%
                              dplyr::full_join(result_2, by = "name")
            }
            result_combine$taxa_level = taxon_level
            result_combine$MiniLength = length
            result_combine$seq_type <- type
            result_combine_2 = rbind(result_combine_2, result_combine)
         }
         
         result_combine_3 = rbind(result_combine_3, result_combine_2)
       }
       
       result_combine_4 = rbind(result_combine_4, result_combine_3)
    }

write.table(result_combine_4, "microfisher/result_combine_4.txt")

check = result_combine_4[which(result_combine_4$taxa_level == "genus" & result_combine_4$seq_type == "DNA"),]




#read the classification result
df=read.table("microfisher/result_combine_4.txt")
#df_mg=df[,c(1, 3:57,173:175)] %>% dplyr::filter(MiniLength == 120 & seq_type == "DNA")
df_mg=df[,c(1, 58:172,173:175)] %>% dplyr::filter(MiniLength == 120 & seq_type == "RNA")
df_mg[is.na(df_mg)] = 0 
df_mg$all = rowMeans(df_mg[2:56])
df_mg = df_mg %>% dplyr::filter(all > 0.001)
df_mg$all
df_mg = df_mg %>% dplyr::filter(taxa_level == "genus")











#read the database taxa info
rm(list=ls())
setwd("F:/PostDoc_dataset/bioinfomatic/1_database")

its1 = read_tsv("ITS1_lineage.txt") 
colnames(its1) = c("taxid", "superkingdom","kingdom","phylum","class", "order", "family", "genus", "species")
its1 = its1 %>% dplyr::mutate(db = "ITS1")
its2 = read_tsv("ITS2_lineage.txt")
colnames(its2) = c("taxid", "superkingdom","kingdom","phylum","class", "order", "family", "genus", "species")
its2 = its2 %>% dplyr::mutate(db = "ITS2")
lsuD1 = read_tsv("lsuD1_lineage.txt")
colnames(lsuD1) = c("taxid", "superkingdom","kingdom","phylum","class", "order", "family", "genus", "species")
lsuD1 = lsuD1 %>% dplyr::mutate(db = "lsuD1")
lsuD2 = read_tsv("lsuD2_lineage.txt") 
colnames(lsuD2) = c("taxid", "superkingdom","kingdom","phylum","class", "order", "family", "genus", "species")
lsuD2 = lsuD2 %>% dplyr::mutate(db = "lsuD2")
df = bind_rows(its1, its2,lsuD1, lsuD2 )
df = df %>% filter(kingdom == "Fungi")

colnames(df)
df_select = df %>% filter(db != "")
print(length(unique(df_select$phylum)))
print(length(unique(df_select$class)))
print(length(unique(df_select$order)))
print(length(unique(df_select$family)))
print(length(unique(df_select$genus)))
print(length(unique(df_select$species)))





summarise(ranks)




ranks = c("superkingdom","kingdom","phylum","class", "order", "family", "genus", "species")
for (rank in ranks) {
  df_rank = df[,c(rank, "db")]  %>% dplyr::filter(df[rank] != "", df[rank] != "NA") %>% dplyr::distinct()

  write.csv(df_rank, paste(rank,"combine_lineage_for_venn.csv", sep = "_"))
}

rank="species"



df_rank = data.frame(df_rank)

venn.plot = venn.diagram(list(df_rank[which(df_rank$db == "ITS1"),][,1], 
                              df_rank[which(df_rank$db == "ITS2"),][,1], 
                              df_rank[which(df_rank$db == "lsuD1"),][,1], 
                              df_rank[which(df_rank$db == "lsuD2"),][,1]), 
                        NULL,fill=c("red", "green", "blue", "grey"), 
                        alpha=c(0.5,0.5,0.5,0.5), 
                        cex = 2,
                        cat.fontface=4, category.names=c("ITS1", "ITS2","LSU D1", "LSU D2"))
grid.draw(venn.plot)
dev.off()


#https://www.r-bloggers.com/2019/04/set-analysis-a-face-off-between-venn-diagrams-and-upset-plots/







df =  read.table("F:/PostDoc_dataset/bioinfomatic/existing_comparison/final_results/result_combine_4.txt")
metadata = read_tsv("F:/PostDoc_dataset/bioinfomatic/existing_comparison/final_results/metadata.txt")


df_1 = df %>% filter(MiniLength == 120) %>% filter(seq_type == "DNA") 
df_1(is.na(df_1)) = 0 
df_1[is.na(df_1)] = 0





df_1[which(df_1 == "NA")] = 0
df_2 = df_1[,c(metadata$File)]
df_2$File = row.names(df_2)

df_3 = df_2 %>% dplyr::left_join(metadata, by = "File")
write.csv(df_1, "F:/PostDoc_dataset/bioinfomatic/existing_comparison/final_results/metaRNA.csv")
a














































#############################################
#generate the plot

Database = list("ITS1","ITS2","LsuD1","LsuD2","ITS","LSU","ITS_LSU")
DBs = c("ITS1", "ITS2","LsuD1","LsuD2","ITS", "LSU","ITS_LSU")

for (database in DBs) {

species <- final_table[which(final_table$taxa_level == "Species" & final_table$Database == database),]
S <- ggplot(data = species, aes(x = length, y = value, colour=Performance)) +
  geom_point() +
  theme_bw() +
  expand_limits(y=c(0,1)) +
  stat_smooth(method ="loess", span = 1, se = TRUE, level = 0.95) +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),# 字体的大小
        axis.title.y = element_blank(),
        axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                    vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                    size=10),
        #axis.ticks.x=element_blank()
  )
genus <- final_table[which(final_table$taxa_level == "Genus" & final_table$Database == database),]
G <- ggplot(data = genus, aes(x = length, y = value, colour=Performance)) +
  geom_point() +
  theme_bw() +
  expand_limits(y=c(0,1)) +
  stat_smooth(method ="loess", span = 1, se = TRUE, level = 0.95) +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),# 字体的大小
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        #axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())
family <- final_table[which(final_table$taxa_level == "Family" & final_table$Database == database),]
F <- ggplot(data = family, aes(x = length, y = value, colour=Performance)) +
  geom_point() +
  theme_bw() +
  expand_limits(y=c(0,1)) +
  stat_smooth(method ="loess", span = 1, se = TRUE, level = 0.95) +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),# 字体的大小
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        #axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())
order <- final_table[which(final_table$taxa_level == "Order" & final_table$Database == database),]
O <- ggplot(data = order, aes(x = length, y = value, colour=Performance)) +
  geom_point() +
  theme_bw() +
  expand_limits(y=c(0,1)) +
  stat_smooth(method ="loess", span = 1, se = TRUE, level = 0.95) +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),# 字体的大小
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        #axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())
class <- final_table[which(final_table$taxa_level == "Class" & final_table$Database == database),]
C <- ggplot(data = class, aes(x = length, y = value, colour=Performance)) +
  geom_point() +
  theme_bw() +
  expand_limits(y=c(0,1)) +
  stat_smooth(method ="loess", span = 1, se = TRUE, level = 0.95) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),# 字体的大小
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        #axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())
phylum <- final_table[which(final_table$taxa_level == "Phylum" & final_table$Database == database),]
P <- ggplot(data = phylum, aes(x = length, y = value, colour=Performance)) +
  geom_point() +
  theme_bw() +
  expand_limits(y=c(0,1)) +
  stat_smooth(method ="loess", span = 1, se = TRUE, level = 0.95) +
  theme(legend.position = c(0.8,0.2) ,
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        #axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())

assign(paste("HitLength_",database,sep = ""),
  ggdraw() +
  draw_plot(S, 0,0,0.22,1) +
  draw_plot(G, 0.22,0,0.19,1) +
  draw_plot(F, 0.41,0,0.19,1) +
  draw_plot(O, 0.60,0,0.19,1) +
  draw_plot(C, 0.79,0,0.19,1)
)

}




pdf("plots/hitlength_DBs_2.pdf", height=8, width=8)
gridExtra::grid.arrange(HitLength_ITS1,HitLength_ITS2,HitLength_LsuD1,HitLength_LsuD2,HitLength_ITS,HitLength_LSU,HitLength_ITS_LSU,nrow=7)
dev.off()







########################################################################
#line plot
##########################################################################
library(dplyr)
library(ggplot2)
library(cowplot)

rm(list=ls())
setwd("F:/PostDoc_dataset/bioinfomatic/performance evaluation")
list.files(getwd())

##########
#Read the data 1

combine <- read.table("speciesNum_DBs_kreport_combine.txt",header = T)
separate <- read.table("speciesNum_DBs_separated.txt",header = T)
df <- rbind(separate, combine)

colnames(df)
df$Sensitivity <- df$Ture_Positive/(df$Ture_Positive + df$False_negative)
df$Accuracy <- df$Ture_Positive/(df$Ture_Positive + df$False_negative + df$False_Positive)
df$Precision <- df$Ture_Positive/(df$Ture_Positive  + df$False_Positive)
head(df)

Sensitivity <- df[,c("group_id","speciesNum","Database","length","taxa_level", "Sensitivity")]
colnames(Sensitivity) <- c("group_id","speciesNum","Database","length","taxa_level","value"  )
Sensitivity$Performance <- "Sensitivity"
Accuracy <- df[,c("group_id","speciesNum","Database","length","taxa_level","Accuracy"  )]
colnames(Accuracy) <- c("group_id","speciesNum","Database","length","taxa_level","value"  )
Accuracy$Performance <- "Accuracy"
Precision <- df[,c("group_id","speciesNum","Database","length","taxa_level","Precision"  )]
colnames(Precision) <- c("group_id","speciesNum","Database","length","taxa_level","value"  )
Precision$Performance <- "Precision"

final_table <- rbind(Sensitivity, Accuracy, Precision)
colnames(final_table)



#############################################
#generate the plot


DBs = c("ITS1", "ITS2","LsuD1","LsuD2","ITS", "LSU","ITS_LSU")
for (database in DBs){
  print(database)

Sensitivity <- final_table[which(final_table$Performance == "Sensitivity" & final_table$Database == database & final_table$length == 130),]

meandata = group_by(Sensitivity, speciesNum, taxa_level) %>%
             summarise(value=mean(value))
meandata$taxa_level <- factor(meandata$taxa_level, levels = c("Species","Genus","Family","Order","Class","Phylum"))
Sensivity_plot <- ggplot(meandata[which(meandata$taxa_level != "Phylum"),], aes(x = taxa_level, y = value, group = speciesNum,shape = speciesNum)) +
  geom_line(aes(color = speciesNum)) +
  geom_point(aes(color = speciesNum)) +
  theme_bw() +
  expand_limits(y=c(0,1)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),# 字体的大小
        axis.title.y = element_blank(),
        axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                    vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                    size=10),
        #axis.ticks.x=element_blank(),
        #axis.ticks.y=element_blank()
        )
  
  
  

Accuracy <- final_table[which(final_table$Performance == "Accuracy" & final_table$Database == database & final_table$length == 130),]

meandata = group_by(Accuracy, speciesNum, taxa_level) %>%
  summarise(value=mean(value))
meandata$taxa_level <- factor(meandata$taxa_level, levels = c("Species","Genus","Family","Order","Class","Phylum"))
Accuracy_plot <- ggplot(meandata[which(meandata$taxa_level != "Phylum"),], aes(x = taxa_level, y = value, group = speciesNum,shape = speciesNum)) +
  geom_line(aes(color = speciesNum)) +
  geom_point(aes(color = speciesNum)) +
  theme_bw() +
  expand_limits(y=c(0,1)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),# 字体的大小
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        #axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()
  )



Precision <- final_table[which(final_table$Performance == "Precision" & final_table$Database == database & final_table$length == 130),]

meandata = group_by(Precision, speciesNum, taxa_level) %>%
  summarise(value=mean(value))
meandata$taxa_level <- factor(meandata$taxa_level, levels = c("Species","Genus","Family","Order","Class","Phylum"))
Precision_plot <- ggplot(meandata[which(meandata$taxa_level != "Phylum"),], aes(x = taxa_level, y = value, group = speciesNum,shape = speciesNum)) +
  geom_line(aes(color = speciesNum)) +
  geom_point(aes(color = speciesNum)) +
  theme_bw() +
  expand_limits(y=c(0,1)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),# 字体的大小
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        #axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()
  )


assign(paste("speciesNum_",database,sep = ""),
       ggdraw() +
       draw_plot(Sensivity_plot, 0,0,0.34,1) +
       draw_plot(Accuracy_plot, 0.34,0,0.32,1) +
       draw_plot(Precision_plot, 0.66,0,0.32,1)
      )
}


pdf("plots/SpeciesNum_DBs_2.pdf", height=7, width=7)
gridExtra::grid.arrange(speciesNum_ITS1,speciesNum_ITS2,speciesNum_LsuD1,speciesNum_LsuD2,speciesNum_ITS,speciesNum_LSU,speciesNum_ITS_LSU,nrow=7)
dev.off()



 








################################################################################################

#read the database taxa info
rm(list=ls())
setwd("/home/microbiome/data_storage/temperal/R")

its1 = read_tsv("ITS1_lineage.txt") 
colnames(its1) = c("taxid", "superkingdom","kingdom","phylum","class", "order", "family", "genus", "species")
its1 = its1 %>% dplyr::mutate(db = "ITS1")
its2 = read_tsv("ITS2_lineage.txt")
colnames(its2) = c("taxid", "superkingdom","kingdom","phylum","class", "order", "family", "genus", "species")
its2 = its2 %>% dplyr::mutate(db = "ITS2")
lsuD1 = read_tsv("lsuD1_lineage.txt")
colnames(lsuD1) = c("taxid", "superkingdom","kingdom","phylum","class", "order", "family", "genus", "species")
lsuD1 = lsuD1 %>% dplyr::mutate(db = "lsuD1")
lsuD2 = read_tsv("lsuD2_lineage.txt") 
colnames(lsuD2) = c("taxid", "superkingdom","kingdom","phylum","class", "order", "family", "genus", "species")
lsuD2 = lsuD2 %>% dplyr::mutate(db = "lsuD2")
df = bind_rows(its1, its2,lsuD1, lsuD2 )
df = df %>% filter(kingdom == "Fungi")

colnames(df)
df_select = df %>% filter(db != "")
print(length(unique(df_select$phylum)))
print(length(unique(df_select$class)))
print(length(unique(df_select$order)))
print(length(unique(df_select$family)))
print(length(unique(df_select$genus)))
print(length(unique(df_select$species)))





summarise(ranks)




ranks = c("superkingdom","kingdom","phylum","class", "order", "family", "genus", "species")
for (rank in ranks) {
  df_rank = df[,c(rank, "db")]  %>% dplyr::filter(df[rank] != "", df[rank] != "NA") %>% dplyr::distinct()
  
  write.csv(df_rank, paste(rank,"combine_lineage_for_venn.csv", sep = "_"))
}

rank="species"



df_rank = data.frame(df_rank)

venn.plot = venn.diagram(list(df_rank[which(df_rank$db == "ITS1"),][,1], 
                              df_rank[which(df_rank$db == "ITS2"),][,1], 
                              df_rank[which(df_rank$db == "lsuD1"),][,1], 
                              df_rank[which(df_rank$db == "lsuD2"),][,1]), 
                         NULL,fill=c("red", "green", "blue", "grey"), 
                         alpha=c(0.5,0.5,0.5,0.5), 
                         cex = 2,
                         cat.fontface=4, category.names=c("ITS1", "ITS2","LSU D1", "LSU D2"))
grid.draw(venn.plot)
dev.off()


#https://www.r-bloggers.com/2019/04/set-analysis-a-face-off-between-venn-diagrams-and-upset-plots/
library(rJava)
library(tidyverse)
library(venneuler)
library(grid)
library(VennDiagram)
rank="phylum"
df_select = df[, c( rank, "db")] %>% distinct()

df = read_tsv("/home/microbiome/data_storage/temperal/R/DNA_ad/class/venn/combined_table.txt", col_names = F)
v <- venneuler(data.frame(df[,1:2]))

#pdf(paste(rank, "venn_plot.pdf", sep = "_"),height=10, width=12)
par(cex = 0.7) 
plot(v, main = "", cex.main = 5.5)
grid.text(
  "",
  x = 0.52,
  y = 0.15,
  gp = gpar(
    fontsize = 30,
    fontface = 3
  )
)
dev.off()
ggsave( paste(rank, "venn_plot.pdf", sep = "_"),height=10, width=12)

a








#################################################################################
#gene logo :https://cloud.tencent.com/developer/article/1511065
#################################################################################
#install from CRAN
#install.packages("ggseqlogo")
#install from github
#devtools::install.github("omarwagih/ggseqlogo")
#install.packages("spiralize")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require('Biostrings')) BiocManager::install("Biostrings")
if (!require('ggplot2')) install.packages('ggplot2')
if (!require('ggseqlogo')) install.packages('ggseqlogo')
if (!require('seqinr')) install.packages('seqinr')
if (!require('gridExtra')) install.packages('gridExtra')


#read the sequence
my_fasta <- read.fasta(file = "order_seq_aligned.aln", seqtype = "DNA", as.string = T)

#
my_fasta_string = vector(mode = 'character')
for (i in 1:length(my_fasta)){
  my_fasta_string[i] = toupper(c2s(my_fasta[[i]]))
}


seq_len=nchar(my_fasta_string)[1.]

write.table(seq_len, file="aligned_length.txt")


#seqlogo
#seqlogo <- ggseqlogo(my_fasta_string,seq_type = "dna", method="bits") + theme(axis.text.x = element_blank()) 
#ggsave("genelogo.pdf",dpi = 300, 
#       width=300, height=15, unit = "in",limitsize = FALSE)
#  annotate('rect', xmin = 0.5, xmax = 95.5, ymin = -0.05, ymax = 2.2, alpha = .1, col='black', fill='red')


# get sequence based on the start and end site




seq_logo <- function(my_fasta_string,seq_len){
  if(seq_len <= 400){my_fasta_string_1_100 <- subseq(my_fasta_string, start = 1,   end = 100)
  my_fasta_string_101_200 <-  subseq(my_fasta_string, start = 101, end = 200)
  my_fasta_string_201_300 <-  subseq(my_fasta_string, start = 201, end = 300)
  my_fasta_string_301_400 <-  subseq(my_fasta_string, start = 301, end = seq_len)
  
  p1 <- ggseqlogo(my_fasta_string_1_100,seq_type = "dna", method="bits") +
    theme(axis.text.x = element_blank()) +
    labs(title = "Base site 1 - 100") 
  
  p2 <- ggseqlogo(my_fasta_string_101_200,seq_type = "dna") + 
    theme(axis.text.x = element_blank()) +
    labs(title = "Base site 101 - 200")
  
  p3 <- ggseqlogo(my_fasta_string_201_300,seq_type = "dna") + 
    theme(axis.text.x = element_blank())+
    labs(title = "Base site 201 - 300")
  
  p4 <- ggseqlogo(my_fasta_string_301_400,seq_type = "dna") + 
    theme(axis.text.x = element_blank())+
    labs(title = paste("Base site 301 -",seq_len,sep=" ")) +
    annotate('rect', xmin = seq_len-299, xmax = 100, ymin = -0.05, ymax = 2.2, alpha = .1, col='white', fill='white')
  
  p <- gridExtra::grid.arrange(p1,p2,p3,p4,nrow=4)
  
  ggsave("sequence_logo.pdf",plot = p,dpi = 600, 
         width=10, height=10, unit = "in") }
  if(seq_len >= 400  & seq_len <= 500){  my_fasta_string_1_100 <-    subseq(my_fasta_string, start = 1,   end = 100)
  my_fasta_string_101_200 <-  subseq(my_fasta_string, start = 101, end = 200)
  my_fasta_string_201_300 <-  subseq(my_fasta_string, start = 201, end = 300)
  my_fasta_string_301_400 <-  subseq(my_fasta_string, start = 301, end = 400)
  my_fasta_string_401_500 <-  subseq(my_fasta_string, start = 401, end = seq_len)
  
  p1 <- ggseqlogo(my_fasta_string_1_100,seq_type = "dna", method="bits") +
    theme(axis.text.x = element_blank()) +
    labs(title = "Base site 1 - 100") 
  
  p2 <- ggseqlogo(my_fasta_string_101_200,seq_type = "dna") + 
    theme(axis.text.x = element_blank()) +
    labs(title = "Base site 101 - 200")
  
  p3 <- ggseqlogo(my_fasta_string_201_300,seq_type = "dna") + 
    theme(axis.text.x = element_blank())+
    labs(title = "Base site 201 - 300")
  
  p4 <- ggseqlogo(my_fasta_string_301_400,seq_type = "dna") + 
    theme(axis.text.x = element_blank())+
    labs(title = "Base site 301 - 400")
  
  p5 <- ggseqlogo(my_fasta_string_401_500,seq_type = "dna") + 
    theme(axis.text.x = element_blank())+
    labs(title = paste("Base site 401 -",seq_len,sep=" ")) +
    annotate('rect', xmin = seq_len-399, xmax = 100, ymin = -0.05, ymax = 2.2, alpha = .1, col='white', fill='white')
  p <- gridExtra::grid.arrange(p1,p2,p3,p4,p5,nrow=5)
  
  ggsave("sequence_logo.pdf",plot = p,dpi = 600, 
         width=10, height=10, unit = "in")}
  if(seq_len >= 500  & seq_len <= 600) {  
    my_fasta_string_1_100 <-    subseq(my_fasta_string, start = 1,   end = 100)
    my_fasta_string_101_200 <-  subseq(my_fasta_string, start = 101, end = 200)
    my_fasta_string_201_300 <-  subseq(my_fasta_string, start = 201, end = 300)
    my_fasta_string_301_400 <-  subseq(my_fasta_string, start = 301, end = 400)
    my_fasta_string_401_500 <-  subseq(my_fasta_string, start = 401, end = 500)
    my_fasta_string_501_600 <-  subseq(my_fasta_string, start = 501, end = seq_len)
    
    p1 <- ggseqlogo(my_fasta_string_1_100,seq_type = "dna", method="bits") +
      theme(axis.text.x = element_blank()) +
      labs(title = "Base site 1 - 100") 
    
    p2 <- ggseqlogo(my_fasta_string_101_200,seq_type = "dna") + 
      theme(axis.text.x = element_blank()) +
      labs(title = "Base site 101 - 200")
    
    p3 <- ggseqlogo(my_fasta_string_201_300,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = "Base site 201 - 300")
    
    p4 <- ggseqlogo(my_fasta_string_301_400,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = "Base site 301 - 400")
    
    p5 <- ggseqlogo(my_fasta_string_401_500,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = "Base site 401 - 500")
    
    p6 <- ggseqlogo(my_fasta_string_501_600,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = paste("Base site 501 -",seq_len,sep=" ")) +
      annotate('rect', xmin = seq_len-499, xmax = 100, ymin = -0.05, ymax = 2.2, alpha = .1, col='white', fill='white')
    
    p <- gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,nrow=6)
    
    ggsave("sequence_logo.pdf",plot = p,dpi = 600, 
           width=10, height=10, unit = "in")}
  if(seq_len >= 600  & seq_len <= 700) {
    my_fasta_string_1_100 <-    subseq(my_fasta_string, start = 1,   end = 100)
    my_fasta_string_101_200 <-  subseq(my_fasta_string, start = 101, end = 200)
    my_fasta_string_201_300 <-  subseq(my_fasta_string, start = 201, end = 300)
    my_fasta_string_301_400 <-  subseq(my_fasta_string, start = 301, end = 400)
    my_fasta_string_401_500 <-  subseq(my_fasta_string, start = 401, end = 500)
    my_fasta_string_501_600 <-  subseq(my_fasta_string, start = 501, end = 600)
    my_fasta_string_601_700 <-  subseq(my_fasta_string, start = 601, end = seq_len)
    
    
    p1 <- ggseqlogo(my_fasta_string_1_100,seq_type = "dna", method="bits") +
      theme(axis.text.x = element_blank()) +
      labs(title = "Base site 1 - 100") 
    
    p2 <- ggseqlogo(my_fasta_string_101_200,seq_type = "dna") + 
      theme(axis.text.x = element_blank()) +
      labs(title = "Base site 101 - 200")
    
    p3 <- ggseqlogo(my_fasta_string_201_300,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = "Base site 201 - 300")
    
    p4 <- ggseqlogo(my_fasta_string_301_400,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = "Base site 301 - 400")
    
    p5 <- ggseqlogo(my_fasta_string_401_500,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = "Base site 401 - 500")
    
    p6 <- ggseqlogo(my_fasta_string_501_600,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = "Base site 501 - 600")
    
    p7 <- ggseqlogo(my_fasta_string_601_700,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = paste("Base site 601 -",seq_len,sep=" ")) +
      annotate('rect', xmin = seq_len-599, xmax = 100, ymin = -0.05, ymax = 2.2, alpha = .1, col='white', fill='white')
    
    p <- gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=7)
    
    ggsave("sequence_logo.pdf",plot = p,dpi = 600, 
           width=10, height=10, unit = "in")
  }
  if(seq_len >= 700  & seq_len <= 800) {
    my_fasta_string_1_100 <-    subseq(my_fasta_string, start = 1,   end = 100)
    my_fasta_string_101_200 <-  subseq(my_fasta_string, start = 101, end = 200)
    my_fasta_string_201_300 <-  subseq(my_fasta_string, start = 201, end = 300)
    my_fasta_string_301_400 <-  subseq(my_fasta_string, start = 301, end = 400)
    my_fasta_string_401_500 <-  subseq(my_fasta_string, start = 401, end = 500)
    my_fasta_string_501_600 <-  subseq(my_fasta_string, start = 501, end = 600)
    my_fasta_string_601_700 <-  subseq(my_fasta_string, start = 601, end = 700)
    my_fasta_string_701_800 <-  subseq(my_fasta_string, start = 701, end = seq_len)
    
    p1 <- ggseqlogo(my_fasta_string_1_100,seq_type = "dna", method="bits") +
      theme(axis.text.x = element_blank()) +
      labs(title = "Base site 1 - 100") 
    
    p2 <- ggseqlogo(my_fasta_string_101_200,seq_type = "dna") + 
      theme(axis.text.x = element_blank()) +
      labs(title = "Base site 101 - 200")
    
    p3 <- ggseqlogo(my_fasta_string_201_300,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = "Base site 201 - 300")
    
    p4 <- ggseqlogo(my_fasta_string_301_400,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = "Base site 301 - 400")
    
    p5 <- ggseqlogo(my_fasta_string_401_500,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = "Base site 401 - 500")
    
    p6 <- ggseqlogo(my_fasta_string_501_600,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = "Base site 501 - 600")
    
    p7 <- ggseqlogo(my_fasta_string_601_700,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = "Base site 601 - 700")
    
    p8 <- ggseqlogo(my_fasta_string_701_800,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = paste("Base site 701 -",seq_len,sep=" ")) +
      annotate('rect', xmin = seq_len-699, xmax = 100, ymin = -0.05, ymax = 2.2, alpha = .1, col='white', fill='white')
    
    p <- gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,nrow=8)
    
    ggsave("sequence_logo.pdf",plot = p,dpi = 600, 
           width=10, height=10, unit = "in")
  }
  if(seq_len >=  800) {
    my_fasta_string_1_100 <-    subseq(my_fasta_string, start = 1,   end = 100)
    my_fasta_string_101_200 <-  subseq(my_fasta_string, start = 101, end = 200)
    my_fasta_string_201_300 <-  subseq(my_fasta_string, start = 201, end = 300)
    my_fasta_string_301_400 <-  subseq(my_fasta_string, start = 301, end = 400)
    my_fasta_string_401_500 <-  subseq(my_fasta_string, start = 401, end = 500)
    my_fasta_string_501_600 <-  subseq(my_fasta_string, start = 501, end = 600)
    my_fasta_string_601_700 <-  subseq(my_fasta_string, start = 601, end = 700)
    my_fasta_string_701_800 <-  subseq(my_fasta_string, start = 701, end = 800)
    my_fasta_string_801_end <-  subseq(my_fasta_string, start = 701, end = seq_len)
    
    p1 <- ggseqlogo(my_fasta_string_1_100,seq_type = "dna", method="bits") +
      theme(axis.text.x = element_blank()) +
      labs(title = "Base site 1 - 100") 
    
    p2 <- ggseqlogo(my_fasta_string_101_200,seq_type = "dna") + 
      theme(axis.text.x = element_blank()) +
      labs(title = "Base site 101 - 200")
    
    p3 <- ggseqlogo(my_fasta_string_201_300,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = "Base site 201 - 300")
    
    p4 <- ggseqlogo(my_fasta_string_301_400,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = "Base site 301 - 400")
    
    p5 <- ggseqlogo(my_fasta_string_401_500,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = "Base site 401 - 500")
    
    p6 <- ggseqlogo(my_fasta_string_501_600,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = "Base site 501 - 600")
    
    p7 <- ggseqlogo(my_fasta_string_601_700,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = "Base site 601 - 700")
    
    p8 <- ggseqlogo(my_fasta_string_701_800,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = paste("Base site 701 - 800")) 
    
    
    p9 <- ggseqlogo(my_fasta_string_801_end,seq_type = "dna") + 
      theme(axis.text.x = element_blank())+
      labs(title = paste("Base site 801 -",seq_len,sep=" ")) +
      annotate('rect', xmin = seq_len-799, xmax = 100, ymin = -0.05, ymax = 2.2, alpha = .1, col='white', fill='white')
    
    
    p <- gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,nrow=9)
    
    ggsave("sequence_logo.pdf",plot = p,dpi = 600, 
           width=10, height=10, unit = "in")
  }
}

seq_logo(my_fasta_string = my_fasta_string,seq_len = seq_len)






















########################################################################
#boxplot
##########################################################################
library(dplyr)
library(ggplot2)
library(cowplot)

rm(list=ls())
setwd("F:/PostDoc_dataset/bioinfomatic/performance evaluation")
list.files(getwd())

##########
#Read the data 1

combine <- read.table("speciesNum_DBs_combine.txt",header = T)
separate <- read.table("speciesNum_DBs_separated.txt",header = T)
df <- rbind(separate, combine)
df=combine
colnames(df)
df$Sensitivity <- df$Ture_Positive/(df$Ture_Positive + df$False_negative)
df$Accuracy <- df$Ture_Positive/(df$Ture_Positive + df$False_negative + df$False_Positive)
df$Precision <- df$Ture_Positive/(df$Ture_Positive  + df$False_Positive)
head(df)

Sensitivity <- df[,c("group_id","speciesNum","Database","length","taxa_level", "Sensitivity")]
colnames(Sensitivity) <- c("group_id","speciesNum","Database","length","taxa_level","value"  )
Sensitivity$Performance <- "Sensitivity"
Accuracy <- df[,c("group_id","speciesNum","Database","length","taxa_level","Accuracy"  )]
colnames(Accuracy) <- c("group_id","speciesNum","Database","length","taxa_level","value"  )
Accuracy$Performance <- "Accuracy"
Precision <- df[,c("group_id","speciesNum","Database","length","taxa_level","Precision"  )]
colnames(Precision) <- c("group_id","speciesNum","Database","length","taxa_level","value"  )
Precision$Performance <- "Precision"

final_table <- rbind(Sensitivity, Accuracy, Precision)
colnames(final_table)



#############################################
#generate the plot

{
species <- final_table[which(final_table$taxa_level == "Species" & final_table$Database == database & final_table$length == 130),]
colnames(species)
S <- ggplot(species, aes(speciesNum, value, fill=interaction(speciesNum,Performance), dodge=speciesNum)) +
  stat_boxplot(geom ='errorbar')+
  geom_boxplot()+
  theme_bw() +
  expand_limits(y=c(0,1)) +
  scale_fill_manual(values=c("red","red","red",
                             "sienna","sienna","sienna",
                             "royalblue2","royalblue2","royalblue2")) +
  theme(
    legend.position = 'none',
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),# 字体的大小
    axis.title.y = element_blank(),
    axis.text.y  = element_text(angle=0,# 设置旋转的角度
                                vjust=0,# 设置纵向廉价距离 hjust为横向偏移距离
                                size=12),
    #axis.ticks.x=element_blank(),
    #axis.ticks.y=element_blank()
  )




genus <- final_table[which(final_table$taxa_level == "Genus" & final_table$Database == database & final_table$length == 130),]
colnames(genus)
G <- ggplot(genus, aes(speciesNum, value, fill=interaction(speciesNum,Performance), dodge=Performance)) +
  stat_boxplot(geom ='errorbar')+
  geom_boxplot()+
  theme_bw() +
  expand_limits(y=c(0,1)) +
  scale_fill_manual(values=c("red","red","red",
                             "sienna","sienna","sienna",
                             "royalblue2","royalblue2","royalblue2")) +
  theme(
    legend.position = 'none',
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),# 字体的大小
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    #axis.ticks.x=element_blank(),
    axis.ticks.y=element_blank())




family <- final_table[which(final_table$taxa_level == "Family" & final_table$Database == database & final_table$length == 130),]
colnames(family)
F <- ggplot(family, aes(speciesNum, value, fill=interaction(speciesNum,Performance), dodge=Performance)) +
  stat_boxplot(geom ='errorbar')+
  geom_boxplot()+
  theme_bw() +
  expand_limits(y=c(0,1)) +
  scale_fill_manual(values=c("red","red","red",
                             "sienna","sienna","sienna",
                             "royalblue2","royalblue2","royalblue2")) +
  theme(
    legend.position = 'none',
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),# 字体的大小
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    #axis.ticks.x=element_blank(),
    axis.ticks.y=element_blank())



order <- final_table[which(final_table$taxa_level == "Order" & final_table$Database == database & final_table$length == 130),]
colnames(order)
O <- ggplot(order, aes(speciesNum, value, fill=interaction(speciesNum,Performance), dodge=Performance)) +
  stat_boxplot(geom ='errorbar')+
  geom_boxplot()+
  theme_bw() +
  expand_limits(y=c(0,1)) +
  scale_fill_manual(values=c("red","red","red",
                             "sienna","sienna","sienna",
                             "royalblue2","royalblue2","royalblue2")) +
  theme(
    legend.position = 'none',
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),# 字体的大小
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    #axis.ticks.x=element_blank(),
    axis.ticks.y=element_blank())



class <- final_table[which(final_table$taxa_level == "Class" & final_table$Database == database & final_table$length == 130),]
colnames(class)
C <- ggplot(class, aes(speciesNum, value, fill=interaction(speciesNum,Performance), dodge=Performance)) +
  stat_boxplot(geom ='errorbar')+
  geom_boxplot()+
  theme_bw() +
  expand_limits(y=c(0,1)) +
  scale_fill_manual(values=c("red","red","red",
                             "sienna","sienna","sienna",
                             "royalblue2","royalblue2","royalblue2")) +
  theme(
    legend.position = 'none',
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),# 字体的大小
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    #axis.ticks.x=element_blank(),
    axis.ticks.y=element_blank())



ITS_LSU = ggdraw() +
  draw_plot(S, 0,0.75,0.18,0.25) +
  draw_plot(G, 0.18,0.75,0.15,0.25) +
  draw_plot(F, 0.33,0.75,0.15,0.25) +
  draw_plot(O, 0.48,0.75,0.15,0.25) +
  draw_plot(C, 0.63,0.75,0.15,0.25)


}













