library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)
library(tibble)
library(RColorBrewer)
library(car)
library(purrr)
library(cowplot)

#directories used. change as needed
TCGA_datasets <- "C:\\Users\\thoma\\Desktop\\TNBC_project\\datasets\\figure_data\\TCGA_data"
our_datasets <- "C:\\Users\\thoma\\Desktop\\TNBC_project\\datasets\\figure_data\\our_data"
FIGURE_IMAGE_WD <- "C:\\Users\\thoma\\Desktop\\TNBC_project\\figures"

#named constants for our plot size:

STAT_COMPARE_SIZE <- 15
LABELSIZE <- 12
HLINE_SIZE <- 1
FONTSIZE <- 48


#common functions:
setwd(TCGA_datasets)
read_mutation_df <- function(df_name, variantname)
  #read in file format of mutation calling output
  #params: df_name: path to file of mutation calls
  #variantName: mutation type called in file, "snp" or "indel"
{
  df <- read.table(df_name)
  df <- df %>% mutate(V2=substr(V2, 3,5))
  colnames(df) <- c(variantname, "sample")
  return(df)
}

add_subtype_our_dataset <- function(df)
  #add lehmann subtype and BRCA labels to each sample in df 
  #returns labeled df 
{
  translation <- read.csv(paste(our_datasets, "clinical_RNAseq_translation.csv", sep="\\")) %>% mutate(V3=as.character(V3), V1=as.character(V1))
  if("20" %in% df$sample){
    #if we are using collaborator's ID-ing system.
    translation <- translation %>% dplyr::select(V1, V2, V4)
    df<- inner_join(translation, df, by=c("V1"= "sample"))
  } else {
    translation <- translation %>% dplyr::select(V2,V3, V4)
    df<- inner_join(translation, df, by=c("V3"= "sample"))
  }
  return(df)
}

my_theme <- function(plot, f=FONTSIZE){
  #common theme used in all plots
  #params: plot: plot to apply theme to. f: font size to be used in plot
  plot <- plot+theme_classic()+theme(text = element_text(face="bold",size=f, color="black"), legend.text=element_text(face="bold",size=f, color="black")) + theme(legend.position="none", axis.text.y = element_text(margin=ggplot2::margin(l = 6), size=f-10), axis.title.y = element_text(margin=ggplot2::margin(r=30), size=f-12), axis.text.x=element_text(angle=-45, vjust = .3, hjust=.5))
  return(plot)
}

my_theme_without_x <- function(plot, f=FONTSIZE){
  #common theme used in all plots
  #params: plot: plot to apply theme to. f: font size to be used in plot
  plot <- my_theme(plot) 
  plot <- plot + theme(axis.title.x =element_blank())
  return(plot)
}


add_subtypes_breast_cbioportal<- function(df, df_id_colnum=1){
  #add subtype and BRCA labels for each sample to breast df
  #Params: df: dataframe to add labels to. df_id_column: col number of column with TCGA ids
  try(if(df_id_colnum < 0 || df_id_colnum > ncol(df)) stop("df id column should be col number of column with TCGA data"))
  subtypes <- read.csv("TNBC_subtypes_summarized.csv", stringsAsFactors = FALSE)
  #lehmanns list of subtypes
  colnames(df)[df_id_colnum] <- "BARCODE"
  df <- full_join(df, subtypes, by="BARCODE")
  brca1_carriers <- read.table("TNBC_union_brca1_carriers.txt", stringsAsFactors = FALSE) %>% unlist()
  brca2_carriers <- read.table("TNBC_union_brca2_carriers.txt", stringsAsFactors = FALSE) %>% unlist()
  df <- df %>% mutate(TNBCtype_seperated=ifelse(BARCODE %in% brca1_carriers, "BRCA1", ifelse(BARCODE %in% brca2_carriers, "BRCA2", TNBCtype)))
  df<- df %>% mutate(TNBCtype_seperated=factor(TNBCtype_seperated, levels=c("BRCA1","BRCA2", "IM", "BL1", "BL2", "LAR", "M", "MSL", "ER", "HER2")))
  df <- df %>% mutate(carrier= ifelse(TNBCtype_seperated %in% c("BRCA1", "BRCA2"), "carrier", "noncarrier"))
  df <- df %>% filter(!(TNBCtype_seperated %in% c("HER2", "ER")))
  df <- arrange(df, TNBCtype_seperated)
  return(df)
}


add_subtypes_ovarian <- function(df, df_id_colnum=1)
{ 
  #add subtype and BRCA labels for each sample to ovarian df
  #Params: df: dataframe to add labels to. df_id_column: col number of column with TCGA ids
  try(if(df_id_colnum < 0 || df_id_colnum > ncol(df)) stop("df id column should be col number of column with TCGA data"))
  colnames(df)[df_id_colnum] <- "ID"
  subtypes <- read.csv("CLOVAR_classifications_TCGA.csv", stringsAsFactors = FALSE)
  subtypes <- subtypes %>% dplyr::select(ID, SUBTYPE)
  df <- full_join(df, subtypes, by="ID")
  brca1 <- read.csv("ovarian_brca1.csv", header=FALSE)$V1 %>% substr(1,12)
  brca2 <- read.csv("ovarian_brca2.csv", header=FALSE)$V1 %>% substr(1,12)
  df <- df %>% mutate(SUBTYPE=ifelse(ID %in% brca1, "brca1", ifelse(ID %in% brca2, "brca2", SUBTYPE)))
  df<- df %>% mutate(SUBTYPE=factor(SUBTYPE, levels=c("brca1","brca2","Immunoreactive", "Differentiated", "Mesenchymal", "Proliferative"))) %>% arrange(SUBTYPE)
  df$SUBTYPE <- dplyr::recode(.x = df$SUBTYPE,Immunoreactive= "IMR", Differentiated= "DIF", Mesenchymal= "MES", Proliferative= "PRO", brca1= "BRCA1", brca2= "BRCA2")
  df <- df %>% mutate(carrier=ifelse(SUBTYPE %in% c("BRCA1", "BRCA2"), "carrier", "noncarrier"))
  return(df)
}


save_plot_list <- function(plot_list, set_width=13.9, set_height=10.9){
  setwd(FIGURE_IMAGE_WD)
  for (i in 1:length(plot_list)){
    ggsave(filename = paste("plot",letters[i], ".png", sep=""), plot=plot_list[[i]], width=set_width, height=set_height)
  }
}

