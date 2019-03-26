#Jenny Smith 
#3/6/19
#Purpose: Clean the manifest for SMRT-seq 


#Set-up environment
suppressPackageStartupMessages(library(dplyr, quietly = TRUE))
suppressPackageStartupMessages(library(magrittr, quietly = TRUE))
suppressPackageStartupMessages(library(tidyr, quietly = TRUE))
library(methods, quietly = TRUE)

options(stringsAsFactors = FALSE)
path <- Sys.getenv("OUTDIR")
sample.info.dir <- "~/scripts/Isoseq3_workflow/samples" #not sure if I should hard-code this?

#check that outdir is not empty. 
if(path == ""){
  print("ERROR: No Directory Specified")
  quit(save = "no", status = 1, runLast = FALSE)
}


#Sample ID map from the bash script manifest.sh
samples <- read.csv(paste(path, "Sample_ID_Map.csv", sep="/"), header=FALSE)
# head(samples)
# dim(samples)

#Sample molecular annotations
info <- read.csv(paste(sample.info.dir,"libraries_submitted_to_FHCRC_genomics_core_for_PacBio.csv", sep="/"))
# head(info)
# dim(info)

#Clean up the manifest and merge in the sample molecular annotations.
update <- samples %>%
  set_colnames(c("Filepath", "Dataset", "Sample_ID_Original")) %>%
  separate(col = Sample_ID_Original,
           into = c("Reg.", "AmpureBeads_Conc", "Description"), sep = "[\\._-]",
           remove = FALSE, extra="merge") %>%
  mutate(Reg.=case_when(
    Reg.=="Sample10" ~ "NBM",
    Reg.=="84631" ~ "846361",
    TRUE ~ Reg.)) %>%
  mutate_at(vars(AmpureBeads_Conc), funs(case_when(
    AmpureBeads_Conc == "4" ~ "0.4x",
    AmpureBeads_Conc == "1" ~ "1x",
    grepl("1X", Description) ~ "1x",
    grepl("0.4X", Description) ~ "0.4x"))) %>%
  full_join(., info, by=c("Reg."="Sample_Name")) %>%
  mutate(Comments=ifelse(grepl("84631.4",Sample_ID_Original),
                         "NOTE: the Sample ID was mis-typed in the subreads.xml file. Need to Follow-up and check on sample.", "")) %>%
  select(Filepath, Dataset, Reg.,  Sample_Number=X, Sample_Info, AmpureBeads_Conc, Sample_ID_Original,Comments) %>%
  arrange(Sample_Number,Sample_Info)


#Save to file
write.csv(update, paste(sample.info.dir,"Sample_ID_Map.csv", sep="/"), row.names = FALSE,quote = FALSE)
