# Required packages #

library(Biostrings)
library(stringr)
library(stringi)
library(seqinr)
library(progress)
library(plyr)
library(ggplot2)
library(ggpubr)

# Install if you want the BLAST functionality #
install.packages('rBLAST', repos = 'https://mhahsler.r-universe.dev')
library(rBLAST)

####################################################################################################################################

#       These steps can be performed for any bacterial genome, as long as reference genomes and full CDS annotations are provided. #
#       All steps taken for Xcc below can be swapped for other bacteria                                                            #

####################################################################################################################################

# Make database for blast #

ref_genome_Xcc8004 <- readDNAStringSet('Xcc_8004_assembly.fna')

makeblastdb("Xcc_8004_assembly.fna", dbtype='nucl')
blast_db_Xcc8004 <- blast(db = "Xcc_8004_assembly.fna", type = "blastn")


# Load fasta file with sequences #

Xcc8004_fasta <- readDNAStringSet("Xcc_8004_genomic_cds.fa")

# Select species to find targets for (can be a list of species to run in sequence) # 

species_names <- c("Xcc8004")

# Adjust according to species names above #
species_list <- list(Xcc8004_fasta)

####################################################################################################################################

#                                                   gRNA target finder function                                                    #

####################################################################################################################################

gRNA_find <- function(data, targets, plot_histogram, blast, blast_data){
  oligo_df <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), 
                       c("species", "gene", "protein", "protein_ID", "target_strand", "STOP_position","oligo_fw","oligo_rv","cds"))
  oligo_df2 <- setNames(data.frame(matrix(ncol = 10, nrow = 0)), 
                        c("species", "gene", "protein", "protein_ID", "target_strand", "STOP_position","oligo_fw","oligo_rv","cds","target_nr"))
  oligo_ph <- setNames(data.frame(matrix(ncol = 9, nrow = 1)), 
                       c("species", "gene", "protein", "protein_ID", "target_strand", "STOP_position","oligo_fw","oligo_rv","cds"))
  oligo_rc_ph <- setNames(data.frame(matrix(ncol = 9, nrow = 1)), 
                          c("species", "gene", "protein", "protein_ID", "target_strand", "STOP_position","oligo_fw","oligo_rv","cds"))
  
  species <- species_names[sp]
  
  pb1 <- progress_bar$new(format = "[:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta] Finding targets in :total genes...",
                          total = length(data), complete = "=", incomplete = " ", current = ">", clear = FALSE, width = 110)
  
  ### Find targets on coding strand ###
  
  for(i in 1:length(data)){
    pb1$tick()
    
    # Add gene ID and protein to dataset #
    
    elements <- unlist(strsplit(names(data)[i], "\\["))
    elements_clean <- str_trim(str_replace_all(elements, "\\[|\\]", ""))
    
    # Extract column header names and IDs #
    
    element_names <- sapply(strsplit(elements_clean[2:(length(elements_clean))], "="), function(x) x[1])
    element_IDs <- sapply(strsplit(elements_clean[2:(length(elements_clean))], "="), function(x) x[2])
    
    # Select columns for locus_tag, protein and protein_id #
    
    gene_ID <- element_IDs[which(element_names == "locus_tag")]
    protein <- element_IDs[which(element_names == "protein")]
    protein_ID <- element_IDs[which(element_names == "protein_id")]
    
    # Find codons that can be edited into STOP codons along the CDS #
    cod <- as.character(codons(DNAString(paste(data[i]))))
    cds <- as.character(data[i])
    cod_nr <- which(cod %in% c("CAG","CAA","CGA"))
    cod_pos <- (cod_nr-1)*3+1
    cod_pos <- cod_pos[cod_pos <= (nchar(cds) - 22)]
    if(length(cod_pos) > 0){
      pam_reg <- substring(cds, cod_pos + 17, cod_pos + 22)
      cod_sel <- which(grepl("GG", pam_reg))
      if(length(cod_sel > 0)){
        cod_pos <- cod_pos[cod_sel]
        pam_pos <- unlist(lapply(gregexpr2("GG",pam_reg[cod_sel]), function(x) min(x)))
        
        oligo <- paste("TAGC",unique(substring(cds, cod_pos + (pam_pos - 5), cod_pos + (pam_pos + 14))), sep = "")
        oligo_rv <- c()
        for(j in 1:length(oligo)){
          oligo_rv <- c(oligo_rv, paste("AAAC", as.character(reverseComplement(DNAString(substring(oligo, 5, 25)[j]))), sep = ""))
        }
        for(k in 1:targets){
          oligo_ph[k,] <- c(species, gene_ID, protein, protein_ID, "Coding", (cod_pos[k]+2)/3, oligo[k], oligo_rv[k], cds)
        }
        oligo_df <- rbind(oligo_df, oligo_ph)
      }
    }
    # Targets on template strand #
    rc_cod_nr <- which(cod == "TGG")
    rc_cod_pos <- (rc_cod_nr-1)*3+1
    rc_cod_pos <- rc_cod_pos[rc_cod_pos > 21]
    if(length(rc_cod_pos) > 0){
      rc_pam_reg <- substring(cds, rc_cod_pos - 21, rc_cod_pos - 15)
      rc_cod_sel <- which(grepl("CC", rc_pam_reg))
      if(length(rc_cod_sel > 0)){
        rc_cod_pos <- rc_cod_pos[rc_cod_sel]
        rc_pam_pos <- unlist(lapply(gregexpr2("CC",rc_pam_reg[rc_cod_sel]), function(x) max(x)))
        rc_pam_loc <- rc_cod_pos - 20 + rc_pam_pos
        rc_oligo <- paste("AAAC", unique(substring(cds, rc_pam_loc + 1 , rc_pam_loc + 20)),sep = "")
        rc_oligo_rv <- c()
        for(j in 1:length(rc_oligo)){
          rc_oligo_rv <- c(rc_oligo_rv, paste("TAGC", as.character(reverseComplement(DNAString(substring(rc_oligo, 5, 25)[j]))), sep = ""))
        }
        for(k in 1:targets){
          oligo_rc_ph[k,] <- c(species, gene_ID, protein, protein_ID, "Template", (rc_cod_pos[k]+2)/3, rc_oligo[k], rc_oligo_rv[k], cds)
        }
        oligo_df <- rbind(oligo_df, oligo_rc_ph)
      }
    }
  }
  
  ## Create final dataframe with oligo candidates ##
  
  colnames(oligo_df) <-  c("species","gene","protein","protein_ID","target_strand","STOP_position","oligo_fw","oligo_rv", "cds")
  for(i in unique(oligo_df$gene)){
    oligo_df_sub <- subset(oligo_df, gene == i)
    oligo_df_sub$STOP_position <- as.numeric(oligo_df_sub$STOP_position)
    oligo_df_sub <- oligo_df_sub[order(oligo_df_sub$STOP_position, decreasing = F, na.last = T),]
    oligo_df_sub$target_nr <- seq(1:min(c(nrow(oligo_df_sub),(targets*2))))
    oligo_df_sub <- subset(oligo_df_sub, target_nr <= targets)
    oligo_df2 <- rbind(oligo_df2, oligo_df_sub)
  }
  
  ## Histogram ##
  
  if(plot_histogram == "Yes"){
    targets_df <- data.frame(width=width(data), seq=as.character(data), names=names(data), row.names = 1:length(data))
    targets_df$gene <- substring(targets_df$names,unlist(gregexpr2("locus_tag=", targets_df$names)) + 10, nchar(targets_df$names))
    targets_df$gene <- substring(targets_df$gene,1, sapply(gregexpr2("]", targets_df$gene), "[[", 1)-1)                             
    
    yes_targets <- as.numeric(length(which(targets_df$gene %in% unique(oligo_df2$gene))))
    no_targets <- nrow(targets_df) - yes_targets
    
    targets_df$target <- paste("No targets (", no_targets, "/", nrow(targets_df), " genes)", sep = "")
    targets_df$target[targets_df$gene %in% unique(oligo_df2$gene)] <- paste("At least one target (", yes_targets, "/", nrow(targets_df), " genes)", sep = "")
    
    # Plot histogram of genes with targets vs. no targets, divided by ORF length
    target_hist <- gghistogram(targets_df, 
                               x = "width", color = "target", fill = "target", bins = 200, rug = TRUE, combine = FALSE, facet.by = "target")+
      scale_fill_manual(values = c("green","red"))+
      scale_color_manual(values = c("darkgreen","darkred"))+
      scale_x_continuous(breaks = seq(0,round_any(max(targets_df$width), 1000, f = ceiling),1000))+
      labs(y = "Count", x = "ORF length (bp)")+
      theme(strip.text.x = element_text(face = "bold"),
            strip.text.y = element_text(face = "bold"),
            axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), size = 14),
            axis.text.y = element_text(size = 10),
            axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 14),
            axis.text.x = element_text(angle = 45, vjust = 0.75, hjust = 1, size = 10),
            legend.position = "none",
            panel.grid.major.y = element_line(color = "grey90"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white", color = "black", linewidth = 1),
            plot.margin = unit(c(1,1,1,1), "cm"))
    
  }
  
  ## BLAST ##
  
  if(blast == "Yes"){
    oligo_df2 <- na.omit(oligo_df2)
    pb2 <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta] Blasting targets...",
                            total = nrow(oligo_df2), complete = "=", incomplete = " ", current = ">", clear = FALSE, width = 110) 
    for(i in 1:nrow(oligo_df2)){
      pb2$tick()
      query <- DNAStringSet(substring(oligo_df2$oligo_fw[i],5,24))
      best_off_target <- predict(blast_data, query , BLAST_args = "-task 'blastn-short' -perc_identity 90")
      oligo_df2$off_target_match[i] <- best_off_target[2,4]
      oligo_df2$off_target_identity[i] <- best_off_target[2,3]
      oligo_df2$off_target_start[i] <- best_off_target[2,9]
      oligo_df2$off_target_end[i] <- best_off_target[2,10]
      oligo_df2$e_value[i] <- best_off_target[2,11]
    }
  }
  print(paste("CBE gRNA targets found for ", length(unique(oligo_df$gene)), " out of ", length(data), " genes", sep = ""))
  print(target_hist)
  return(list(oligo_df2, target_hist))
}

####################################################################################################################################

#                                                    RUN THE CBE FINDER SCRIPT                                                     #

####################################################################################################################################

# Run finder: arguments (data = data, targets = nr. of targets you want to find per gene, plot_histogram = plot "Yes/No", blast = "Yes/No")

for(sp in 1:length(species_names)){
  assign(paste0("results_", species_names[sp]),
         gRNA_find(data = species_list[[sp]],
                   targets = 3,
                   plot_histogram = "Yes",
                   blast = "No",
                   blast_data = blast_db_Xcc8004))
}

data <- as.data.frame(results_Xcc8004[1])

                      