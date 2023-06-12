# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

# function that returns csv names and number of rows to skip in each csv

library(readr)
library(dplyr)

read_n_skip_fun <- function(folder_address_of_csv_files) {
  file_list <- list.files(folder_address_of_csv_files, pattern = "*.csv", full.names = TRUE)

  skip_lines <- sapply(file_list, function(x) grep("prot_hit_num", readr::read_lines(x))[1] - 1, USE.NAMES = FALSE)

  list(csv_location = file_list, no_lines_to_skip = skip_lines)
}

check_columns <- function(file_location, skip_lines, columns) {
  dat <- read_csv(file_location, skip = skip_lines, col_types = cols())

  dat <- select(dat, all_of(columns))

  if (sum(names(dat) == columns) == length(columns)) {
    return("good")
  } else {
    return("bad")
  }
}

combined_function <- function(folder_address) {
  list1 <- read_n_skip_fun(folder_address)
  col_to_check <- c('prot_hit_num', 'prot_acc', 'pep_seq', 'prot_seq', 'pep_var_mod', 'pep_var_mod_pos')

  count_good <- 0

  for (file in seq_along(list1$no_lines_to_skip)) {
    col_check_result <- check_columns(list1$csv_location[file], list1$no_lines_to_skip[file], col_to_check)

    if (col_check_result == "good") {
      count_good <- count_good + 1
    }
  }

  print(paste("Number of files:", length(list1$csv_location)))
  print(paste("Number of 'good' results:", count_good))
}

raw_to_count <- function(folder_address) {
  list1 <- read_n_skip_fun(folder_address)

  # Extract the folder name from the folder_address
  folder_name <- basename(folder_address)

  # Create the output folder path
  output_folder <- file.path(folder_address, 'new after conversion', 'data')

  # Create the subfolders if they don't exist
  dir.create(file.path(output_folder, 'peptide info'), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_folder, 'protein sequence with modification'), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_folder, 'whole protein sequence with modification'), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_folder, 'ratio data'), recursive = TRUE, showWarnings = FALSE)

  for (file in seq_along(list1$no_lines_to_skip)) {
    #################
    ## now reading files
    final_list_of_prot_seq <- list()
    complete_mod_prot_seq <- list()
    list_peptide_info <- list()

    dat <- read_csv(list1$csv_location[file], skip = list1$no_lines_to_skip[file], col_types = cols())

    # removing all duplicate and na
    dat %>%
      select(prot_hit_num, prot_acc, pep_seq, prot_seq, pep_var_mod, pep_var_mod_pos) %>%
      na.omit() %>%
      distinct(prot_hit_num, pep_seq, prot_seq, pep_var_mod, pep_var_mod_pos, .keep_all = TRUE) %>%
      mutate(end_of_pep_in_prot = str_locate(pattern = pep_seq, prot_seq)[,2]) %>%
      filter(!is.na(end_of_pep_in_prot)) %>%
      mutate(pep_var_mod_pos = str_extract(pep_var_mod_pos,"(?<=..).+(?=..)"),
             pep_var_mod_pos = gsub('2', 0, pep_var_mod_pos),
             pep_var_mod_pos = gsub('4', 0, pep_var_mod_pos)) -> dat

    # after removing c and m
    # removing strings which only contains 0

    a <- vector()
    for (i in 1:nrow(dat)) {
      a[i] <- sum(str_split(dat$pep_var_mod_pos[i],'',simplify = T) == '0') / (length(str_split(dat$pep_var_mod_pos[i],'',simplify = T)))
    }

    dat$new_pep_var_mod_pos <- a
    dat$new_pep_var_mod_pos <- ifelse(dat$new_pep_var_mod_pos == 1, NA, dat$pep_var_mod_pos)

    dat %>% na.omit() -> dat

    # finding unique number of protein seq (via prot_hit_num)
    uniq_ID <- data.frame(prot_acc = unique(dat$prot_acc))

    for (n in 1:nrow(uniq_ID)) {

      dat %>%
        filter(prot_acc == uniq_ID[n,1]) %>%
        mutate(length_of_prot = nchar(prot_seq),
               length_of_pep = nchar(pep_seq)) %>%
        na.omit() -> xyz5

      # start and end position of peptide in protein sequence
      #xyz5$end_of_pep_in_prot <- str_locate(pattern = xyz5$pep_seq, xyz5$prot_seq)[,2]


      # calculate the postion of identifier in peptide
      final_list <- lapply(lapply((strsplit(xyz5$pep_var_mod_pos,'')), as.numeric), function(x) which(x != 0))
      find_new_pep_pos <- plyr::ldply(final_list, rbind)

      # calculate the number of modifications
      xyz5$number_of_pos_of_mod <- unlist(lapply(final_list, function (x) length(x)))

      # calculating postion of peptide modification in protein
      final_pos_of_pep_in_pro <- xyz5$end_of_pep_in_prot - (xyz5$length_of_pep - find_new_pep_pos)

      # which identifier modified in peptide
      which_pep_mod <- plyr::ldply(lapply(lapply((strsplit(xyz5$pep_var_mod_pos,'')), as.numeric), function(x) x[x != 0]), rbind)

      names(final_pos_of_pep_in_pro) <- paste0('pos_', names(final_pos_of_pep_in_pro))
      names(which_pep_mod) <- paste0('which_amino_', names(which_pep_mod))

      # making matrix of protein seq as column and peptide as row

      start_n_end_of_pep_in_prot <- str_locate(pattern = xyz5$pep_seq, xyz5$prot_seq)
      pep_start <- start_n_end_of_pep_in_prot[,1]
      pep_end <- start_n_end_of_pep_in_prot[,2]

      list_peptide_info[[n]] <- data.frame(pep_seq = xyz5$pep_seq, final_pos_of_pep_in_pro, pep_start=pep_start, pep_end=pep_end, which_pep_mod, pep_var_mod_pos=xyz5$pep_var_mod_pos) %>% arrange(pos_1)
      names(list_peptide_info) <- uniq_ID$prot_acc[n]

      matching_mat <- matrix(NA, ncol = (nchar(unique(xyz5$prot_seq))), nrow = nrow(xyz5))
      colnames(matching_mat) <- unlist(strsplit(unique(xyz5$prot_seq),''), use.names = FALSE)
      rownames(matching_mat) <- xyz5$pep_seq

      for (i in 1:nrow(start_n_end_of_pep_in_prot)) {
        #matching_mat[i,pep_start[i]:pep_end[i]] <- 'A'
        matching_mat[i,(as.numeric(final_pos_of_pep_in_pro[i,1:xyz5$number_of_pos_of_mod[i]]))] <- as.numeric(which_pep_mod[i,1:xyz5$number_of_pos_of_mod[i]])
      }

      # check unique values in all columns

      uniq_list <- apply(matching_mat, 2,function(x) unique(x))

      #sum(unlist(lapply(uniq_list, function(x) length(x) > 2)))

      if (sum(unlist(lapply(uniq_list, function(x) length(x) > 2))) == 0) {
        print('No different modification in same columns: Good ')
      } else{
        print('Error: different modification in same columns: Bad')
      }

      #unlist(lapply(lapply(uniq_list, function(x) x[!(is.na(x))]), as.numeric))
      freq_table <- t(as.matrix(unlist(lapply(lapply(uniq_list, function(x) x[!(is.na(x))]), function (x) ifelse(length(x) > 0, x, NA)))))
      freq_table[is.na(freq_table)] <- 0


      #final_list_of_prot_seq <- list()
      final_list_of_prot_seq[[n]] <- freq_table
      complete_mod_prot_seq[[n]] <- append(paste0(colnames(freq_table)[apply(freq_table, MARGIN = 1, function(x) x!=0)],':',which(freq_table != 0)),paste0('prot_len',':',ncol(freq_table)))

      print(n)
    }

    names(list_peptide_info) <- uniq_ID$prot_acc
    names(complete_mod_prot_seq) <- uniq_ID$prot_acc
    names(final_list_of_prot_seq) <- uniq_ID$prot_acc

    # Save the results in the output folder
    saveRDS(list_peptide_info, file = file.path(output_folder, 'peptide info', paste0(str_split(list1$csv_location[file], pattern = "/", simplify = TRUE)[4], '.rds')))
    saveRDS(complete_mod_prot_seq, file = file.path(output_folder, 'protein sequence with modification', paste0(str_split(list1$csv_location[file], pattern = "/", simplify = TRUE)[4], '.rds')))
    saveRDS(final_list_of_prot_seq, file = file.path(output_folder, 'whole protein sequence with modification', paste0(str_split(list1$csv_location[file], pattern = "/", simplify = TRUE)[4], '.rds')))

    # complete_protein_seq_with_mod[[file]] <- complete_mod_prot_seq
    # peptide_info[[file]] <- list_peptide_info

    # names(complete_protein_seq_with_mod)[file] <- str_split(use_for_loop$csv_location[file], pattern = "/", simplify = TRUE)[5]
    # names(peptide_info)[file] <- str_split(use_for_loop$csv_location[file], pattern = "/", simplify = TRUE)[5]

    total_R_prot <- unlist(lapply(final_list_of_prot_seq, function(x) sum(colnames(x) == 'R')))
    total_C_prot <- unlist(lapply(final_list_of_prot_seq, function(x) sum(colnames(x) == 'C')))
    total_K_prot <- unlist(lapply(final_list_of_prot_seq, function(x) sum(colnames(x) == 'K')))
    total_M_prot <- unlist(lapply(final_list_of_prot_seq, function(x) sum(colnames(x) == 'M')))
    total_P_prot <- unlist(lapply(final_list_of_prot_seq, function(x) sum(colnames(x) == 'P')))
    total_T_prot <- unlist(lapply(final_list_of_prot_seq, function(x) sum(colnames(x) == 'T')))

    total_R_pep <- unlist(lapply(final_list_of_prot_seq, function(x) sum(x == 1)))
    total_C_pep <- unlist(lapply(final_list_of_prot_seq, function(x) sum(x == 2)))
    total_K_pep <- unlist(lapply(final_list_of_prot_seq, function(x) sum(x == 3)))
    total_M_pep <- unlist(lapply(final_list_of_prot_seq, function(x) sum(x == 4)))
    total_P_pep <- unlist(lapply(final_list_of_prot_seq, function(x) sum(x == 5)))
    total_T_pep <- unlist(lapply(final_list_of_prot_seq, function(x) sum(x == 6)))

    ratio_R <- (total_R_pep/total_R_prot)
    ratio_C <- (total_C_pep/total_C_prot)
    ratio_K <- (total_K_pep/total_K_prot)
    ratio_M <- (total_M_pep/total_M_prot)
    ratio_P <- (total_P_pep/total_P_prot)
    ratio_T <- (total_T_pep/total_T_prot)

    data.frame(uniq_ID,
               total_R_prot ,
               total_R_pep ,
               ratio_R ,

               total_C_prot ,
               total_C_pep ,
               ratio_C ,

               total_K_prot ,
               total_K_pep ,
               ratio_K ,

               total_M_prot ,
               total_M_pep ,
               ratio_M ,

               total_P_prot ,
               total_P_pep ,
               ratio_P ,

               total_T_prot ,
               total_T_pep ,
               ratio_T ) -> ratio_matrix

    write.csv(ratio_matrix, paste0(output_folder, '/ratio data', '/ratio',sub(".*/", "_", list1$csv_location)[file]), row.names = F)
    #prot_names[[file]] <- ratio_matrix$prot_acc

    rm(matching_mat, ratio_matrix, uniq_list, freq_table, list_peptide_info, final_list_of_prot_seq, complete_mod_prot_seq)
  }
}
