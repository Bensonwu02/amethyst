
# Function to cluster CNV profiles with k-means. Useful to separate out malignant and non-malignant cells. 
# Will probably add additional functionality to resolve subclonal structure.
cluster_cnv_profiles <- function(cnv_matrix, k = 2) {

  # k-means cannot accept NA values. Set all NA values to 1 i.e no CNV
  cnv_matrix[is.na(cnv_matrix)] <- 1

  # use two centroids to separate malignant and non-malignant
  clust <- kmeans(cnv_matrix, centers = k)

  clust_assignments <- clust$cluster
  clust_df <- data.frame(clust_assignments = clust_assignments)
  clust_df$clust_assignments <- factor(clust_df$clust_assignments, levels = c(1,2))

  return (clust_df)
  
}

# Function to normalize and log-transform each column
normalize_log_transform <- function(mat) {
  # scale counts such that they sum to 10,000
  scaled_matrix <- apply(mat, 2, function(x) {
    scaled_x <- (x / sum(x, na.rm = TRUE)) * 10000
    return(scaled_x)
  })
  
  # Log transform
  log_matrix <- log2(1 + scaled_matrix)
  
  return(log_matrix)
}

# Modified makeWindows to extract coverage
makeWindows_cov <- function(
    obj,
    type,
    genes = NULL,
    promoter = FALSE,
    stepsize = NULL,
    bed = NULL,
    index = paste0("chr_", tolower(type)),
    groupBy = NULL,
    threads = 1,
    futureType = "multicore",
    nmin = 2,
    save = FALSE) {

  # check parameters
  if (sum(!is.null(bed), !is.null(stepsize), !is.null(genes)) > 1) {
    stop("Please only specify a fixed step size, gene list, or input bed file.")
  }

  # set up multithreading
  if (threads > 1) {
    if (futureType == "multicore") {
      future::plan(future::multicore, workers = threads)
    }
    if (futureType == "multisession") {
      future::plan(future::multisession, workers = threads)
    } else if (!(futureType %in% c("multicore", "multisession"))) {
      stop("Options for parallelizing include multicore or multisession.")
    }
  }

  # check paths exist
  if (is.null(obj@h5paths)) {
    stop("Please generate the path list for each barcode and store in the obj@h5paths slot.")
  }

  # check index
  if (is.null(obj@index[[index]])) {
    stop("Please index which rows in the h5 files correspond to each chromosome using indexChr.")
  }

  # Define chromosome groups
  chr_groups <- as.list(names(obj@index[[index]]))

  by_chr <- list()

  if (!is.null(bed)) {
    # Processing input bed file
    bed <- if (is.character(bed)) {
      data.table::data.table(read.table(bed))
    } else {
      data.table::setDT(bed)
    }
    setnames(bed, c("chr", "start", "end"))
    chr_groups <- as.list(unique(as.character(bed$chr)))
    bed <- split(bed, bed$chr)

    for (chr in chr_groups) {

      # Get index positions
      sites <- obj@index[[index]][[chr]]

      # Get paths
      barcodes <- rownames(obj@h5paths[rownames(obj@h5paths) %in% sites$cell_id, , drop = FALSE])
      paths <- obj@h5paths[rownames(obj@h5paths) %in% sites$cell_id, , drop = FALSE]$paths

      # Sum c and t values for each barcode
      windows <- furrr::future_pmap(.l = list(paths, barcodes), .f = function(path, barcode) {
        tryCatch({
          bed_tmp <- copy(bed[[chr]])

          # Read data from HDF5 file for the given cell (barcode)
          h5 <- data.table::data.table(rhdf5::h5read(path, name = paste0(type, "/", barcode),
                                                     start = sites$start[sites$cell_id == barcode],
                                                     count = sites$count[sites$cell_id == barcode]))

          # Summing CpG counts (methylated + unmethylated)
          meth_window <- h5[bed_tmp, on = .(chr = chr, pos >= start, pos <= end), nomatch = 0L,
                            .(chr, start, end, total_cpg = sum(c + t, na.rm = TRUE))]
          meth_window <- meth_window[total_cpg >= nmin]
          data.table::setnames(meth_window, "total_cpg", barcode)
        }, error = function(e) {
          cat("Error processing data for barcode", barcode, ":", conditionMessage(e), "\n")
          return(NULL)
        })
      }, .progress = TRUE)

      # Merge data.tables for each barcode in chunks
      windows <- split(windows, ceiling(seq_along(windows)/1000))
      windows_merged <- Reduce(function(x, y) merge(x, y, by = c("chr", "start", "end"), all = TRUE, sort = FALSE),
                               furrr::future_map(windows, ~ Reduce(function(x, y) merge(x, y, by = c("chr", "start", "end"), all = TRUE, sort = FALSE), .x), .progress = TRUE))

      by_chr[[chr]] <- windows_merged
      cat("\nCompleted ", chr, "\n")
    }

    # Combine data from all chromosomes and format
    windows_merged <- data.table::rbindlist(by_chr, fill = TRUE, use.names = TRUE)
    data.table::setorder(windows_merged, chr, start)
    windows_merged <- windows_merged[, window := paste0(chr, "_", start, "_", end)]
    windows_merged[, c("chr", "start", "end") := NULL]
    windows_merged <- windows_merged |> tibble::column_to_rownames(var = "window")

  }

  # If stepsize is specified, create windows
  if (!is.null(stepsize)) {
    for (chr in chr_groups) {

      # Get index positions
      sites <- obj@index[[index]][[chr]]

      # Get paths
      barcodes <- rownames(obj@h5paths[rownames(obj@h5paths) %in% sites$cell_id, , drop = FALSE])
      paths <- obj@h5paths[rownames(obj@h5paths) %in% sites$cell_id, , drop = FALSE]$paths

      windows <- furrr::future_pmap(.l = list(paths, barcodes), .f = function(path, barcode) {
        tryCatch({
          # Read data from HDF5 file for the given cell (barcode)
          h5 <- data.table::data.table(rhdf5::h5read(path, name = paste0(type, "/", barcode),
                                                     start = sites$start[sites$cell_id == barcode],
                                                     count = sites$count[sites$cell_id == barcode]))

          # Summing CpG counts (methylated + unmethylated)
          h5 <- h5[pos %% stepsize == 0, pos := pos + 1]
          meth_window <- h5[, window := paste0(chr, "_", plyr::round_any(pos, stepsize, floor), "_", plyr::round_any(pos, stepsize, ceiling))]
          meth_window <- meth_window[, .(total_cpg = sum(c + t, na.rm = TRUE)), by = window]
          meth_window <- meth_window[total_cpg >= nmin]
          data.table::setnames(meth_window, "total_cpg", barcode)
        }, error = function(e) {
          cat("Error processing data for barcode", barcode, ":", conditionMessage(e), "\n")
          return(NULL)
        })
      }, .progress = TRUE)

      # Merge data.tables per cell in chunks
      windows <- split(windows, ceiling(seq_along(windows)/1000))
      windows_merged <- Reduce(function(x, y) merge(x, y, by = "window", all = TRUE, sort = FALSE),
                               furrr::future_map(windows, ~ Reduce(function(x, y) merge(x, y, by = "window", all = TRUE, sort = FALSE), .x), .progress = TRUE))

      if (save) {
        saveRDS(windows_merged, paste0("tmp_", type, "_", chr, "_", stepsize, "kb_windows_nmin", nmin, ".RData"))
      }
      by_chr[[chr]] <- windows_merged
      cat("\nCompleted ", chr, "\n")
    }

    windows_merged <- data.table::rbindlist(by_chr, fill = TRUE, use.names = TRUE)

    # Order and clean up the data
    windows_merged <- windows_merged[, c("chr", "start", "end") := {
      split_parts <- data.table::tstrsplit(window, "_", fixed = TRUE)
      list(split_parts[[1]], as.numeric(split_parts[[2]]), as.numeric(split_parts[[3]]))
    }]
    windows_merged <- windows_merged[(end - start) == stepsize]
    data.table::setorder(windows_merged, chr, start)
    windows_merged[, c("chr", "start", "end") := NULL]
    windows_merged <- windows_merged |> tibble::column_to_rownames(var = "window")
  }

  if (threads > 1) {
    future::plan(future::sequential)
    gc()
  }

  return(windows_merged)
}


#' CNV inference using amethyst object
#' 
#' @param obj Amethyst object
#' @param reference_cells Vector containing names of reference cells
#' @param query_cells Vector containing names of query cells
#' @param step_size Size of genomic window
#' @param num_threads Number of threads for multithreading
#' @param plot_dir Output directory
#' @param nmin Minimum number of observations to keep genomic window
#' @param k Number of clusters to extract from CNV profiles. k = 2 is good for malignant vs non-malignant
#' @param cluster_cells Whether to cluster cells or not
#' @param save Whether to save CNV matrix
#' 
#' @return Heatmap of query cell CNV profiles
#' 
cnvInference <- function(obj, reference_cells, query_cells, step_size = 50*1e6, num_threads = 4, plot_dir, nmin = 3, k = 2, cluster_cells = FALSE, save = FALSE, matrix_name = "cnv_matrix") {

  require(ComplexHeatmap)
  require(colorRamp2)
  require(RColorBrewer)
  
  # Get CpG coverage matrix 
  obj@genomeMatrices[["raw_cpg_cov"]] <- makeWindows_cov(
    obj,
    type = "CG",
    genes = NULL,
    promoter = FALSE,
    stepsize = step_size,
    bed = NULL,
    index = "chr_cg",
    groupBy = NULL,
    threads = num_threads,
    futureType = "multicore",
    nmin = nmin,
    save = FALSE)

  print(head(obj@genomeMatrices[["raw_cpg_cov"]]))
  saveRDS(obj, file.path(obj_dir, "cnv_test.rds"))

  sum_matrix <- obj@genomeMatrices[["raw_cpg_cov"]]

  sum_matrix <- as.data.frame(normalize_log_transform(sum_matrix))

  # Select reference and query cells based on log-transformed total coverage
  ref_matrix <- sum_matrix %>% 
    select(contains(paste0(reference_cells)))

  query_matrix <- sum_matrix %>% 
    select(contains(paste0(query_cells)))

  # Compute log2 ratios between query and reference cells
  log2_ratios <- query_matrix

  for (cell in query_cells) {

    query_col <- query_matrix[[paste0(cell)]]

    # Compute log2 ratio between query cell and reference mean
    log2_ratio <- query_col - rowMeans(ref_matrix, na.rm = TRUE)

    log2_ratios[[paste0("ratio_", cell)]] <- log2_ratio
  }

  # Prepare matrix for heatmap (log2 ratios of query cells)
  ratio_matrix <- log2_ratios %>%
    select(contains("ratio"))

  # Reverse log transform 
  ratio_matrix <- 2^ratio_matrix

  chromosomes <- gsub("_.*", "", rownames(ratio_matrix))
  ratio_matrix$chr <- chromosomes
  
  # Remove sex chromosomes
  ratio_matrix <- ratio_matrix %>% 
    filter(!(chr %in% c("chrX", "chrY"))) 

  # Order chromosomes
  chromosome_numbers <- as.numeric(gsub("chr", "", ratio_matrix$chr))
  ratio_matrix$chr <- chromosome_numbers
  ratio_matrix <- ratio_matrix %>% 
    arrange(chr) %>% 
    select(-chr)

  windows <- rownames(ratio_matrix)
  print(windows)

  # Transpose so rows are cells and columns are genomic windows
  cnv_matrix <- as.data.frame(t(ratio_matrix))
  colnames(cnv_matrix) <- windows

  if (cluster_cells == TRUE) {

    set.seed(111)
    # Cluster CNVs to classify tumor cells
    clust_df <- cluster_cnv_profiles(cnv_matrix, k = k)

    # Export cluster assignments
    write.csv(clust_df, file.path(plot_dir, "cnv_cluster_assignment.csv"))

    row_annot <- rowAnnotation(`Cluster` = clust_df$clust_assignments,
      show_annotation_name = FALSE,
      border = TRUE,
      col = list(`Cluster` = c("1" = "black", "2" = "grey")))
  }
  
  if (save == TRUE) {
    # Export cnv matrix
    write.csv(cnv_matrix, file.path(plot_dir, paste0(matrix_name, ".csv")))  
  }

  color.scheme <- rev(brewer.pal(9,"RdBu"))
  colors <- colorRamp2(c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2), color.scheme)

  chromosome_split <- factor(sort(chromosome_numbers), levels = c(seq(1,22)))
  # Export chromosome split
  saveRDS(chromosome_split, file.path(plot_dir, "chromosome_split.rds"))

  chr_cols <- brewer.pal(22, "Set3")
  chr_cols <- c(chr_cols, chr_cols[1:10])
  names(chr_cols) <- unique(chromosome_split)
  col_annot <- HeatmapAnnotation(`Chromosome` = chromosome_split,
    show_annotation_name = FALSE,
    border = TRUE,
    col = list(`Chromosome` = chr_cols),
    height = unit(0.5, "cm"))

  if (cluster_cells == TRUE) {
    ht <- Heatmap(as.matrix(cnv_matrix),
              column_split = chromosome_split,
              left_annotation = row_annot,
              top_annotation = col_annot,
              cluster_rows = TRUE, 
              cluster_columns = FALSE,
              show_row_names = FALSE,
              clustering_method_rows = "ward.D2",
              show_column_names = FALSE,
              col = colors,
              name = "CNV",
              row_dend_width = unit(1, "cm"),
              border_gp = gpar(col = "black", lty = "solid"),
              column_gap = unit(0, "cm"),
              width = unit(30, "cm"),
              height = unit(15, "cm")
              )
  } else {
    ht <- Heatmap(as.matrix(cnv_matrix),
              column_split = chromosome_split,
              top_annotation = col_annot,
              cluster_rows = TRUE, 
              cluster_columns = FALSE,
              show_row_names = FALSE,
              clustering_method_rows = "ward.D2",
              show_column_names = FALSE,
              col = colors,
              name = "CNV",
              row_dend_width = unit(1, "cm"),
              border_gp = gpar(col = "black", lty = "solid"),
              column_gap = unit(0, "cm"),
              width = unit(30, "cm"),
              height = unit(15, "cm")
              )
  }
  
  return(ht)
}