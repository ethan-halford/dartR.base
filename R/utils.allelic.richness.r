# # bootstrapping function
all.rich <- function(df,
                    indices,
                    boot_method = "loc") {

  # all.rich_fun <- function(x){
  #   # reshape matrix into long format
  #   long_data <- melt(x, varnames = c("ind", "site"), value.name = "genotype")
  # 
  #   # Tally and compute reference/alternate allele counts
  #   allele_counts <- long_data %>%
  #     group_by(site, genotype) %>%
  #     tally(name = "n") %>%
  #     na.omit() %>%
  #     mutate(
  #       ref_allele = case_when(
  #         genotype == 0 ~ n * 2,
  #         genotype == 1 ~ n,
  #         TRUE          ~ 0
  #       ),
  #       alt_allele = case_when(
  #         genotype == 2 ~ n * 2,
  #         genotype == 1 ~ n,
  #         TRUE          ~ 0
  #       )
  #     ) %>%
  #     group_by(site) %>%
  #     summarise(
  #       ref_total  = sum(ref_allele),
  #       alt_total  = sum(alt_allele),
  #       raw_count  = sum(ref_allele + alt_allele),
  #       .groups    = "drop"
  #     )
  # 
  #   # Rarefaction calculation
  #   min_sample_size <- min(allele_counts$raw_count, na.rm = TRUE)
  # 
  #   # Compute mean corrected richness across sites
  #   allele_counts_fin <- allele_counts %>%
  #     mutate(
  #       r_ref = 1 - choose(raw_count - ref_total, min_sample_size) / choose(raw_count, min_sample_size),
  #       r_alt = 1 - choose(raw_count - alt_total, min_sample_size) / choose(raw_count, min_sample_size)
  #     ) %>%
  #     summarise(mean_richness = mean(r_ref + r_alt, na.rm = TRUE)) %>%
  #     pull(mean_richness) %>%
  #     round(6)
  # 
  #   return(allele_counts_fin)
  # }

  Rcpp::cppFunction(
    '
double allelicRichness(NumericMatrix x) {
  int nrow = x.nrow();
  int ncol = x.ncol();
  
  // Vectors to store per-site values
  std::vector<int> rawCounts;      // total allele count per site
  std::vector<double> refTotals;   // sum of reference allele counts per site
  std::vector<double> altTotals;   // sum of alternate allele counts per site
  
  // Loop over sites (columns)
  for (int j = 0; j < ncol; j++) {
    int count0 = 0, count1 = 0, count2 = 0;
    // Loop over individuals (rows)
    for (int i = 0; i < nrow; i++) {
      double val = x(i, j);
      if (R_IsNA(val)) continue; // skip NA values
      if (val == 0) count0++;
      else if (val == 1) count1++;
      else if (val == 2) count2++;
    }
    int nonNA = count0 + count1 + count2;
    if(nonNA == 0) continue; // skip sites with no data
    int raw_count = 2 * nonNA; // each individual contributes 2 alleles
    double ref_total = 2.0 * count0 + 1.0 * count1;
    double alt_total = 2.0 * count2 + 1.0 * count1;
    
    rawCounts.push_back(raw_count);
    refTotals.push_back(ref_total);
    altTotals.push_back(alt_total);
  }
  
  int m = rawCounts.size();
  if (m == 0) return NA_REAL; // if no sites with data
  
  // Determine the minimum raw count (minimum sample size) across sites
  int min_sample_size = rawCounts[0];
  for (int i = 1; i < m; i++) {
    if (rawCounts[i] < min_sample_size)
      min_sample_size = rawCounts[i];
  }
  
  double sumRichness = 0.0;
  
  // For each site, compute r_ref and r_alt using the hypergeometric rarefaction formula:
  // r = 1 - choose(raw - allele_total, min_sample_size) / choose(raw, min_sample_size)
  // Instead of using choose() directly, we compute the ratio as a product.
  for (int i = 0; i < m; i++) {
    int r_count = rawCounts[i];
    double ref_total = refTotals[i];
    double alt_total = altTotals[i];
    
    double prod_ref = 1.0;
    // If (r_count - ref_total) is less than min_sample_size, then choose(...) is zero
    if ((r_count - ref_total) < min_sample_size) {
      prod_ref = 0.0;
    } else {
      for (int k = 0; k < min_sample_size; k++) {
        prod_ref *= ((r_count - ref_total - k) / static_cast<double>(r_count - k));
      }
    }
    double r_ref = 1.0 - prod_ref;
    
    double prod_alt = 1.0;
    if ((r_count - alt_total) < min_sample_size) {
      prod_alt = 0.0;
    } else {
      for (int k = 0; k < min_sample_size; k++) {
        prod_alt *= ((r_count - alt_total - k) / static_cast<double>(r_count - k));
      }
    }
    double r_alt = 1.0 - prod_alt;
    
    double site_richness = r_ref + r_alt;
    sumRichness += site_richness;
  }
  
  double meanRichness = sumRichness / m;
  // Round the result to 6 decimal places
  double factor = 1e6;
  double rounded = std::round(meanRichness * factor) / factor;
  
  return rounded;
}
    '
  )


  df <- df[indices,]

  if(boot_method == "loc"){
    df <- t(df)
  }

  res <- allelicRichness(df)

  return(res)

}

# Rcpp::cppFunction(
# '
#   inline double round_n(double x, int n) {
#     double factor = std::pow(10.0, n);
#     return std::floor(x * factor + 0.5) / factor;
#   }
# 
# double all_rich_cpp(NumericMatrix df,
#                     IntegerVector indices,
#                     std::string boot_method = "loc")
# {
#   int n_idx = indices.size();
#   int ncol  = df.ncol();
# 
#   NumericMatrix subdf(n_idx, ncol);
#   for(int i = 0; i < n_idx; i++) {
#       int ridx = indices[i] - 1;
#       for(int j = 0; j < ncol; j++) {
#         subdf(i, j) = df(ridx, j);
#       }
#   }
# 
# NumericMatrix data;
# if(boot_method == "loc") {
#   NumericMatrix subdfT(ncol, n_idx);
#   for(int i = 0; i < n_idx; i++) {
#     for(int j = 0; j < ncol; j++) {
#       subdfT(j, i) = subdf(i, j);
#     }
#   }
#   data = subdfT;
# } else {
#   data = subdf;
# }
# 
# int nrow_final = data.nrow();
# int ncol_final = data.ncol();
# 
#   std::vector<double> ref_allele(nrow_final, 0.0);
# std::vector<double> alt_allele(nrow_final, 0.0);
# std::vector<double> total_allele(nrow_final, 0.0);
# 
# for(int i = 0; i < nrow_final; i++) {
#   for(int j = 0; j < ncol_final; j++) {
#     double val = data(i, j);
#     if(NumericMatrix::is_na(val)) continue;
# 
#     int genotype = static_cast<int>(val);
#     switch(genotype) {
#       case 0:
#       ref_allele[i]   += 2.0;
#       total_allele[i] += 2.0;
#       break;
#       case 1:
#       ref_allele[i]   += 1.0;
#       alt_allele[i]   += 1.0;
#       total_allele[i] += 2.0;
#       break;
#       case 2:
#       alt_allele[i]   += 2.0;
#       total_allele[i] += 2.0;
#       break;
#       default:
#         break;
#       }
#     }
#   }
# 
#   double min_sample_size = R_PosInf;
#   for(int i = 0; i < nrow_final; i++) {
#     if(total_allele[i] > 0.0 && total_allele[i] < min_sample_size) {
#       min_sample_size = total_allele[i];
#     }
#   }
#   if(min_sample_size == R_PosInf) {
#     return NA_REAL;
#   }
# 
#   std::vector<double> site_richness(nrow_final, NA_REAL);
# 
#   for(int i = 0; i < nrow_final; i++) {
#     double raw_count = total_allele[i];
#     double r_ref = NA_REAL;
#     double r_alt = NA_REAL;
# 
#       if(raw_count >= min_sample_size) {
#         double denom = Rf_choose(raw_count, min_sample_size);
# 
#         if(!ISNAN(denom) && denom != 0.0) {
#           double num_ref = Rf_choose(raw_count - ref_allele[i], min_sample_size);
#           double num_alt = Rf_choose(raw_count - alt_allele[i], min_sample_size);
# 
#           if(!ISNAN(num_ref)) {
#             r_ref = 1.0 - (num_ref / denom);
#           }
#           if(!ISNAN(num_alt)) {
#             r_alt = 1.0 - (num_alt / denom);
#           }
#         }
#       }
# 
#       if(!ISNAN(r_ref) && !ISNAN(r_alt)) {
#         site_richness[i] = r_ref + r_alt;
#       }
#     }
# 
# double sum_val = 0.0;
# int count_val  = 0;
# for(int i = 0; i < nrow_final; i++) {
#   if(!ISNAN(site_richness[i])) {
#     sum_val += site_richness[i];
#     count_val++;
#   }
# }
# 
# if(count_val == 0) {
#   return NA_REAL;
# }
# 
# double mean_richness = sum_val / count_val;
# mean_richness = round_n(mean_richness, 6);
# 
# return mean_richness;
#   }'
# )
