#This script imputes missing data using random forests.

if (!suppressPackageStartupMessages(require("pacman"))) { install.packages("pacman") } #if pacman isn't already installed, install it.

pacman::p_load(
  tidyverse, 
  missForest, 
  reshape2
)

data <- read_tsv("data/Sahul_structure_wide.tsv",col_types = cols()) %>% 
  dplyr::rename(ID = glottocode) #this column is already aggregated for dialects in make_wide.R

#removing non-binary traits and traits with all 0
binary_feats <- data %>% 
  reshape2::melt(id.vars = "ID") %>% 
  distinct(variable, value) %>% 
  group_by(variable) %>% 
  summarise(sum = sum(value, na.rm = T)) %>% 
  filter(sum == 1) %>% 
  distinct(variable) %>% 
  as.matrix() %>% 
  as.vector()

data <- data %>% 
  dplyr::select(ID, all_of(binary_feats))

missing_prop <- data %>% 
  dplyr::select(-ID) %>%
  is.na() %>%
  mean()

cat(sprintf("%0.2f%% of the data read was missing before imputation.\n",
            missing_prop * 100
))


# remove languages with less than LG_NA_PROP_CUTOFF data
# 0 = good (no missing data), 1 = bad (entirely missing)
LG_NA_PROP_CUTOFF = 0.5

# remove features present in less than FEAT_NA_PROP_CUTOFF languages
FEAT_NA_PROP_CUTOFF = 0.5

#counting the missing data per col, i.e. per feature
col_na_means <- data %>% 
  as.matrix() %>% 
  is.na() %>% 
  colMeans() 

cols_to_keep <- tibble(Feature_ID = colnames(data),
                       feature_na_prop = col_na_means) %>% 
  filter(Feature_ID != "ID") %>% 
  arrange(-feature_na_prop) %>% 
  distinct(Feature_ID, .keep_all = T) %>% 
  filter(feature_na_prop < FEAT_NA_PROP_CUTOFF) %>% 
  dplyr::select(Feature_ID) %>% 
  as.matrix() %>% 
  as.vector()

#Some languages have a large amount of missing values. Let's exclude those with a large amount of missing data

data$na_prop <- data %>% 
  column_to_rownames("ID") %>% 
  apply(1, function(x) mean(is.na(x)))

#df with the features and observations with the fewest missing as per cut offs.
data_few_missing_values <- data %>% 
  filter(na_prop < LG_NA_PROP_CUTOFF ) %>% 
  dplyr::select(-na_prop) %>% 
  dplyr::select(ID, all_of(cols_to_keep))


#In order to make sure the imputation is doing categorically and not for continuous integers the values are replaces by strings. yes, I know this is silly but I've found it simplfies things, unfortunately.
data_prepped_for_imputation <- data_few_missing_values %>% 
  reshape2::melt(id.vars = "ID") %>% 
  mutate(value = str_replace_all(value, "0", "0 - absent")) %>%
  mutate(value = str_replace_all(value, "1", "1 - present")) %>%
  mutate(value = str_replace_all(value, "2", "2 - multistate")) %>%
  mutate(value = str_replace_all(value, "3", "3 - multistate")) %>%
  mutate(value = str_replace_all(value, "4", "4 - multistate")) %>%
  mutate(value = as.factor(value)) %>%
  dplyr::select(ID, variable, value) %>%
  spread(key = variable, value, drop = FALSE) 

missing_value_percentage_after_cropping <- mean(is.na(data_prepped_for_imputation))
cat(sprintf(
  "The data has now been cropped, and there remains %0.2f%% missing values.\n",
  missing_value_percentage_after_cropping * 100
))

cat("These are the values that will be imputed by random forests.\n")

#Saving the language IDs separately, it's better to leave them out of the imputation workflow
IDs <- data_prepped_for_imputation$ID

#actual imputation
imputed_data <- data_prepped_for_imputation %>%
  dplyr::select(-ID) %>%
  as.matrix() %>%
  data.frame() %>%
  mutate_all(as.factor) %>% 
  missForest() 

imputed_data$OOBerror

df_imputed <- imputed_data$ximp
df_imputed$ID <- IDs

df_imputed_num <- df_imputed %>%
  mutate_all(as.character) %>% 
  reshape2::melt(id = "ID")  %>%
  mutate(value =  str_replace_all(value, "0 - absent", "0")) %>%
  mutate(value =  str_replace_all(value,  "1 - present", "1")) %>%
  mutate(value =  str_replace_all(value,  "2 - multistate", "2")) %>%
  mutate(value =  str_replace_all(value,  "3 - multistate", "3")) %>%
  mutate(value =  str_replace_all(value,  "4 - multistate", "4")) %>%
  mutate(value = as.numeric(value)) %>% 
  dplyr::select(ID, Parameter_ID = variable, value) %>%
  spread(key = Parameter_ID, value = value, drop = FALSE) 

cat("The imputation of missing values is done. ", round(mean(is.na(data_prepped_for_imputation))*100,2), "% of the data was imputed. 
    Note: we did not impute all missing values in the entire set. Before imputation the dataset was cropped to a subset where languages with more than ",   LG_NA_PROP_CUTOFF *100 ,"% missing values were removed and features with more than",  FEAT_NA_PROP_CUTOFF*100,"% missing languages were removed. This left us with", df_imputed_num %>% nrow(), "languages and", df_imputed_num[,-1] %>% ncol(),"features

In the full dataset, there was",round(missing_prop, 2)*100, "% of the data missing.

After this cropping the dataset contained ", round(mean(is.na(data_prepped_for_imputation))*100,2), "% missing data, this is what was imputed. The Out of Bag error or for the imputation is", round(imputed_data$OOBerror, 2), "

After imputation, the dataset contains",round(mean(is.na(df_imputed_num)*100,2)), "% missing values and is ready for analysis that requires a complete dataset (PCA or otherwise).")

df_imputed_num %>% 
  write_tsv("data/Sahul_structure_wide_imputed.tsv")
