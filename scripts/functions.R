### all functions
### _ Loading in functions
cat("Load in the functions")

# ! TODO: remove all occurences where a function first merges it with some df
## ! instead, have seperate merging functions, e.g., merge_annotations, that do
## ! this first and then use function to carry out task

## - functions

merge_cmap <- function(
    df,
    cmap = CMAP_coloc,
    by_join = "SAMPLEID") {
  merged_df <- df %>%
    left_join(
      cmap,
      by = by_join
    )

  return(merged_df)
}

merge_annotations <- function(
    df,
    annotation_table = annoNifHDB_updt,
    by_join = "AUID") {
  merged_df <- df %>%
    left_join(
      annotation_table,
      by = by_join
    )

  return(merged_df)
}



# remove_variables <- function(pattern_to_remove) {
#    # List all variables in the environment
# all_vars <- ls()

# # Find variables containing "pattern" in their names and remove them
# to_remove <- all_vars[grep(pattern_to_remove, all_vars)]
# rm(list = to_remove)
# }

# remove_variables("temp");ls()

transform_data_lng <- function(
    input_df,
    starts_with_col,
    names_to_col,
    values_to_col) {
  if ("group_id" %in% names(input_df)) {
    result <- input_df %>%
      select(-group_id) %>%
      pivot_longer(
        cols = starts_with(starts_with_col),
        names_to = starts_with_col,
        values_to = values_to_col
      )
  } else {
    result <- input_df %>%
      pivot_longer(
        cols = starts_with(starts_with_col),
        names_to = names_to_col,
        values_to = values_to_col
      )
  }
  return(result)
}


add_total_row <- function(
    df,
    summary_column = NULL, # name of column summarizing
    column_name,
    pull_name = NULL,
    all_columns) {
  if (all_columns == FALSE) {
    df_with_total_count <- df %>%
      bind_rows(
        tibble(
          !!column_name := "Total",
          {{ summary_column }} := sum(df %>%
            pull({{ summary_column }}))
        )
      )
  } else {
    df_with_total_count <- df %>%
      bind_rows(
        tibble(
          !!column_name := "Total",
          df %>%
            summarise(across(where(is.numeric), ~ mean(.)))
        )
      )
  }

  return(df_with_total_count)
}

cat(
  'USEAGE:
add_total_row(
    column_name = "studyID",
    summary_column = n,
    all_columns = FALSE
    )'
)

### adds group id using the unique_sample_id_key
###
add_group_id <- function(df, var_select = SAMPLEID, by_var = "SAMPLEID") {
  df_with_group_id <- df %>%
    left_join(unique_sample_id_key %>%
      select({{ var_select }}, group_id), by = by_var)

  return(df_with_group_id)
}

dedup_by_group <- function(df, ...) {
  contains_group_id <- "group_id" %in% names(df)

  if (contains_group_id) {
    df_deup <- df %>%
      # filter(!duplicated(group_id))
      distinct(..., .keep_all = TRUE)
  } else {
    df_deup <- df %>%
      add_group_id() %>%
      # filter(!duplicated(group_id))
      # distinct(distinct_var1, distinct_var2, distinct_var3)
      distinct(..., .keep_all = TRUE)
  }

  return(df_deup)
}

main_average_ra_dedup_by_group <- function(df_lng, ..., mean_by) {
  df_lng_mean_ra_depup_by_group <- df_lng %>%
    add_group_id() %>%
    group_by(...) %>%
    mutate(
      # Calculate averages for the relevant columns
      average_value_AUID = mean({{ mean_by }}, na.rm = T)
      # across(starts_with("AUID"), ~ mean(., na.rm = T))
    ) %>%
    ungroup() %>%
    dedup_by_group(...)

  return(df_lng_mean_ra_depup_by_group)
}

remove_samples_nucleic_acid <- function(df, nucleic_acid_type, nucleic_acid_samples_key = DNA_samples_key) {
  nucleic_acid_flag <- c("DNA", "RNA") %in% nucleic_acid_type
  if (nucleic_acid_flag[1] == TRUE) {
    nucleic_acid <- df %>%
      filter(SAMPLEID %in% {{ nucleic_acid_samples_key }})
    cat("Only", nucleic_acid_type, "will be returned

    ")
  } else if (nucleic_acid_flag[2] == TRUE) {
    nucleic_acid <- df %>%
      filter(!SAMPLEID %in% {{ nucleic_acid_samples_key }})
    cat("Only", nucleic_acid_type, "will be returned

    ")
  } else {
    cat("You supplied '", nucleic_acid_type, ",' which is neither DNA or RNA. Please make sure it is in ALL CAPS


    ", sep = "")
    nucleic_acid <- NULL
  }

  return(nucleic_acid)
}

count_and_arrange <- function(data, group_vars, arrange_var = n) {
  df_count_and_arrange <- data %>%
    count(across(all_of(group_vars))) %>%
    arrange(desc({{ arrange_var }}))

  return(df_count_and_arrange)
}

add_percentage <- function(
    df,
    sum_by_percent,
    percentage_id,
    grouping_by = NULL,
    remove_columns = NULL) {
  # total_column_name <- paste0("total_", {{sum_by_percent}})
  df_with_total <- df %>%
    group_by({{ grouping_by }}) %>%
    # group_by(across(all_of(grouping_by ))) %>%
    mutate(
      total = sum({{ sum_by_percent }}),
      {{ percentage_id }} := {{ sum_by_percent }} / total * 100
    ) %>%
    ungroup()

  if (!is.null({{ remove_columns }})) {
    df_with_total <- df_with_total %>%
      select(-{{ remove_columns }})
  }
  # select(select_by, sum_by_percent, percentage)
  # select({{ select_by }})

  return(df_with_total)
}


sum_tax <- function(
    abundance_table,
    annotation_table,
    joining_by_col,
    total_counts_id,
    grouping_for_percentage = NULL,
    percentage_id,
    grouping_by_1 = NULL,
    grouping_by_2 = NULL,
    sum_by,
    ...) {
  result <- abundance_table %>%
    left_join(annotation_table, joining_by_col) %>%
    group_by({{ grouping_by_1 }}, {{ grouping_by_2 }}) %>%
    # group_by({{ grouping_by_1 }}) %>%
    mutate({{ total_counts_id }} := sum({{ sum_by }}, na.rm = T)) %>%
    ungroup() %>%
    distinct({{ total_counts_id }}, {{ grouping_by_1 }}, {{ grouping_by_2 }}, .keep_all = T) %>%
    # distinct({{ total_counts_id }}, ..., .keep_all = T) %>%
    arrange(desc({{ total_counts_id }})) %>%
    group_by({{ grouping_for_percentage }}) %>%
    mutate(
      total = sum({{ total_counts_id }}),
      {{ percentage_id }} := {{ total_counts_id }} / total * 100
    ) %>%
    ungroup() %>%
    # add_percentage(
    #   sum_by_percent = total_counts,
    #   percentage_id = "study_id_percentage",
    # ) %>%
    select(..., {{ total_counts_id }}, total, {{ percentage_id }}, everything()) # %>%
  # ungroup() # %>%
  # select(..., percentage)

  return(result)
}

clean_percentages <- function(
    df,
    grouping_by,
    clean_var,
    threshold = 1.0,
    # sum_by,
    percentage_id,
    column_var_mutate,
    column_var) {
  cleaned_df <- df %>%
    # group_by({{ grouping_by }}) %>%
    # summarise(sum(study_id_percentage))  %>%
    mutate(
      {{ column_var_mutate }} := ifelse({{ clean_var }} < threshold, "other", {{ column_var }}),
      # percentage_id = ifelse(percentage_id < 1.0, sum(percentage_id), percentage_id)
    ) %>%
    group_by({{ column_var_mutate }}, {{ grouping_by }}) %>%
    summarise({{ percentage_id }} := sum({{ clean_var }})) %>%
    mutate(
      # percentage_id = sprintf("%.1f", percentage_id)
      {{ percentage_id }} := round({{ percentage_id }}, 1)
    ) %>%
    ungroup() %>%
    # arrange(desc(percentage_id))  %>%
    arrange(desc({{ percentage_id }}))

  return(cleaned_df)
}



### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
cat("Done loading script!!!")
cat("Woooooooohooooooo!!!")
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
### _ Finished loading in the data ### _ Finished loading in the data
