# Import MICs and Covariates

Import MICs and Covariates

## Usage

``` r
import_mics_with_metadata(
  data,
  mic_column,
  metadata_columns = NULL,
  code_column = NULL,
  combination_agent = 0,
  log_reg_value = FALSE,
  scale = "log",
  round = FALSE,
  include_mic_bounds = FALSE,
  low_con = NULL,
  high_con = NULL,
  concentration_by_covariate = NULL
)
```

## Arguments

- data:

  Data frame containing the MICs and covariates

- mic_column:

  String, name of column in data corresponding to the MIC values

- metadata_columns:

  Vector of column names (as strings) for covariates to be included in
  the data frame produced

- code_column:

  String, name of column containing any symbols for MICs (if data is in
  a 2 column format)

- combination_agent:

  Numerical, if the MIC is not a combination agent, 0. If combination
  agent must be separated by '/', use 1 to select the value before the
  '/', or 2 for the value after. Recommend a value where log2(value) is
  an integer.

- log_reg_value:

  Logical, TRUE if a column for logistic regression model should be
  included in output (MICs with \`\>\` will be doubled, MICs with
  \`\<\`, \`\<=\`, or \`≤\` are halved)

- scale:

  String, "log" if MIC value should be converted to log2 scale
  (preferred for subsequent fitting of the model using the EM algorithm)

- round:

  Set to true if log2(MIC values) are integers, but decimal MIC values
  are rounded (e.g. 0.12 in place of 0.125)

- include_mic_bounds:

  Logical, if TRUE includes left and right boundaries of interval on MIC
  scale (in addition to on log2 scale if scale is "log")

- low_con:

  Numerical, the lowest concentration tested, specify on the same scale
  as the data. If left null and concentration_by_covariate table is not
  supplied, will be set based on the data

- high_con:

  Numerical, the highest concentration tested, specify on the same scale
  as the data. If left null and concentration_by_covariate table is not
  supplied, will be set based on the data

- concentration_by_covariate:

  Data frame. Table, columns "low_con" and "high_con" are matched to
  data using covariates. Include one row for each combination of
  covariates.

## Examples

``` r
import_mics_with_metadata(data = tibble::tibble(MIC_A = c("≤0.12", ">16", 4, 2), t = runif(4, 0, 10)),
                         mic_column = "MIC_A",
                         metadata_columns = "t",
                         log_reg_value = TRUE,
                         scale = "log",
                         round = TRUE)
#> # A tibble: 4 × 8
#>   obs_id left_bound right_bound mic_column lr_column     t low_con high_con
#>    <int>      <dbl>       <dbl> <chr>          <dbl> <dbl>   <dbl>    <dbl>
#> 1      1       -Inf          -3 ≤0.12           0.12  5.93      -3        4
#> 2      2          4         Inf >16            32     6.87      -3        4
#> 3      3          1           2 4               4     7.35      -3        4
#> 4      4          0           1 2               2     4.61      -3        4

import_mics_with_metadata(data = tibble::tibble(MIC_A = c(0.125, 16, 4, 2), code_A = c("≤", ">", NA, NA), t = runif(4, 0, 10)),
                         mic_column = "MIC_A",
                         metadata_columns = "t",
                         code_column = "code_A",
                         log_reg_value = FALSE,
                         scale = "log",
                         round = FALSE,
                         include_mic_bounds = TRUE)
#> # A tibble: 4 × 10
#>   obs_id left_bound right_bound mic_column code_column left_bound_mic
#>    <int>      <dbl>       <dbl> <chr>      <chr>                <dbl>
#> 1      1       -Inf          -3 0.125      ≤                        0
#> 2      2          4         Inf 16         >                       16
#> 3      3          1           2 4          NA                       2
#> 4      4          0           1 2          NA                       1
#> # ℹ 4 more variables: right_bound_mic <dbl>, t <dbl>, low_con <dbl>,
#> #   high_con <dbl>

import_mics_with_metadata(data = tibble::tibble(MIC_A = c("≤10/1", ">80/8", "40/4", "20/2"), t = runif(4, 0, 10)),
                         mic_column = "MIC_A",
                         metadata_columns = "t",
                         combination_agent = 2,
                         log_reg_value = FALSE,
                         scale = "log",
                         round = FALSE)
#> # A tibble: 4 × 8
#>   obs_id left_bound right_bound mic_column code_column     t low_con high_con
#>    <int>      <dbl>       <dbl> <chr>      <chr>       <dbl>   <dbl>    <dbl>
#> 1      1       -Inf           0 1          <=          3.06        0        3
#> 2      2          3         Inf 8          >           7.32        0        3
#> 3      3          1           2 4          NA          1.34        0        3
#> 4      4          0           1 2          NA          0.998       0        3
```
