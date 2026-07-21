# simulate_mics

function that wraps together all the functions that determine t's
distribution, pi and its trends, trends in the mean, component draws,
epsilon, covariates, and censors the data

## Usage

``` r
simulate_mics(
  n = 300,
  t_dist = function(n) {
     runif(n, min = 0, max = 16)
 },
  pi = function(t) {
z <- 0.17 + 0.025 * t - 0.00045 * t^2
     tibble(`1` = 1 - z,
    `2` = z)
 },
  `E[X|T,C]` = function(t, c) {
     case_when(c == "1" ~ -4 + (0.24 * t) - (0.0055 *
    t^2), c == "2" ~ 3 + 0.001 * t, TRUE ~ NaN)
 },
  sd_vector = c(`1` = 1, `2` = 1.05),
  covariate_list = NULL,
  covariate_effect_vector = c(0),
  conc_limits_table = NULL,
  low_con = -3,
  high_con = 6,
  scale = "log"
)
```

## Arguments

- n:

  Number of observations

- t_dist:

  A function of n for drawing values of t

- pi:

  A function of time that returns a vector of weights that sum to 1.

- sd_vector:

  A vector with length equal to the number of components, with the
  elements named "1", "2",...

- covariate_list:

  List of covariates, each one has its own format, see examples of
  numeric and categorical covariates

- covariate_effect_vector:

  Vector of covariate effects corresponding to the covariates listed
  above

- conc_limits_table:

  If concentration limits vary by some covariate use this table to
  specify limits for each value of the covariate. Is right-joined to
  data by the covariate values

- low_con:

  If concentration limits are constant for all observations, used to set
  the lowest tested concentration on the log2(MIC) scale

- high_con:

  If concentration limits are constant for all observations, used to set
  the highest tested concentration on the log2(MIC) scale

- scale:

  What scale ("log" or "MIC") the data returned by simulate_mics is.
  Default is "log" which corresponds to log2(MIC)

- \`E\[X\|T, C\]\`:

  A function of time and component that returns a value of mu (component
  mean) for any given time and component

## Examples

``` r
#Covariate List
covariate_list = list(c("numeric", "normal", 0, 1), c("categorical", c(0.3, 0.4)), c("numeric", "uniform", 0, 5))
#Covariate Effect Vector
covariate_effect_vector = c(2, #intercept for all covariates combined
                                                                10, #slope for covariate_1
                                                                100, #effect of level b vs a of covariate 2
                                                                3 #slope for covariate_3
                                                                )
#Concentration Limits Table
conc_limits_table = tibble::as_tibble(rbind(c("a", -3, 3),
                                    c("b", -4, 4),
                                    c("c", -4, 4)),`.name_repair` = "unique"
) |> dplyr::rename("covariate_2" = 1, "low_cons" = 2, "high_cons" = 3)
#> New names:
#> • `` -> `...1`
#> • `` -> `...2`
#> • `` -> `...3`

simulate_mics(
n = 300,
t_dist = function(n){runif(n, min = 0, max = 16)},
pi = function(t) {
  z <- 0.17 + 0.025 * t - 0.00045 * t ^ 2
  tibble::tibble("1" = 1 - z, "2" = z)
},
`E[X|T,C]` = function(t, c)
{
  dplyr::case_when(c == "1" ~ -4.0 + (0.24 * t) - (0.0055 * t ^ 2),
            c == "2" ~ 3 + 0.001 * t,
            TRUE ~ NaN)
},
sd_vector = c("1" = 1, "2" = 1.05),
covariate_list = covariate_list,
covariate_effect_vector = covariate_effect_vector,
conc_limits_table = conc_limits_table,
scale = "log")
#> New names:
#> • `` -> `...1`
#> • `` -> `...2`
#> • `` -> `...3`
#> Joining with `by = join_by(covariate_2)`
#> Joining with `by = join_by(t, p, comp, x, sd, epsilon, covariate_1,
#> covariate_2, covariate_3, total_cov_effect, observed_value, low_cons,
#> high_cons)`
#> # A tibble: 300 × 17
#>        t p         comp      x    sd epsilon covariate_1 covariate_2 covariate_3
#>    <dbl> <list>    <chr> <dbl> <dbl>   <dbl>       <dbl> <chr>             <dbl>
#>  1 12.2  <dbl [2]> 1     -1.89  1     -0.350      0.363  b                  2.63
#>  2 10.7  <dbl [2]> 2      3.01  1.05   1.47       1.44   b                  4.15
#>  3  9.74 <dbl [2]> 2      3.01  1.05   0.179      1.08   b                  3.34
#>  4 11.4  <dbl [2]> 1     -1.98  1      0.704      0.140  b                  3.53
#>  5  2.40 <dbl [2]> 1     -3.46  1      0.608     -0.550  a                  3.81
#>  6 15.8  <dbl [2]> 1     -1.58  1     -0.533     -0.0287 a                  2.62
#>  7  8.44 <dbl [2]> 1     -2.37  1     -0.903      0.553  b                  1.78
#>  8  3.28 <dbl [2]> 1     -3.27  1      2.51      -0.692  a                  1.03
#>  9  4.06 <dbl [2]> 1     -3.12  1     -0.586      0.774  b                  3.53
#> 10 15.8  <dbl [2]> 2      3.02  1.05  -0.191      1.17   b                  4.82
#> # ℹ 290 more rows
#> # ℹ 8 more variables: total_cov_effect <dbl[,1]>, observed_value <dbl[,1]>,
#> #   tested_concentrations <list>, left_bound <dbl>, right_bound <dbl>,
#> #   indicator <dbl>, low_con <dbl>, high_con <dbl>

```
