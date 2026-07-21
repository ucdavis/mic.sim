# Title

Title

## Usage

``` r
log_reg(
  data,
  split_by = "ecoff",
  data_type,
  drug,
  date_col,
  date_type,
  first_year,
  s_breakpoint,
  r_breakpoint,
  ecoff,
  visual_split,
  k = NULL
)
```

## Arguments

- data:

  either possible_data or a set of mic data ready you want to run
  import_mics on

- data_type:

  use either "possible_data" if passing in a possible_data object or
  "import" if you want to import mics and run logistic regression

- drug:

  NULL if using possible_data, if importing data should be the name of
  the column of the mics or a vector of the mic column and the sign
  colomn

- date_col:

  string, what is the name of the column in the data that corresponds to
  time of sampling

- date_type:

  string, either "decimal", "date", or "year" use decimal if using t
  from possible data, date or year if importing mic data and the date
  column is a date or just a year respectively

- first_year:

  NULL if date_type is "decimal", otherwise a numeric year or decimal
  year value if using "year" or "date" for date_type respectively

- s_breakpoint:

  string, the breakpoint on the MIC scale for what constitutes a
  susceptible isolate, e.g. \<=8 (ug/mL, do not include units)

- r_breakpoint:

  string, the breakpoint on the MIC scale for what constitutes a
  resistant isolate, e.g. \>=128 (ug/mL, do not include units)

- ecoff:

  string or numeric, see plot_fm()

- visual_split:

  string or numeric, see plot_fm()

- k:

  variable passed into logistic regression GAM
