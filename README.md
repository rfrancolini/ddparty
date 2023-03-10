ddparty
================

## DDPARTY

This is for managing and understanding your ddpcr data from QuantaLife.

## Requirements

- [R v4+](https://www.r-project.org/)
- [dplyr](https://CRAN.R-project.org/package=dplyr)
- [readr](https://CRAN.R-project.org/package=readr)
- [stringr](https://CRAN.R-project.org/package=stringr)
- [ggplot2](https://CRAN.R-project.org/package=ggplot2)
- [tidyverse](https://CRAN.R-project.org/package=tidyverse)

## Installation

    remotes::install_github("rfrancolini/ddparty")

## Function: read_ddpcr()

This function will take in csv data from the QuantaLife program and
output it as a usable csv. When downloading csv data from QuantaLife,
the automatically generated csv is formatting differently than the csv
you will obtain if you manually call samples and export it afterwards.
This function allows you to read in either type of csv, add sample
metadata to it, and it will output it as a standard csv with gene copy
number calculated.

### Read Auto Example Data

“Auto” data in this instance is the automatically outputted csv file
when all sample wells have been called automatically. The variable
“calltype” should be set to “auto”.

``` r
library(ddparty)
library(readr)

AutoData <- RawData_example()
SampleMetadata <- metadata_example()

RawAutoData <- read_csv(AutoData)
```

    ## Rows: 26 Columns: 63
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (9): Well, ExptType, Experiment, Sample, TargetType, Target, Status, Co...
    ## dbl (14): CopiesPer20uLWell, PoissonConfMax, PoissonConfMin, Positives, Nega...
    ## lgl (40): TotalConfMax, TotalConfMin, Ch1+Ch2+, Ch1+Ch2-, Ch1-Ch2+, Ch1-Ch2-...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#This is what the raw csv output looks like
head(RawAutoData)
```

    ## # A tibble: 6 × 63
    ##   Well  ExptType    Exper…¹ Sample Targe…² Target Status Conce…³ Super…⁴ Copie…⁵
    ##   <chr> <chr>       <chr>   <chr>  <chr>   <chr>  <chr>  <chr>   <chr>     <dbl>
    ## 1 A05   Absolute Q… ABS     042_2… Ch1Unk… Mem    OK     2       ddPCR …      40
    ## 2 A08   Absolute Q… ABS     075_2… Ch1Unk… Mem    OK     12.7    ddPCR …     254
    ## 3 B08   Absolute Q… ABS     076_2… Ch1Unk… Mem    OK     21.8    ddPCR …     436
    ## 4 C08   Absolute Q… ABS     077_2… Ch1Unk… Mem    OK     2.2     ddPCR …      44
    ## 5 D08   Absolute Q… ABS     078_2… Ch1Unk… Mem    CHECK  No Call ddPCR …      NA
    ## 6 G04   Absolute Q… ABS     040_2… Ch1Unk… Mem    OK     1.1     ddPCR …      22
    ## # … with 53 more variables: TotalConfMax <lgl>, TotalConfMin <lgl>,
    ## #   PoissonConfMax <dbl>, PoissonConfMin <dbl>, Positives <dbl>,
    ## #   Negatives <dbl>, `Ch1+Ch2+` <lgl>, `Ch1+Ch2-` <lgl>, `Ch1-Ch2+` <lgl>,
    ## #   `Ch1-Ch2-` <lgl>, Linkage <lgl>, AcceptedDroplets <dbl>, CNV <lgl>,
    ## #   TotalCNVMax <lgl>, TotalCNVMin <lgl>, PoissonCNVMax <lgl>,
    ## #   PoissonCNVMin <lgl>, ReferenceCopies <lgl>, UnknownCopies <lgl>,
    ## #   Ratio <lgl>, TotalRatioMax <lgl>, TotalRatioMin <lgl>, …

``` r
library(ddparty)

AutoData <- RawData_example()
SampleMetadata <- metadata_example()

AutoExample <- read_ddpcr(filename = AutoData, metadata = SampleMetadata, target = "Mem", calltype = "auto")

#This is our processed data with our metadata attached
head(AutoExample)
```

    ## # A tibble: 6 × 28
    ##   Well  Sample    Status Conce…¹ Copie…² Poiss…³ Poiss…⁴ Posit…⁵ Negat…⁶ Accep…⁷
    ##   <chr> <chr>     <chr>  <chr>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ## 1 A05   042_21SP… OK     2            40     2.8     1.4      35   20199   20234
    ## 2 A08   075_21SM… OK     12.7        254    14.5    10.9     191   17547   17738
    ## 3 B08   076_21SM… OK     21.8        436    24.1    19.4     322   17245   17567
    ## 4 C08   077_21SM… OK     2.2          44     3       1.6      37   19606   19643
    ## 5 D08   078_21SM… CHECK  No Call      NA    NA      NA         0   20386   20386
    ## 6 G04   040_21SP… OK     1.1          22     1.6     0.6      17   18917   18934
    ## # … with 18 more variables: Threshold <dbl>, MeanAmplitudeofPositives <dbl>,
    ## #   MeanAmplitudeofNegatives <dbl>, MeanAmplitudeTotal <dbl>,
    ## #   PoissonConfMax68 <dbl>, PoissonConfMin68 <dbl>, ConcentrationNum <dbl>,
    ## #   CalcCopies <dbl>, year <dbl>, month <dbl>, day <dbl>, season <chr>,
    ## #   site <chr>, site_code <chr>, sample_name <chr>, sample_letter <chr>,
    ## #   tube_number <dbl>, `concentration_(ng/ul)` <dbl>, and abbreviated variable
    ## #   names ¹​Concentration, ²​CopiesPer20uLWell, ³​PoissonConfMax, …

### Read Manual Example Data

“Manual” data in this instance is the csv file that is outputted when at
least one sample well has been manually called, and the csv has been
exported from the QuantaLife GUI. The variable “calltype” should be set
to “manual”.

``` r
library(ddparty)
library(readr)

ManData <- ManuallyEditedData_example()
SampleMetadata <- metadata_example()

ManAutoData <- read_csv(ManData)
```

    ## Rows: 41 Columns: 65
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (10): Well, Sample, Target, Conc(copies/µL), Status, Experiment, SampleT...
    ## dbl (16): Copies/20µLWell, PoissonConfMax, PoissonConfMin, Accepted Droplets...
    ## lgl (39): TotalConfMax, TotalConfMin, Linkage, CNV, TotalCNVMax, TotalCNVMin...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#This is what the raw csv output looks like (note it is different than the "Auto" csv)
head(ManAutoData)
```

    ## # A tibble: 6 × 65
    ##   Well  Sample     Target Conc(…¹ Status Exper…² Sampl…³ Targe…⁴ Super…⁵ DyeNa…⁶
    ##   <chr> <chr>      <chr>  <chr>   <chr>  <chr>   <chr>   <chr>   <chr>   <chr>  
    ## 1 A05   042_21SP_… Mem    2.0367… OK     DQ      Unknown Unknown ddPCR … FAM    
    ## 2 A08   075_21SM_… Mem    12.736… OK     DQ      Unknown Unknown ddPCR … FAM    
    ## 3 B08   076_21SM_… Mem    21.764… OK     DQ      Unknown Unknown ddPCR … FAM    
    ## 4 C08   077_21SM_… Mem    2.2181… OK     DQ      Unknown Unknown ddPCR … FAM    
    ## 5 D08   078_21SM_… Mem    No Call CHECK  DQ      Unknown Unknown ddPCR … FAM    
    ## 6 G04   040_21SP_… Mem    1.0567… OK     DQ      Unknown Unknown ddPCR … FAM    
    ## # … with 55 more variables: `Copies/20µLWell` <dbl>, TotalConfMax <lgl>,
    ## #   TotalConfMin <lgl>, PoissonConfMax <dbl>, PoissonConfMin <dbl>,
    ## #   `Accepted Droplets` <dbl>, Positives <dbl>, Negatives <dbl>,
    ## #   `Ch1+Ch2+` <dbl>, `Ch1+Ch2-` <dbl>, `Ch1-Ch2+` <dbl>, `Ch1-Ch2-` <dbl>,
    ## #   Linkage <lgl>, CNV <lgl>, TotalCNVMax <lgl>, TotalCNVMin <lgl>,
    ## #   PoissonCNVMax <lgl>, PoissonCNVMin <lgl>, ReferenceCopies <lgl>,
    ## #   UnknownCopies <lgl>, Threshold1 <dbl>, Threshold2 <lgl>, …

``` r
library(ddparty)

ManData <- ManuallyEditedData_example()
SampleMetadata <- metadata_example()

ManExample <- read_ddpcr(filename = ManData, metadata = SampleMetadata, target = "Mem", calltype = "manual")

#This is our processed data with our metadata attached
head(ManExample)
```

    ## # A tibble: 6 × 28
    ##   Well  Sample    Status Conce…¹ Copie…² Poiss…³ Poiss…⁴ Posit…⁵ Negat…⁶ Accep…⁷
    ##   <chr> <chr>     <chr>  <chr>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ## 1 A05   042_21SP… OK     2.0367…    40.7    2.79   1.43       35   20199   20234
    ## 2 A08   075_21SM… OK     12.736…   255.    14.5   10.9       191   17547   17738
    ## 3 B08   076_21SM… OK     21.764…   435.    24.1   19.4       322   17245   17567
    ## 4 C08   077_21SM… OK     2.2181…    44.4    3.01   1.58       37   19606   19643
    ## 5 D08   078_21SM… CHECK  No Call    NA     NA     NA           0   20386   20386
    ## 6 G04   040_21SP… OK     1.0567…    21.1    1.65   0.629      17   18917   18934
    ## # … with 18 more variables: Threshold <dbl>, MeanAmplitudeofPositives <dbl>,
    ## #   MeanAmplitudeofNegatives <dbl>, MeanAmplitudeTotal <dbl>,
    ## #   PoissonConfMax68 <dbl>, PoissonConfMin68 <dbl>, ConcentrationNum <dbl>,
    ## #   CalcCopies <dbl>, year <dbl>, month <dbl>, day <dbl>, season <chr>,
    ## #   site <chr>, site_code <chr>, sample_name <chr>, sample_letter <chr>,
    ## #   tube_number <dbl>, `concentration_(ng/ul)` <dbl>, and abbreviated variable
    ## #   names ¹​Concentration, ²​CopiesPer20uLWell, ³​PoissonConfMax, …
