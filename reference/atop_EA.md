# East Asian consultation agreements from ATOP

Consultation alliance agreements involving East and Central Asian
countries from 2010 to 2018. Sliced from the ATOP dyad-year dataset
(atop5_1dy). The dataset includes 18 countries: an East Asian core (CHN,
IND, JPN, KOR, PRK, PHL, THA, VNM) and their consultation treaty
partners (ARM, KAZ, KGZ, PAK, RUS, TJK, TKM, UKR, USA, UZB).

## Usage

``` r
data(atop_EA)
```

## Format

A data frame with 156 rows and 3 columns:

- year:

  Year of consultation agreement (2010-2018)

- country1, country2:

  3-letter ISO country codes for agreement partners

## Source

[Alliance Treaty Obligations and Provisions (ATOP) dyad-year
dataset](http://www.atopdata.org)

## References

Leeds, Brett Ashley, Jeffrey M. Ritter, Sara McLaughlin Mitchell, and
Andrew G. Long. 2002. Alliance Treaty Obligations and Provisions,
1815-1944. International Interactions 28: 237-260.
([ATOP](http://www.atopdata.org))

## Examples

``` r
data(atop_EA)
head(atop_EA)
#> # A tibble: 6 × 3
#>    year country1 country2
#>   <dbl> <chr>    <chr>   
#> 1  2010 USA      KOR     
#> 2  2011 USA      KOR     
#> 3  2012 USA      KOR     
#> 4  2013 USA      KOR     
#> 5  2014 USA      KOR     
#> 6  2015 USA      KOR     
```
