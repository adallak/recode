#' @title The global bank stock return volatility data
#'
#' @description A subset of data from the Demirer, M., F. X. Diebold, L. Liu, and K. Yilmaz. 2018. Estimating global bank
#'  network connectedness. Journal of Applied Econometrics 33(1).
#' Original data contains 96 banks from 29 developed and emerging economies (countries) from
#' September 12, 2003, to February 7, 2014. The subset contains only
#' economies where the number of banks in each economy is greater than 4, total of 54
#' banks.
#' Below we provide details on bank names and the country they belong:
#' Variable       Bank        Country
#' hsba.ln        HSBC        Great Britain (gb)
#' X8306.to       MTBH        Japan (jp)
#' jpm            JPM         USA (us)
#' barc.ln        BARC        Great Britain (gb)
#' bac            BAC         USA (us)
#' c              Citi        USA (us)
#' X8411.to       MZH         Japan (jp)
#' rbs.ln         RBS         Great Britain (gb)
#' X8316.to       SMTM        Japan (jp)
#' wfc            WFC         USA (us)
#' lloy.ln        LLOY        Great Britain (gb)
#' ucg.mi         UCG         Italy (it)
#' gs             GS          USA (us)
#' isp.mi         ISP         Italy (it)
#' MS             MS          USA (us)
#' td.t           TD          Canada (ca)
#' ry.t           RBC         Canada (ca)
#' nab.au         NAB         Australia (au)
#' bns.t          BNS         Canada (ca)
#' cba.au         CBA         Australia (au)
#' stan.ln        STAN        Great Britain (gb)
#' X600036.sh     CMB         China (cn)
#' anz.au         ANZ         Australia (au)
#' wbc.au         WBC         Australia (au)
#' X600000.sh     SHGO        China (cn)
#' X600016.sh     CMB         China (cn)
#' bmo.t          BMO         Canada (ca)
#' X8308.to       RSNH        Japan (jp)
#' X8604.to       NMRH        Japan (jp)
#' X8309.to       SMTM        Japan (jp)
#' cm.t           CBC         Canada (ca)
#' bk.us          BNYM        USA (us)
#' usb            USB         USA (us)
#' pnc.us         PNC         USA (us)
#' X000001.sz     PAB         China (cn)
#' cof            COF         USA (us)
#' X600015.sh     HXB         China (cn)
#' bmps.mi        BMPS        Italy (it)
#' stt.us         SSC         USA (us)
#' BBT            BBT         USA (us)
#' NA.T           NBC         Canada (ca)
#' STI.US         STB         USA (us)
#' bp.mi          BP          Italy (it)
#' mqg.au         MQG         Australia (au)
#' X8354.to       FFG         Japan (jp)
#' X8332.to       BOY         Japan (jp)
#' fitb.us        FTB         USA (us)
#' rf.us          RF          USA (us)
#' X8331.to       CBB         Japan (jp)
#' uni.mi         UGF         Italy (it)
#' X8377.to       HFG         Japan (jp)
#' X8355.to       ShB         Japan (jp)
#' mb.mi          MBCF        Italy (it)
#' X8418.to       YFG         Japan (jp)
#' @docType data
#'
#' @usage data(globalbank)
#'
#' @format ## `globalbank`
#' A data frame with 2672 rows and 54 columns:
#' @source <Demiret et. al. (2018)>
"globalbank"



#' @title Air Pollution Dataset
#'
#' @description Time series dataset from Dahlhaus and Eichler (2003) and
#' Davis et al. (2016). Recorded variables include `CO` and `NO`
#' (pollutants mainly emitted from the cars and industry),
#' `NO` and `O` (generated from different reactions in the atmosphere),
#' and the global solar radiation intensity `R`.
#'
#' @docType data
#'
#' @usage data(airpollute)
#'
#' @format ## `airpollute`
#' A data frame with 8760 rows and 5 columns:
#' @source <Dahlhaus and Eichler (2003)>
"airpollute"


#' @title Oslo 2 data
#'
#' @description Includes n = 280 breast cancer samples collected from hospitals in Oslo
#'  The gene expression for p=100 genes has been quantified by measuring their proteins using RPPA technology.
#'
#' @docType data
#'
#' @usage data(oslo2rppa)
#'
#' @format ## `oslo2rppa`
#' A data frame with 280 rows and 100 columns:
#' @source <Lingjarde, C., Lien, T.G., Borgan, O. et al.
#'  Tailored graphical lasso for data integration in gene network reconstruction.
#'  BMC Bioinformatics 22, 498 (2021). https://doi.org/10.1186/s12859-021-04413-z>
"oslo2rppa"


#' @title Dataset from cardiovascular microarray study
#'
#' @description The dataset has 63 subjects with 44 healthy controls, and
#' 19 cardiovascular patients, and 20436 genes measured for
#' each subject.
#'
#' @docType data
#'
#' @format ## `cardiodata`
#' A data frame with 20436 rows and 63 columns:
#' @source <Efron. B. Large-Scale InferenceL Empirical Bayes
#' Methods for Estimation, Testing, and Prediction (2010)>
"cardiodata"

#' @title Percentage of Body Fat and Body Measurements
#'
#' @description Age, weight, height, and 10 body circumference measurements
#' are recorded for 252 men. Each man's percentage of body fat was
#' accurately estimated by an underwater weighing technique.
#'
#' @docType data
#'
#' @format ## `fat`
#' A data frame with 252 observations on the following 18 variables:
#' @source <Johnson R. Journal of Statistics Education v.4, n.1 (1996)>
"fat"

#' Information related to the COVID-19 infections for each
#' date-region pair
#'
#' @docType data
#'
#' @format ## `covid`
#' Covid data
#' @source <https://github.com/GoogleCloudPlatform/covid-19-open-data/blob/main/docs/table-epidemiology.md>
"covid"

#' @title Log returns of 8 stocks
#'
#' @description Log-returns of eight US stocks (Walmart, Exxon, Ford. GE, Conoco Phillips,Citigroup,
#' IBM, and AIG) from 2017 to 2022. The dataset consists of 314 observations
#'
#' @docType data
#'
#' @usage data(stocks)
#'
#' @format ## `stocks`
#' A data frame with 314 rows and 8 columns:
#' @source YahooFinance
"stocks"


#' @title Car brand perception survey.
#' @description The data is based on a survey of customers
#' from Nissan's target segments, described as people between 25 and 35 with annual household incomes between
#' 50000 and 100000. Participants were asked about their perceptions (from scale 1 to 7) of 10 different cars
#' (Infinity G20, Ford T-bird, Audi 90, Toyota Supra, Eagle Talon, Honda Prelude, Saab 900, Pontiac Firebird,
#' BMW 318i, Mercury Capri) on 15 attributes
#'
#' @docType data
#' @usage data(brand)
#'
#' @format ## A data frame with 1000 rows and 16 variables:
"brand"


#' @title  Stock and Watson (2002) dataset.
#' @description A dataset containing 215 monthly US macroeconomic variables spanning the period 1959-1998.
#' The dataset is borrowed from  Stock and Watson (2002).
#'
#' @docType data
#' @usage data(SW2002)
#'
#' @format ## A data frame with 215 rows and 146 variables:
#' @source http://www.princeton.edu/~mwatson/ddisk/wc00.zip
"SW2002"


#' @title  Part of Million Song dataset.
#' @description A dataset consists of information about 10000 songs, which are mostly western and range from
#' 1922 to 2011. It contains p = 90 timbre features.
#' The dataset is borrowed from  https://archive.ics.uci.edu/dataset/203/yearpredictionmsd.
#'
#' @docType data
#' @usage data(yearpredictionmsd)
#'
#' @format ## A data frame with 10000 rows and 91 variables:
#' @source https://archive.ics.uci.edu/dataset/203/yearpredictionmsd
"yearpredictionmsd"
