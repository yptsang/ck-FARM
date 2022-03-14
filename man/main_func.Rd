\name{An R Package to Discover Big Data Association for Business Intelligence}

\title{ck-FARM}
\description{
  The proposed algorithm, namely ck-FARM, is deployed in the R programming environment with applying five existing   libraries, including Matrix, arules, ca, ggplots, and factoextra. Thus, the statistical significance between process parameters can be examined, while fuzzy association rules (FARs) can be automatically built. Based on the fuzzy relationship of process parameters, corresponding business management/re-engineering strategies can be established to strengthen business intelligence capabilities.
}

\usage{
corrana_data(data)
mem_fun_result(data)
main_func(data, alpha, lambda)
}

\arguments{
  \item{data}{
    Given dataset for discovering the data association (in either .txt/.csv file)
  }
  \item{alpha}{
    Support count threshold
  }
  \item{lambda}{
    Confidence threshold
  }
  \item{max_clu}{
    Maximum number of clusters for tuning fuzzy membership functions
  }
}
