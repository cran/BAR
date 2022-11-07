\name{next_allocation_rate_BAR}
\alias{next_allocation_rate_BAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the allocation probability for the next block of new patients using Bayesian adaptive randomization
}
\description{
Calculate updated allocation probability for each arm based on the accumulative data with binary outcomes
}
\usage{
next_allocation_rate_BAR(n, success_count, tot_num,
                         power_c = "n/2N",
                         lower_bound = .05,
                         control_arm = "",
                         seed = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
the number of patients enrolled for each arm
  }
  \item{success_count}{
the number of responders for each arm
  }
  \item{tot_num}{
the total number of patients enrolled for the trial. If this number cannot be pre-planned, the user can choose argument "power_c" to be numeric instead of "n/2N". In this case, even if the "tot_num" is given a number, this number will not be used
  }
  \item{power_c}{
the power correction of allocation probability. The default value is power_c = "n/2N" and can also be numeric, e.g., power_c = .5
  }
  \item{lower_bound}{
the lower bound of the allocation probability. It must between 0 and \eqn{\frac{1}{K}}. The default value is lower_bound = .05; K indicates total number of arms (including control arm)
  }
  \item{control_arm}{
if this argument is "fixed", then allocation probability of control arm (the first slot)
will be fixed to \eqn{\frac{1}{K}}. The default of this argument will return unfixed results; K indicates total number of arms (including control arm)
  }
  \item{seed}{
the seed. The default value is seed = 100
  }

}
\details{
We show how the updated allocation probabilities for each arm are calculated.\cr \cr Treatments are denoted by \eqn{k = 1,\ldots,K}.\eqn{\;N} is the total sample size. If no burn-in(s), the BAR will be initiated start of a study, that is, for each enrolled patient, \eqn{n = 1,\ldots,N,} the BAR will be used to assign each patient. Denoting the true unknown response rates of \eqn{K} treatments by \eqn{\pi_{1},\ldots,\pi_{K},\;}we can compute \eqn{K} posterior probabilities: \eqn{r_{k,n} = Pr(\pi_{k} = max\{\pi_{1},\ldots,\pi_{K}\}\;|\;Data_{n})}, here, \eqn{n} refers to the \eqn{n}-th patient and \eqn{k} refers to the \eqn{k}-th arm. We calculate the updated probabilities of the BAR algorithm according to the following steps.\cr \eqn{\;}\cr Step 1: (Normalization) Normalize \eqn{r_{k,n}} as \eqn{r_{k,n}^{(c)} = \frac{(r_{k,n})^{c}}{\sum_{j=1}^{K}(r_{j,n})^{c}}}, here \eqn{\;c = \frac{n}{2N}}.\cr \cr
Step 2: (Restriction) To avoid the BAR sticking to very low/high probabilities, a restriction rule to the posterior probability \eqn{r_{k,n}^{(c)}} will be applied: \deqn{Lower\;Bound \le r_{k,n}^{(c)} \le 1 - (K - 1) \times Lower\;Bound,}
\deqn{0 \le Lower\;Bound \le \frac{1}{K}}
After restriction, the posterior probability is denoted as \eqn{r_{k,n}^{(c,re)}}.\cr \cr
Step 3: (Re-normalization) Then, we can have the updated allocation probabilities by the BAR denoted as:\deqn{r_{k,n}^{(f)} = \frac{r_{k,n}^{(c,re)}\times(\frac{r_{k,n}^{(c,re)}}{\frac{n_{k}}{n}})^{2}}{\sum_{j=1}^{K}\{r_{j,n}^{(c,re)}\times(\frac{r_{j,n}^{(c,re)}}{\frac{n_{j}}{n}})^{2}\}}} where \eqn{n_{k}} is the number of patients enrolled on arm \eqn{k} up-to-now.\cr \cr
Step 4: (Re-restriction) Finally, restricts again by using \deqn{Lower\;Bound \le r_{k,n}^{(f)} \le 1 - (K-1) \times Lower\;Bound,}
\deqn{0 \le Lower\;Bound \le \frac{1}{K}}and denote \eqn{r_{k,n}^{(ff)}} as the allocation probability used in the BAR package.
}
\value{
\code{next_allocation_rate_BAR()} returns the updated allocation probability for each arm
}
\references{
Wathen JK, Thall PF. A simulation study of outcome adaptive randomization in multi-arm clinical trials. Clin Trials. 2017 Oct; 14(5): 432-440. doi: 10.1177/1740774517692302.

Xiao, Y., Liu, Z. & Hu, F. Bayesian doubly adaptive randomization in clinical trials. Sci. China Math. 60, 2503-2514 (2017). doi: 10.1007/s11425-016-0056-1.

Hu F, Zhang L X. Asymptotic properties of doubly adaptive biased coin designs for multi-treatment clinical trials. Ann Statist, 2004, 30: 268–301.
}
\author{
Chia-Wei Hsu, Haitao Pan
}
\examples{
## power_c = "n/2N"
next_allocation_rate_BAR(n = c(30, 30, 30),
                         success_count = c(5, 6, 12),
                         tot_num = 150)

## power_c = .5
next_allocation_rate_BAR(n = c(30, 30, 30),
                         success_count = c(5, 6, 12),
                         tot_num = 150, power_c = .5)
}
