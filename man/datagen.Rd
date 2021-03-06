\name{datagen}
\alias{datagen}

\title{
Generate epidemics from distance-based and/or network-based continuous-time ILMs.
}
\description{
Generate epidemics from distance-based and/or network-based ILMs.
}
\usage{

datagen(type, kerneltype, kernelmatrix, distancekernel = NULL,

initialepi = NULL, tmax = NULL, suspar = NULL, transpar = NULL,

powersus = NULL, powertrans = NULL, kernel.par = NULL, spark = NULL,

gamma = NULL, delta, suscov = NULL, transcov = NULL)

}

\arguments{
  \item{type}{
 type of compartmental framework, with the choice of ``SIR'' for Susceptible-Infectious-Removed and ``SINR'' for Susceptible-Infectious-Notified-Removed.
}
  \item{kerneltype}{
type of kernel function with choice of ``distance'' for a distance-based ILM, ``network'' for a contact network-based ILM and ``both'' for a combination of network- and distance-based ILM.
}
  \item{kernelmatrix}{
a matrix of the XY coordinates of individuals when \code{kerneltype} is set to ``distance''. An object of class ``contactnet'', or a matrix of incidence (0 or 1) or weights [0,\eqn{\infty}) representing contact network when \code{kerneltype} is set to ``network''. A list of two matrices (1: the XY coordinates of individuals, 2: contact network matrix or an object of class ``contactnet'') when \code{kerneltype} is set to ``both''.
}

  \item{distancekernel}{
the spatial kernel type when \code{kerneltype} is set to ``distance'' or ``both''. Choices are ``powerlaw'' for a power-law distance kernel or ``Cauchy'' for a Cauchy distance kernel.
}
  \item{initialepi}{
a matrix of the initial \eqn{n_{ini}}{n_ini} infected individuals. When \code{type} is set to ``SIR'', it is an \eqn{n_{ini}\times 4}{n_ini by 4} matrix representing  the id numbers of the individuals, removal times, infectious periods, and infection times. When \code{type} is set to ``SINR'', it is an \eqn{n_{ini}\times 6}{n_ini by 6} matrix representing the id numbers of the individuals, removal times, delay periods, notification times, incubation periods, and infection times. If set to the default = NULL, one randomly chosen individual is infected with infection time set to zero.
}
  \item{tmax}{
the time at which the epidemic simulation ends. If not set the simulation is run until the epidemic dies out or all individuals are infected and then enter the removed state.
}
  \item{suspar}{
parameter(s) for the susceptibility function (>0). Default value(s) is 1.
}
  \item{transpar}{
parameter(s) for the transmissibility function (>0). Default value(s) is 1.
}
  \item{powersus}{
the power parameter(s) of the susceptibility function (>0). Default value(s) is 1.
}
  \item{powertrans}{
the power parameter(s) of the transmissibility function (>0). Default value(s) is 1.
}
  \item{kernel.par}{
a scalar spatial parameter when \code{kerneltype} is set to ``distance'' (>0), or a vector of the spatial and network effect parameters when \code{kerneltype} is set to ``both'' (>0). Default = NULL means that it is not required, as used when the \code{kerneltype} is set to ``network''.
}
  \item{spark}{
spark parameter (>=0), representing random infections that are unexplained by other parts
of the model. Default value is 0.
}
  \item{gamma}{
the notification effect parameter for SINR model. The default value is 1.
}
  \item{delta}{
a vector of the shape and rate parameters of the gamma-distributed infectious period (SIR) or a \eqn{2\times2}{2 by 2} matrix of the shape and rate parameters of the gamma-distributed incubation and delay periods (SINR). Note that an exponential distribution can be assigned as it is a special case of the gamma distribution where the shape parameter is equal to 1.
}
  \item{suscov}{
covariate matrix of the susceptibility function, representing each covariate in a column. If intercept exists or in the case of no covariates, the first column is set to unity.
}
  \item{transcov}{
covariate matrix of the transmissibility function, representing each covariate in a column. If intercept exists or in the case of no covariates, the first column is set to unity.
}
}

\details{

We consider the simulation of epidemics from continuous-time ILMs in the following way. Each infected individual has an infection life history defined by their time of infection and the length of time spent in the infectious state. We are assuming the conditional intensity functions stay constant between events, so that the time to the next infected individual, given that the last infection occurred at time \eqn{t}, is given by the minimum \eqn{W_{j} \sim Exp(\lambda_{jt})}{W_j ~ Exp(\lambda_{jt})}, where \eqn{W_{j}}{W_j} represents the ``waiting time'' for susceptible individual \eqn{j} becoming infected.

Under both SIR and SINR ILMs option, an epidemic is simulated using the specified model starting with a randomly chosen infected individual \eqn{k} at time \eqn{I_{1}}{I_1}, or as specified via the \code{initialepi} option. Under the SIR model, the rate of infection of each susceptible individual \eqn{j} at time \eqn{t} \eqn{\in \rm I\!R^{+}}{in [0,\infty)} follows the framework of Deardon et al. (2010):

\if{html}{\figure{SIR-rate.jpg}{options: width="80\%"}}
\if{latex}{
\deqn{\lambda_{jt} = \left[ \Omega_{S}(j) \sum_{i \in \mathcal{I}{(t)}}{\Omega_{T}(i) \kappa(i,j)} \right]+ \epsilon(j,t),}
}

where \eqn{\mathcal{I}(t)}{I(t)} is the set of infectious individual at time \eqn{t}{t}. Under the SINR model, the rate of infection of each susceptible individual \eqn{j}{j} at time \eqn{t}{t} follows Jewell et al. (2009):

\if{html}{\figure{SINR-rate.jpg}{options: width="80\%"}}
\if{latex}{
\deqn{\lambda_{jt} = \left[ \Omega_{S}(j) \sum_{i \in \mathcal{N}^{-}{(t)}}{\Omega_{T}(i) \kappa(i,j)} + \gamma \left(\Omega_{S}(j) \sum_{i \in \mathcal{N}^{+}{(t)}}{\Omega_{T}(i)\kappa(i,j)}\right)\right] + \epsilon(j,t),}
}

where \eqn{\mathcal{N}^{-}{(t)}}{N^{-}(t)} is the set of infectious individual at time \eqn{t} who have been infected but have not reached the notified state; \eqn{\mathcal{N}^{+}{(t)}}{N^{+}(t)} is the set of infectious individual at time \eqn{t} who have been infected and reached the notified state; \eqn{\gamma}{\gamma} is the notification effect parameter usually allowing for infectious risk reduction after notification (Jewell et al., 2009); \eqn{\Omega_{S}(j)} is a function of risk factors associated with susceptible individual \eqn{j}{j} (i.e., susceptibility); and \eqn{\Omega_{T}(i)} is a function of risk factors associated with infectious individual \eqn{i}{i} (i.e., transmissibility). The \eqn{\Omega_{S}(j)} and \eqn{\Omega_{T}(i)} are defined as:

\if{html}{\figure{sus-trans-covariate.jpg}{options: width="80\%"}}
\if{latex}{
\deqn{ \Omega_{S}(j) = \mathbf{S} \mathbf{X}_{j}^{\phi}; \ \ \textrm{and}  \ \ \  \Omega_{T}(i) = \mathbf{T} \mathbf{Z}_{i}^{\xi},}
}

where \eqn{\mathbf{S}}{S} and \eqn{\mathbf{T}}{T} are the parameters vectors of the susceptibility and transmissibility functions, respectively; \eqn{\mathbf{X}}{X} and \eqn{\mathbf{Z}}{Z} are the susceptibility and transmissibility risk factor matrices, respectively; and \eqn{\phi}{\phi} and \eqn{\xi}{\xi} are vectors of the power parameters of the susceptibility and transmissibility functions, respectively. These power parameters allow for nonlinear relationships between the susceptibility and transmissibility risk factors and the infection rate (Deardon et al., 2010).

The infection kernel \eqn{\kappa(i,j)}{k(i,j)} allows for distance (spatial) and/or contact network-based mechanisms. In distance-based models, the kernel function \eqn{\kappa(i,j)}{k(i,j)} is a power-law or Cauchy function of distance. In network-based models, the kernel function \eqn{\kappa(i,j)}{k(i,j)} is the contact network elements (binary (\eqn{c_{ij}}) or weighted (\eqn{w_{ij}})). In both distance-network-based models, the kernel function \eqn{\kappa(i,j)}{k(i,j)} is represented by two terms: (1) spatial function, and (2) the contact network.

The optional kernel function formulae are as follows:


\itemize{
\item \strong{Distance-based ILMs}:
	\itemize{
		\item Power-law spatial kernel:

			\eqn{\kappa(i,j) = d_{ij}^{-\beta}}{k(i,j) = d_{ij}^{-\beta}}
		\item Cauchy spatial kernel:

			\eqn{\kappa(i,j) = \frac{\beta}{d_{ij}^{2}+\beta^{2}}}{k(i,j) = \beta/(d_{ij}^{2}+\beta^{2})}
		}

\item \strong{Network-based ILMs}:
	\itemize{
		\item \eqn{\kappa(i,j) = c_{ij}}{k(i,j) = c_{ij}} for binary contact network or,
		\item \eqn{\kappa(i,j) = w_{ij}}{k(i,j) = w_{ij}} for weighted contact network.
	}
\item \strong{Both distance- and network-based ILMs}:
	\itemize{
		\item With power-law spatial kernel:
		\itemize{
			\item \eqn{\kappa(i,j) = d_{ij}^{-\beta_{1}} +\beta_{2}c_{ij}}{k(i,j) = d_{ij}^{- \beta_{1}} + \beta_{2}c_{ij}} with binary contact network or,
			\item \eqn{\kappa(i,j) = d_{ij}^{-\beta_{1}} +\beta_{2}w_{ij}}{k(i,j) = d_{ij}^{- \beta_{1}} + \beta_{2}w_{ij}} with weighted contact network.
		}
		\item With Cauchy spatial kernel:
		\itemize{
			\item \eqn{\kappa(i,j) = \frac{\beta_{1}}{(d_{ij}^{2}+\beta_{1}^{2})} + \beta_{2}c_{ij}}{k(i,j) = \beta_{1}/(d_{ij}^{2}+\beta_{1}^{2}) + \beta_{2}c_{ij}} with binary contact network or,
			\item \eqn{\kappa(i,j) = \frac{\beta_{1}}{(d_{ij}^{2}+\beta_{1}^{2})} + \beta_{2}w_{ij}}{k(i,j) = \beta_{1}/(d_{ij}^{2}+\beta_{1}^{2}) + \beta_{2}w_{ij}} with weighted contact network.
		}
	}
}

At time \eqn{I_{s}}{I_s}, the waiting time until being infected is drawn for each susceptible individual \eqn{j}{j} from \eqn{W_{j} \sim Exp(\lambda_{jI_{s}})}{W_j ~ Exp(\lambda_{jI_s})}. The individual with minimum \eqn{W}{W} is chosen as the next infected individual, and assigned infection time \eqn{I_{s+1} = I_{s}+min(W_{j})}{I_{s+1} = I_s + min(W_j)}. Under SIR model, this newly-infected individual is also assigned an infectious period \eqn{D_{j}}{D_j} generated from the infectious period (exponential  or gamma) distribution \eqn{f(D_{j};\delta)}{f(D_j;\delta)}. Thus, its removal time becomes \eqn{R_{s+1} = I_{s+1}+D_{j}}{R_{s+1} = I_{s+1}+D_j}. Under SINR model, each infected individual is assigned an incubation period defining the time from infection to notification, and a delay period defining the time from notification to removal. Thus, this newly-infected individual is assigned an incubation and delay period \eqn{D^{(inc)}_{j}}{D^{(inc)}_j} and \eqn{D^{(delay)}_{j}}{D^{(delay)}_j}, generated from the (exponential  or gamma) incubation and delay periods distributions \eqn{f(D^{(inc)}_{j};\delta^{(inc)})}{f(D^{(inc)}_j;\delta^{(inc)})} and \eqn{f(D^{(delay)}_{j};\delta^{(delay)})}{f(D^{(delay)}_j;\delta^{(delay)})}, respectively. Thus, its notification and removal times become \eqn{N_{s+1} = I_{s+1}+D^{(inc)}_{j}}{N_{s+1} = I_{s+1}+D^{(inc)}_j} and \eqn{R_{s+1} = N_{s+1}+D^{(delay)}_{j}}{R_{s+1} = N_{s+1}+D^{(delay)}_j}, respectively. The process is repeated until no infectives remain in the population, or until \eqn{I_{s+1} >}{I_s+1 > } \code{tmax} if it is specified.

}

\value{
The function \code{datagen} returns objects of class ``datagen'' containing the following:
\if{html}{
	 \describe{
    \item{type:}{The specified type of the compartmental framework.}
    \item{kerneltype:}{The specified type of the kernel function.}
    \item{epidat:}{Under an SIR ILM, it is a matrix with four columns representing  the id numbers of the individuals, removal times, infectious periods, and infection times, respectively. Under an SINR ILM, it is a matrix with six columns: the id numbers of the individuals, removal times, delay periods, notification times, incubation periods, and infection times. Uninfected individuals are assigned infinity values \eqn{(\infty)} for both their removal and infection times. Note that the rows of \code{epidat} is sorted according to the infection times in ascending order. The id numbers of the individuals is used here to respresnt the corresponding row number of \code{location}, contact \code{network}, and covariate matrices.}
    \item{location:}{A matrix of the XY coordinates of individuals if it is given in the function, otherwise, a zero matrix.}
    \item{network:}{The contact network matrix if it is given in the function.}
  }
}
\if{latex}{
	 \describe{
	    \item{type:}{The specified type of the compartmental framework.}
	    \item{kerneltype:}{The specified type of the kernel function.}
    \item{epidat:}{Under an SIR ILM, it is a matrix with four columns representing  the id numbers of the individuals, removal times, infectious periods, and infection times, respectively. Under an SINR ILM, it is a matrix with six columns: the id numbers of the individuals, removal times, delay periods, notification times, incubation periods, and infection times. Uninfected individuals are assigned infinity values \eqn{(\infty)} for both their removal and infection times. Note that the rows of \code{epidat} are sorted according to the infection times in ascending order. The id numbers of the individuals are used here to represent the corresponding row number of the \code{location}, contact \code{network}, and covariate, matrices.}
	    \item{location:}{A matrix of the XY coordinates of individuals if it is given in the function, otherwise, a zero matrix.}
	    \item{network:}{The contact network matrix if it is given in the function, otherwise, a zero matrix.}
	  }
}
\if{text}{
	 \describe{
    \item{type:}{The specified type of the compartmental framework.}
    \item{kerneltype:}{The specified type of the kernel function.}
    \item{epidat:}{Under an SIR ILM, it is a matrix with four columns representing  the id numbers of the individuals, removal times, infectious periods, and infection times, respectively. Under an SINR ILM, it is a matrix with six columns: the id numbers of the individuals, removal times, delay periods, notification times, incubation periods, and infection times. Uninfected individuals are assigned infinity values \eqn{(\infty)} for both their removal and infection times. Note that the rows of \code{epidat} are sorted according to the infection times in ascending order. The id numbers of the individuals are used here to represent the corresponding row number of \code{location}, contact \code{network}, and covariate, matrices.}
    \item{location:}{A matrix of the XY coordinates of individuals if it is given in the function, otherwise, a zero matrix.}
    \item{network:}{The contact network matrix if it is given in the function, otherwise, a zero matrix.}
  }
}
}

\references{
  Jewell, C. P., Kypraios, T., Neal, P., and Roberts, G. O. (2009). Bayesian analysis for emerging infectious diseases. \emph{Bayesian Analysis}, 4(3):465-496.

  Deardon, R., Brooks, S. P., Grenfell, B. T., Keeling, M. J., Tildesley, M. J., Savill, N. J., Shaw, D. J., and Woolhouse, M. E. (2010). Inference for individual-level models of infectious diseases in large populations. \emph{Statistica Sinica}, 20(1):239.
}

\examples{
set.seed(3245)
loc<- matrix(cbind(runif(100, 1, 100),runif(100, 1, 100)), ncol = 2, nrow = 100)
net1<- contactnet(type = "powerlaw", location = loc, beta = 0.9)
cov <- cbind(rep(1, 100), runif(100, 1, 50))
#To simulate the epidemic using the network kernel:
epi.net <- datagen(type = "SIR", kerneltype = "network",
               kernelmatrix = net1, suspar = c(0.08, 0.05),
               delta = c(4, 2), suscov = cov)
#To simulate the epidemic using the distance kernel:
epi.dis <- datagen(type = "SIR", kerneltype = "distance",
               distancekernel = "powerlaw", kernelmatrix = loc,
               suspar = c(0.8, 0.5), kernel.par = 2,
               delta = c(4, 2), suscov = cov)
}

\seealso{
  \code{\link{loglikelihoodepiILM}, \link{contactnet}, \link{plot.datagen}}.
}

\keyword{Simulation}
