#' Simulation the data for DmBTC Model
#'
#' @param mu1 the mean value of non-diseased individuals.
#' @param mu2 the mean value of diseased individuals, mu2 greater than mu1.
#' @param DOR Magnitude of diagnostic accuracy.
#' @param prevalence Prevalence of disease.
#' @param nstudy Number of studies.
#' @param is_symmetric Symmetry or asymmetry in the summary receiver operating characteristics (SROC) curve.
#' @return A list containing the results of the simulations.
#' @export

DmBTC_Simulation<-function(mu1,
                          mu2,
                          DOR,
                          prevalence,
                          nstudy,
                          is_symmetric,
                          seed){

  ## generate data
  rm()

  seed=seed

  mu1=mu1

  mu2=mu2

  DOR=DOR

  sigma1=((mu2-mu1)*0.75*3.14159265358979)/(sqrt(3)*log(DOR))

  sigma2=is_symmetric*sigma1

  prevalence=prevalence

  threshold=(mu1+mu2)/2

  nstudys=nstudy

  sensitivity=exp((sqrt((3.14)^2/3))*((mu2-threshold)/sigma2))/(1+exp(sqrt((3.14)^2/3)*((mu2-threshold)/sigma2)))

  sensitivity

  specificity=1-(exp((sqrt((3.14)^2/3))*((mu1-threshold)/sigma1))/(1+exp(sqrt((3.14)^2/3)*((mu1-threshold)/sigma1))))

  specificity

  ### Sampling data
  samplesize=sample(seq(20,300,1),nstudys)

  ndisease=array(0,nstudys)

  TP=array(0,nstudys)

  FN=array(0,nstudys)

  TN=array(0,nstudys)

  FP=array(0,nstudys)

  for(i in 1:nstudys){

          ndisease[i]=as.integer(samplesize[i]*prevalence)

          test_disease=array(0,ndisease[i])

          test_disease_binary=array(0,ndisease[i])

          test_nondisease=array(0,samplesize[i]-ndisease[i])

          test_nondisease_binary=array(0,samplesize[i]-ndisease[i])

          for(j in 1:ndisease[i]){

            test_disease[j]=rnorm(1,mu2,sigma2)

            test_disease_binary[j]=ifelse(test_disease[j]>threshold,1,0)

            }

          TP[i]=sum(test_disease_binary)

          FN[i]=ndisease[i]-TP[i]

          for(k in 1:samplesize[i]-ndisease[i]){

            test_nondisease[k]=rnorm(1,mu1,sigma1)

            test_nondisease_binary[k]=ifelse(test_nondisease[k]<threshold,1,0)

            }

          TN[i]=sum(test_nondisease_binary)

          FP[i]=samplesize[i]-ndisease[i]-TN[i]

          }

  ID=rep(1:nstudys)

  metadata=cbind(ID,TP,FN,TN,FP)

  metadata=as.data.frame(metadata)

  return(c(sensitivity,specificity,metadata))
}
