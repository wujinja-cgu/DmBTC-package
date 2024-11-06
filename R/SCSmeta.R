#' Run the Split Component Synthesis (SCS) Model
#'
#' @param tp The number of true-positive.
#' @param fp The number of false-positive.
#' @param fn The number of false-negative.
#' @param tn The number of true-negative.
#' @return A list containing the results of the model fitting.
#' @export

SCSmeta=function(tp,fp,fn,tn){

  require(arm)

  rm()

  #data correction
  studies=length(tp)

  tp_c=tp
  fp_c=fp
  fn_c=fn
  tn_c=tn

  continuity = rep(0, studies)

  for (i in 1: studies){
    if (tp[i]==0|fp[i]==0|fn[i]==0|tn[i]==0) continuity[i]=1}

  for (i in 1: studies){  if (continuity[i]==1)    tp_c[i]=tp[i]+0.5}
  for (i in 1: studies){  if (continuity[i]==1)    fp_c[i]=fp[i]+0.5}
  for (i in 1: studies){  if (continuity[i]==1)    fn_c[i]=fn[i]+0.5}
  for (i in 1: studies){  if (continuity[i]==1)    tn_c[i]=tn[i]+0.5}

  #the ivhet estimator

  n_var=1/tp_c+1/fn_c+1/tn_c+1/fp_c
  dor_rand=(tp_c*tn_c)/(fp_c*fn_c)
  l_dor<-log(dor_rand)

  w=1/n_var
  l_dor_ihvet=sum(w*l_dor)/sum(w)

  Q=sum(w*(l_dor-l_dor_ihvet)^2)

  nom=Q-(studies-1)
  denom=sum(w)-sum(w^2)/sum(w)
  bst_h=nom/denom
  bst=max(0,bst_h)

  var_ihv=sum((w/sum(w))^2*(n_var+bst))
  sd_ihv=sqrt(var_ihv)

  l_low_ihvet=l_dor_ihvet-1.96*sd_ihv
  l_up_ihvet=l_dor_ihvet+1.96*sd_ihv

  dor_ihvet=exp(l_dor_ihvet)
  dor_ihvet_lo=exp(l_low_ihvet)
  dor_ihvet_up=exp(l_up_ihvet)

  #pooled Se, Sp estimation

  centerlogDOR=l_dor-l_dor_ihvet

  logitSe <- lm(logit(tp_c/(tp_c+fn_c)) ~ centerlogDOR)
  logitSp <- lm(logit(tn_c/(tn_c+fp_c)) ~ centerlogDOR)

  Se_ihvet=invlogit(as.numeric(logitSe$coefficients[1]))
  Sp_ihvet=invlogit(as.numeric(logitSp$coefficients[1]))

  var_Se_ihvet=var_ihv*(0.5-((Se_ihvet*(1-Se_ihvet))-0.25)+((Sp_ihvet*(1-Sp_ihvet))-0.25))
  var_Sp_ihvet=var_ihv*(0.5+((Se_ihvet*(1-Se_ihvet))-0.25)-((Sp_ihvet*(1-Sp_ihvet))-0.25))

  sd_se=sqrt(var_Se_ihvet)
  sd_sp=sqrt(var_Sp_ihvet)

  l_low_Se_ihvet=as.numeric(logitSe$coefficients[1])-1.96*sd_se
  l_up_Se_ihvet=as.numeric(logitSe$coefficients[1])+1.96*sd_se

  Se_ihvet_l=invlogit(l_low_Se_ihvet)
  Se_ihvet_u=invlogit(l_up_Se_ihvet)

  l_low_Sp_ihvet=as.numeric(logitSp$coefficients[1])-1.96*sd_sp
  l_up_Sp_ihvet=as.numeric(logitSp$coefficients[1])+1.96*sd_sp

  Sp_ihvet_l=invlogit(l_low_Sp_ihvet)
  Sp_ihvet_u=invlogit(l_up_Sp_ihvet)

  #AUC estimation

  if (dor_ihvet>=1) AUC=invlogit(l_dor_ihvet/2) else AUC=1-invlogit(log(1/dor_ihvet)/2)
  if (dor_ihvet>=1) AUC_l=invlogit(l_low_ihvet/2) else AUC_l=1-invlogit(log(1/dor_ihvet_lo)/2)
  if (dor_ihvet>=1) AUC_u=invlogit(l_up_ihvet/2) else AUC_u=1-invlogit(log(1/dor_ihvet_up)/2)

  #pLR & nLR

  pLR=Se_ihvet/(1-Sp_ihvet)
  ln_pLR = log(pLR)

  nLR=(1-Se_ihvet)/Sp_ihvet
  ln_nLR = log(nLR)

  prop1 = abs(ln_pLR) / (abs(ln_pLR)+abs(ln_nLR))
  prop2 = abs(ln_nLR) / (abs(ln_pLR)+abs(ln_nLR))

  se_plr = sqrt(prop1*var_ihv)
  if (dor_ihvet==1) se_plr = sqrt(0.5*var_ihv)

  pLR_l = exp(ln_pLR-1.96*se_plr)
  pLR_u = exp(ln_pLR+1.96*se_plr)

  se_nlr = sqrt(prop2*var_ihv)
  if (dor_ihvet==1) se_nlr = sqrt(0.5*var_ihv)

  nLR_l = exp(ln_nLR-1.96*se_nlr)
  nLR_u = exp(ln_nLR+1.96*se_nlr)

  #ROC curve
  sp_roc=seq(0,1,0.01)
  se_roc=(dor_ihvet*(1-sp_roc))/((dor_ihvet*(1-sp_roc))+sp_roc)
  se_roc_l=(dor_ihvet_lo*(1-sp_roc))/((dor_ihvet_lo*(1-sp_roc))+sp_roc)
  se_roc_up=(dor_ihvet_up*(1-sp_roc))/((dor_ihvet_up*(1-sp_roc))+sp_roc)

  fp_roc=1-sp_roc
  #plot(se_roc~fp_roc, type="l", col="red", xlab="1-Specificity", ylab="Sensitivity")
  #text(0.79,0.0, "Weights are from Doi's IVhet model", cex=0.7, font=3)
  #lines(se_roc_l~fp_roc, lty=3)
  #lines(se_roc_up~fp_roc, lty=3)
  #abline(0, 1, lty=2)

  #points((fp/(tn+fp)),(tp/(tp+fn)),cex=w/5)
  #points((1-Sp_ihvet),Se_ihvet, pch=18,col="blue", cex=2)

  #print(paste("The DOR is",round(dor_ihvet,3),"(",round(dor_ihvet_lo,3),"; ",round(dor_ihvet_up,3),")"))
  #print(paste("The AUC is",round(AUC,3),"(",round(AUC_l,3),"; ",round(AUC_u,3),")"))
  #print(paste("The Se is",round(Se_ihvet,3),"(",round(Se_ihvet_l,3),"; ",round(Se_ihvet_u,3),")"))
  #print(paste("The Sp is",round(Sp_ihvet,3),"(",round(Sp_ihvet_l,3),"; ",round(Sp_ihvet_u,3),")"))
  #print(paste("The pLR is",round(pLR,3),"(",round(pLR_l,3),"; ",round(pLR_u,3),")"))
  #print(paste("The nLR is",round(nLR,3),"(",round(nLR_l,3),"; ",round(nLR_u,3),")"))
  return(c(Se_ihvet,Sp_ihvet))

}
