
#########################################################################################
#         Functions for performing mediation analysis from Kim YM et al. 2020    ´      #
#########################################################################################

## Causal mediation analysis for a binary mediator
# (1) using cox or clr for mediator model will cause a problem since no intercept (th0). 
#     However, in KIM's study fro tables 12-14, glm was used as mediator model. We followed KIM for out study. 
#     For the R mediation package, cox is also not applied as mediator model.
# (2) when no interaction between exposure and mediator, add 0 to the last position of "beta"

CDMA <- function(beta,th,u){   
	beta1 <- beta[1]; beta2 <- beta[2]; beta3 <- beta[length(beta)] 
	th0 <- th[1]; th1 <- th[2]; th2 <- th[-c(1,2)]
	x1 <- 1; x2 <- 0
	
	# equation (17) for conditional log OR of direct effects ---------
	NDEN <- exp(beta1*x1)*(1+exp(beta2+beta3*x1+th0+th1*x2+sum(th2*u)))  # exposed 
	NDED <- exp(beta1*x2)*(1+exp(beta2+beta3*x2+th0+th1*x2+sum(th2*u)))  # not exposed
	NDE <- log(NDEN/NDED)
	
	# equation (18) for conditional log OR of indirect effects ---------
	NIEN <- (1+exp(th0+th1*x2+sum(th2*u)))*(1+exp(beta2+beta3*x1+th0+th1*x1+sum(th2*u)))
	NIED <- (1+exp(th0+th1*x1+sum(th2*u)))*(1+exp(beta2+beta3*x1+th0+th1*x2+sum(th2*u)))
	NIE <- log(NIEN/NIED)
	
	# equation (19) for MP on the log OR difference scale ---------
	MP <- NIE/(NDE+NIE)
	
	# output direct and indirect effects on OR scales but MP on difference in log ORs scales -------
	list(NIE=exp(NIE),NDE=exp(NDE),MP=MP) 
}


## Variance function for a binary mediator (revised based on KIM YM et al., 2020 and our model)
VAR.D <- function(beta, th, sig.beta, sig.th, u, interaction=T){
	## u : baseline covariates and row vectors
	## beta3: the interaction term between an exposure and a mediator on the outcome
	## sig.bet : variance-covariance matrix of beta estimators
	## sig.th : variance-covariance matrix of theta estimators
	p <- length(beta); q <- length(th)
	if(interaction==T){
		beta1 <- beta[1]; beta2 <- beta[2]; beta3 <- beta[3]; beta4 <- beta[4:p]
	} else {
		beta1 <- beta[1]; beta2 <- beta[2]; beta3 <- 0; beta4 <- beta[3:p]
		sig.beta1 <- matrix(0,nrow=p+1,ncol=p+1)
		sig.beta1[1:2,1:2] <- sig.beta[1:2,1:2]; 
		sig.beta1[1:2,4:(p+1)] <- sig.beta[1:2,3:p]
		sig.beta1[4:(p+1),1:2] <- sig.beta[3:p,1:2]
		sig.beta1[4:(p+1),4:(p+1)] <- sig.beta[3:p,3:p]
		sig.beta <- sig.beta1
		p <- p+1 
	}
	sighc <- matrix(0, nrow=p+q,ncol=p+q)
	th0 <- th[1]; th1 <- th[2]; th2 <- th[3:q]
	sighc[1:q,1:q]<-sig.th; 
	sighc[(q+1):(p+q),(q+1):(p+q)]<-sig.beta
	
	x1 <- 1; x2 <- 0
	
	d1e<-exp(beta2+beta3*x1+th0+th1*x2+th2%*%t(u))
	d2e<-exp(beta2+beta3*x2+th0+th1*x2+th2%*%t(u))
	d1<-d1e/(1+d1e); d2<-d2e/(1+d2e)
	# bd1<-d1-d2; bd2<-x2*(d1-d2); bd3<-u*(d1-d2); bd4<-(x1-x2); bd5<-d1-d2; bd6<-0; bd7<-x1*d1-x2*d2
	bd1<-d1-d2; bd2<-x2*(d1-d2); bd3<-u*c(d1-d2); bd4<-(x1-x2); bd5<-d1-d2; bd6<-x1*d1-x2*d2; bd7 <- c(0,0)
	
	i1e<-exp(beta2+beta3*x1+th0+th1*x1+th2%*%t(u))
	i2e<-exp(beta2+beta3*x1+th0+th1*x2+th2%*%t(u))
	i3e<-exp(th0+th1*x1+th2%*%t(u))
	i4e<-exp(th0+th1*x2+th2%*%t(u))
	i1<-i1e/(1+i1e); i2<-i2e/(1+i2e); i3<-i3e/(1+i3e); i4<-i4e/(1+i4e)
	# bi1<-(i4+i1-i3-i2); bi2<-x2*(i4-i2)+x1*(i1-i3); bi3<-u*(i4+i1)-u*(i3+i2); bi4<-0; bi5<-i1-i2; bi6<-0; bi7<-x1*(i1-i2)
	bi1<-(i4+i1-i3-i2); bi2<-x2*(i4-i2)+x1*(i1-i3); bi3<-u*c(i4+i1)-u*c(i3+i2); bi4<-0; bi5<-i1-i2; bi6<-x1*(i1-i2); bi7<- c(0,0)
	
	NDE_D <- t(c(bd1,bd2,bd3,bd4,bd5,bd6,bd7))
	NIE_D <- t(c(bi1,bi2,bi3,bi4,bi5,bi6,bi7))
	
	NDE_V <- NDE_D%*%sighc%*%t(NDE_D) 
	NIE_V <- NIE_D%*%sighc%*%t(NIE_D)
	list(NDE_V=NDE_V, NIE_V=NIE_V)
}



#########################################################
#           Input variables for each scenario    ´      #
#########################################################

######## Only change this part when running for different scenarios ########
# which sex?
sex_n <- 1
sexs <- c("male","female")

# which outcome?
outcomeName <- "childless_NoART"

# which diagnose to use?
# mod_pattern <- "cond_logit_sibmatch_CorrectSpouse"    

# which sibs to include? Here, we compare four strategies to select siblings from each family.
sib_pattern <- "withchildClosest"  # 1_case : 1_control regarding disease outcome




#########################################################
#           Set working environment                     #
#########################################################

# work directory and function
setwd("REGRESSION/")
r_dir <- "r_files/"


library(data.table)
library(dplyr)
options(tibble.width = Inf)
library(survival)

'%!in%' <- function(x,y)!('%in%'(x,y))
invlogit <- function (x) {1/(1+exp(-x))}




#########################################################
#       Read in data and extract for a specific sex     #
#########################################################

## Demographic info, the same for all analyses, born from 1956 to 1982
indexW <- data.frame(get(load(paste0(r_dir, "indexW_4550_fullsib.Rdata")))) %>% filter(index_sex==sex_n) %>% 
      select(index_id, index_sex, parent_id, index_birth_date, endfollowup4550_date, endfollowup4550_age)
nrow(indexW)  # 443,695 for male


outcome_age_dat <- data.frame(get(load(paste0(r_dir, "indexW_5682_everyone_",outcomeName,".Rdata")))) %>% filter(index_id %in% indexW$index_id) %>%
	                     mutate(outcome_age=outcome_age-1)   # change from age at giving birth to age at pregnancy
nrow(outcome_age_dat)  # 317,772


## Date at diagnosis
ry_first_index_1 <- data.frame(get(load(paste0(r_dir,"indexW_5682_everyone_endpoint.Rdata")))) %>% filter(ID %in% indexW$index_id) %>% select(ID, ENDPOINT, EVENT_F_DATE)
ry_first_index_3 <- data.frame(get(load(paste0(r_dir,"indexW_5682_everyone_Diabetes.Rdata")))) %>% filter(ID %in% indexW$index_id) %>% select(ID, ENDPOINT, EVENT_F_DATE)
ry_first_index_2 <- data.frame(get(load(paste0(r_dir,"indexW_5682_everyone_Qendpoint.20220125.Rdata")))) %>% filter(ID %in% indexW$index_id) %>% select(ID, ENDPOINT, EVENT_F_DATE)
ry_first_index <- rbind(ry_first_index_1, ry_first_index_2, ry_first_index_3)
nrow(ry_first_index)  # 13,704 for diabetes, 2,242 for malformation and 986,015 for other categories 


## Correct for spouseless
marriage <- data.frame(get(load(paste0(r_dir, "indexW_5682_everyone_marriage.Rdata"))))
nrow(marriage)  # 674,384


## Disease endpoint list, significant for both childlessness and spouseless 
disease_lst <- read.table("Diseases_for_mediationAnalysis_SigForChildless.tsv", sep="\t", header=T) %>% data.frame() %>% 
      filter(analysis.Finnish=="Yes") %>% 
      mutate(Endpoint=as.character(Endpoint), LONGNAME=as.character(LONGNAME), sex=as.character(sex), analysis.Finnish=as.character(analysis.Finnish), analysis.Swedish=as.character(analysis.Swedish)) %>% 
      mutate(Endpoint=gsub(" \\(before age 40\\)", "", Endpoint))
if (sex_n==1) {disease_lst <- disease_lst %>% filter(sex=="male")}
if (sex_n==2) {disease_lst <- disease_lst %>% filter(sex=="female")}
dim(disease_lst)   # 52  5




###################################################################################
#        stop/status for outcome events only (ignore disease exposure)            #
###################################################################################

outcome_window <- indexW %>% left_join(outcome_age_dat,by="index_id") %>% 
      mutate(stop=ifelse(is.na(outcome_age), endfollowup4550_age, outcome_age), status=ifelse(is.na(outcome_age), 0, 1)) %>% 
      select(index_id, index_sex, index_birth_date, endfollowup4550_age, parent_id, stop, status) 



for (disease in  disease_lst$Endpoint ) {
	# disease <- "F5_ALCOHOL_DEPENDENCE"
	# disease <- "F5_SCHZPHR"
	
	print(paste0(sexs[sex_n],": ",which(disease==disease_lst$Endpoint),": ", disease))
	
	
	###########################################################################
	#                Add age onset for a specific diseas                      #
	###########################################################################

	## Select disease endpoint of interest and merge with outcome_window
	outcome_window_endpoint <- ry_first_index %>% mutate(index_id=as.numeric(as.character(ID))) %>% filter(ENDPOINT==disease) %>% select(index_id, EVENT_F_DATE) %>% 
	      right_join(outcome_window, by="index_id") %>% 
	      mutate(age_onset=round(as.numeric((as.Date(as.character(EVENT_F_DATE),"%Y%m%d") - as.Date(as.character(index_birth_date),"%Y%m%d"))/365.25),3)) %>% 
	      select(-EVENT_F_DATE)
	
	
	###########################################################################
	#         Sibling-pairs disconcordant on both disease and outcome         #
	###########################################################################

	## Keep only families having full-siblings discordant on outcome status, that is at least one sibling with children and one being childless; 
	parent_outcome1 <- outcome_window_endpoint %>% filter(status==1) %>% select(parent_id) %>% unique()  # 177,327
	parent_outcome0 <- outcome_window_endpoint %>% filter(status==0) %>% select(parent_id) %>% unique()  # 97,080
	p <- outcome_window_endpoint %>% filter(parent_id %in% parent_outcome0$parent_id) %>% filter(parent_id %in% parent_outcome1$parent_id)
	dim(p)  # 188911      9
	
	
	## Within each family, randomly select one sibling with children as control, and the childless siblings with closest birth year with the control as case;
	set.seed(123)
	p_outcome0 <- p %>% filter(status==0) %>% group_by(parent_id) %>% sample_n(1) %>% ungroup()   # 79,122
	p_outcome1 <- p %>% filter(status==1) %>% inner_join(p_outcome0[,c("parent_id","index_birth_date")], by="parent_id") %>% mutate(dif=abs(index_birth_date.y-index_birth_date.x)) %>% 
	                   group_by(parent_id) %>% filter(dif==min(dif)) %>% sample_n(1) %>% 
	                   select(-index_birth_date.y, -dif) %>% rename(index_birth_date="index_birth_date.x") %>% ungroup()
	
	dim(p_outcome0)  # 38928     8
	dim(p_outcome1)  # 38928     8
	
	## disease status at the time point when the event occurs to the case that is the age at first birth for the individual with children
	p_outcome1 <- p_outcome1 %>% mutate(age_onset=ifelse(age_onset>stop, NA, age_onset))
	p_outcome0 <- p_outcome0 %>% inner_join(p_outcome1[,c("parent_id","stop")], by="parent_id") %>% mutate(stop=ifelse(stop.x>stop.y, stop.y, stop.x), age_onset=ifelse(age_onset>stop, NA, age_onset)) %>% select(-stop.x, -stop.y)   # otherwise all diagnose before stop
	
	# combine childless samples and samples with children
	pp <- rbind(p_outcome0, p_outcome1) %>% mutate(disease=ifelse(is.na(age_onset),0,1)) %>% mutate(age=2018-as.numeric(substr(as.character(index_birth_date),1,4)), age2=age^2)
	
	## add has_spouse 
	pp <- pp %>% mutate(has_spouse=ifelse(index_id %in% c(marriage$index_id, marriage$spouse_id), 1, 0)) 
	
	## extract one data which disconcordant on both outcome and disease status for reporting N
	pp_dif <- intersect(pp[pp$disease==0,"parent_id"],pp[pp$disease==1,"parent_id"])
	ppp <- pp %>% filter(parent_id %in% pp_dif$parent_id)   # also disconcordant on disease status, for reporting N (e.g N_Record) which reflect the sample size contribute to the estimate of disease status.
	
	
	
	###########################################################################
	#                         regression analysis                             #
	###########################################################################
	
	pp$parent_id <- factor(pp$parent_id)
	pp$disease <- factor(pp$disease)
	pp <- pp %>% mutate(childless=ifelse(status==0,1,0), spouseless=ifelse(has_spouse==0,1,0))
	
	
	## outcome model ----------
	# not consider interaction
	m_outcome <- clogit(childless ~ disease + spouseless + age + age2  + strata(parent_id), data=pp)  
	beta.N <- summary(m_outcome)$coeff[,1]
	beta.N <- c(beta.N,0)
	sigma.beta <- vcov(m_outcome)
	
	# # consider interactions
	# m_outcome_I <- clogit(childless ~ disease + spouseless + age + age2 + disease*spouseless + strata(parent_id), data=pp) 
	# beta_I.N <- summary(m_outcome_I)$coeff[,1]
	
	
	## mediator model ----------
	m_mediator_glm <- glm(spouseless ~ disease + age + age2, family=binomial, data=pp)  
	th_glm.N <- summary(m_mediator_glm)$coefficients[,1]
	sigma.th <- vcov(m_mediator_glm)
	
	# m_mediator <- clogit(spouseless ~ disease + age + age2 + strata(parent_id), data=pp)  
	# th.N <- summary(m_mediator)$coeff[,1]
	
	
	## baseline covariates for the mediator model ----------
	u <- c(mean(pp$age), mean(pp$age2))
	
	## nature direct and indirect effects
	NCC.CMA <- CDMA(beta.N, th_glm.N, u)
	NCC.NIE <- NCC.CMA$NIE
	NCC.NDE <- NCC.CMA$NDE
	NCC.MP <- NCC.CMA$MP
	
	
	## variance for direct and indirect effects ----------
	dim(u) <- c(1,2)   # need to make u as row vector here
	NCC.VAR <- VAR.D(beta=beta.N[-length(beta.N)], th=th_glm.N, sig.beta=sigma.beta, sig.th=sigma.th, u=u, interaction=F)
	NCC.NIE_VAR <- NCC.VAR$NIE_V
	NCC.NDE_VAR <- NCC.VAR$NDE_V
	
	df <- c(disease, sexs[sex_n], NCC.NIE, NCC.NDE, NCC.NIE_VAR, NCC.NDE_VAR, NCC.MP)
	dim(df) <- c(1,7)
	colnames(df) <- c("Endpoint", "sex", "NIE", "NDE", "NIE.variance", "NDE.variance", "MP")
	print(df)
	
	if (!file.exists(paste0("RESULT_MediationAnalysis.tsv"))){
		write.table(df, paste0("RESULT_MediationAnalysis.tsv"), append=F, quote=F, sep="\t", row.names=F, col.names=T)
	} else {
		write.table(df, paste0("RESULT_MediationAnalysis.tsv"), append=T, quote=F, sep="\t", row.names=F, col.names=F)
	}
	
}

