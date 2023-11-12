
### Association between disease and childlessness with a conditional logistic regression model (sibling design) --------------------------

# To avoid the confounding from reverse causation and the higher prevalence of diseases in older people, 
# We use a matched pair case-control study design within full-siblings and then apply a conditional logistic regression model or a coxPH model stratified by full-sib families. 

# The sibling matching processes could be: 
# 1) Remove individuals with the outcome event before the disease exposure. 
#    By removing individuals with the outcome event prior to the disease exposure, we can avoid assigning potentially affected individuals to be unaffected; 
# 2) Keep only families having full siblings discordant on outcome status, that is, at least one sibling with children and one being childless; 
# 3) Within each family, randomly select one sibling with children as control and the childless siblings with the closest birth year with the control as case; 
# 4) If the case and the control are also discordant on disease status at the time point when the event occurs to the case, that is, the age at first birth for the childless individual,  
#    then we enroll this full-sibling pair into our matched pair case-control study.



#########################################################
#           Input variables for each scenario    ´      #
#########################################################

######## Only change this part when running for different scenarios ########
# which sex?
sex_n <- 1
sexs <- c("male","female")

# which outcome?
# outcomeName <- "childless_NoART"  
outcomeName <- "childless_NoART_NoDeath4550"  


# which diagnose to use? (both patterns are available in for loop)
mod_pattern <- "cond_logit_sibmatch"    

# which sibs to include? Here, we compare four strategies to select siblings from each family.
sib_pattern <- "withchildClosest"  




#########################################################
#           Set working environment                     #
#########################################################

# for disease endpoint with too few cases
df_NA <- matrix(NA, ncol=33, nrow=1)
colnames(df_NA) <- c("Estimate","HR","SE","Z","P_val","HR_025","HR_975","avg_disease","avg_age","avg_age2","beta_age","beta_age2","vcv_11","vcv_12","vcv_13","vcv_22","vcv_23","vcv_33","prob_disease0","prob_disease1","prob_disease","probSE_disease","N_Family","N_Sample","N_Record","N_00","N_11","N_01","N_10","model","sex","LRS","Endpoint")


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
if (outcomeName=="childless_NoART_NoDeath4550"){
	indexW <- data.frame(get(load(paste0(r_dir, "indexW_4550_fullsib.Rdata")))) %>% filter(index_sex==sex_n) %>% 
	      filter( (is.na(index_death_date)) | (!is.na(index_death_date) & index_death_date!=endfollowup4550_date) ) %>% 
	      select(index_id, index_sex, parent_id, index_birth_date, endfollowup4550_date, endfollowup4550_age)
	nrow(indexW)  # 208,497 for male
} else {
	indexW <- data.frame(get(load(paste0(r_dir, "indexW_4550_fullsib.Rdata")))) %>% filter(index_sex==sex_n) %>% 
	      select(index_id, index_sex, parent_id, index_birth_date, endfollowup4550_date, endfollowup4550_age)
	nrow(indexW)  # 223,196 for male
}

outcome_age_dat <- data.frame(get(load(paste0(r_dir, "indexW_5682_everyone_",substr(outcomeName,1,15),".Rdata")))) %>% filter(index_id %in% indexW$index_id) %>%
	                     mutate(outcome_age=outcome_age-1)   # change from age at giving birth to age at pregnancy
nrow(outcome_age_dat)  # 162,346


## Date at diagnosis
ry_first_index_1 <- data.frame(get(load(paste0(r_dir,"indexW_5682_everyone_endpoint.Rdata")))) %>% filter(ID %in% indexW$index_id) %>% select(ID, ENDPOINT, EVENT_F_DATE)
ry_first_index_2 <- data.frame(get(load(paste0(r_dir,"indexW_5682_everyone_Qendpoint.Rdata")))) %>% filter(ID %in% indexW$index_id) %>% select(ID, ENDPOINT, EVENT_F_DATE)
ry_first_index_3 <- data.frame(get(load(paste0(r_dir,"indexW_5682_everyone_Diabetes_op1_age40.Rdata")))) %>% filter(ID %in% indexW$index_id) %>% select(ID, ENDPOINT, EVENT_F_DATE)
ry_first_index <- rbind(ry_first_index_1, ry_first_index_2, ry_first_index_3)
nrow(ry_first_index)  # 13,704 for diabetes, 2,242 for malformation and 986,015 for other categories 


## Disease endpoint list (skip endpoints already done)
disease30_lst <- ry_first_index %>% mutate(ENDPOINT=as.character(ENDPOINT)) %>% group_by(ENDPOINT) %>% count() %>% filter(n>=30)  # 1,253
mod_pattern <- "cond_logit_sibmatch"
if (file.exists(paste0("RESULT_",mod_pattern,"_", outcomeName ,"_", sexs[sex_n], ".tsv"))){
	done <- read.table(paste0("RESULT_",mod_pattern,"_", outcomeName ,"_", sexs[sex_n], ".tsv"), sep="\t", header=T)
	disease30_lst <- disease30_lst %>% filter(ENDPOINT %!in% done$Endpoint)
}  # 1,186 




###################################################################################
#        stop/status for outcome events only (ignore disease exposure)            #
###################################################################################

outcome_window <- indexW %>% left_join(outcome_age_dat,by="index_id") %>% 
      mutate(stop=ifelse(is.na(outcome_age), endfollowup4550_age, outcome_age), status=ifelse(is.na(outcome_age), 0, 1)) %>% 
      select(index_id, index_sex, index_birth_date, endfollowup4550_age, parent_id, stop, status) 



for (disease in disease30_lst$ENDPOINT ) {
	# disease <- "F5_SCHIZPER"
	# disease <- "N14_MALEINFERT"
	print(paste0(sexs[sex_n],": ",which(disease==disease30_lst$ENDPOINT),": ", disease))



	###########################################################################
	#                Add age onset for a specific disease                      #
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
	dim(p)  # 93270      9


	## Within each family, randomly select one sibling with children as control, and the childless siblings with closest birth year with the control as case;
	set.seed(123)
	p_outcome0 <- p %>% filter(status==0) %>% group_by(parent_id) %>% sample_n(1) %>% ungroup()   # 79,122
	p_outcome1 <- p %>% filter(status==1) %>% inner_join(p_outcome0[,c("parent_id","index_birth_date")], by="parent_id") %>% mutate(dif=abs(index_birth_date.y-index_birth_date.x)) %>% 
	                   group_by(parent_id) %>% filter(dif==min(dif)) %>% sample_n(1) %>% 
	                   select(-index_birth_date.y, -dif) %>% rename(index_birth_date="index_birth_date.x") %>% ungroup()
	
	dim(p_outcome0)  # 38762     10
	dim(p_outcome1)  # 38762     10

	for (mod_pattern in c("cond_logit_4550","cond_logit_sibmatch")){
	
		if (mod_pattern=="cond_logit_sibmatch"){
			## disease status at the time point when the event occurs to the case that is the age at first birth for the individual with children
			p_outcome1 <- p_outcome1 %>% mutate(age_onset=ifelse(age_onset>stop, NA, age_onset))
			p_outcome0 <- p_outcome0 %>% inner_join(p_outcome1[,c("parent_id","stop")], by="parent_id") %>% mutate(stop=ifelse(stop.x>stop.y, stop.y, stop.x), age_onset=ifelse(age_onset>stop, NA, age_onset)) %>% select(-stop.x, -stop.y)
		}   # otherwise all diagnose before stop
	
		
		pp <- rbind(p_outcome0, p_outcome1) %>% mutate(disease=ifelse(is.na(age_onset),0,1)) %>% mutate(age=2018-as.numeric(substr(as.character(index_birth_date),1,4)), age2=age^2)
		pp_dif <- intersect(pp[pp$disease==0,"parent_id"],pp[pp$disease==1,"parent_id"])
		ppp <- pp %>% filter(parent_id %in% pp_dif$parent_id)   # also disconcordant on disease status, for reporting N (e.g N_Record) which reflect the sample size contribute to the estimate of disease status.
	
	
	
	
		###########################################################################
		#                         regression analysis                             #
		###########################################################################
	
		if (ppp %>% select(index_id) %>% unique() %>% nrow()>=30){  # only analyze disease with more than 30 cases
			pp$parent_id <- factor(pp$parent_id)
			pp$disease <- factor(pp$disease)
			m_mod_sib <- clogit(status ~ disease + age + age2 + strata(parent_id), data=pp)  # use the dataset disconcordant on outcome
			df <- as.data.frame(t(as.data.frame(summary(m_mod_sib)$coeff["disease1",])))
			rownames(df) <- "Endpoint"
			df <- df[,c("coef", "exp(coef)", "se(coef)", "z","Pr(>|z|)")]
			colnames(df) <- c("Estimate", "HR", "SE","Z", "P_val")
			df$HR <- summary(m_mod_sib)$conf.int["disease1","exp(coef)"]
			df$HR_025 <- summary(m_mod_sib)$conf.int["disease1","lower .95"]
			df$HR_975 <- summary(m_mod_sib)$conf.int["disease1","upper .95"]
			
			df$avg_disease <- mean(as.numeric(as.character(pp$disease)))
			df$avg_age <- mean(pp$age)
			df$avg_age2 <- mean(pp$age2)
			
			df$beta_age <- summary(m_mod_sib)$coeff["age","coef"]
			df$beta_age2 <- summary(m_mod_sib)$coeff["age2","coef"]
	
			vcv <- vcov(m_mod_sib)
			df$vcv_11 <- vcv[1,1]
			df$vcv_12 <- vcv[1,2]
			df$vcv_13 <- vcv[1,3]
			df$vcv_22 <- vcv[2,2]
			df$vcv_23 <- vcv[2,3]
			df$vcv_33 <- vcv[3,3]
			
			# probability and marginal effects		
			allmean <- data.frame(disease=c(0,1), age=rep(mean(pp$age),2), age2=rep(mean(pp$age2),2))
			allmean$disease <- as.factor(allmean$disease)
			lp <- predict(m_mod_sib, newdata=allmean, type="lp", reference="sample") 
			df$prob_disease0 <- invlogit(lp[1])   # -0.1284011
			df$prob_disease1 <- invlogit(lp[2])   # -0.1284011
			df$prob_disease <- invlogit(lp[2]) - invlogit(lp[1])   # -0.1284011
	
			# how to calculate SE (following the same way as logistic regression)
			xm <- t(data.frame(disease=df$avg_disease, age=df$avg_age, age2=df$avg_age2))
			be <- as.matrix(na.omit(coef(m_mod_sib))) 
			disx0 <- disx1 <- xm
			disx0["disease",] <- 0
			disx1["disease",] <- 1
			gr <- dlogis(t(be) %*% disx1) %*% t(disx1) - dlogis(t(be) %*% disx0) %*% t(disx0)
			df$probSE_disease <- sqrt(gr %*% vcv %*% t(gr))   # 0.0713484
	
		} else {
			df <- df_NA
		}
		
		df[1,"N_Family"] <- nrow(unique(ppp[,c("parent_id")]))  # 2*2 table for large SE results
		df[1,"N_Sample"] <- nrow(unique(ppp[,c("index_id")]))   
		df[1,"N_Record"] <- nrow(ppp)  # change to N cases
		df[1,"N_00"] <- ppp %>% filter(disease==0 & status==0) %>% nrow()
		df[1,"N_11"] <- ppp %>% filter(disease==1 & status==1) %>% nrow()
		df[1,"N_01"] <- ppp %>% filter(disease==0 & status==1) %>% nrow()
		df[1,"N_10"] <- ppp %>% filter(disease==1 & status==0) %>% nrow()
		df[1,"model"] <- mod_pattern
		df[1,"sex"] <- sexs[sex_n] 
		df[1,"LRS"] <- outcomeName
		df[1,"Endpoint"] <- disease
	
		print(df)
		if (!file.exists(paste0("RESULT_",mod_pattern,"_", outcomeName ,"_", sexs[sex_n], ".tsv"))){
			write.table(df, paste0("RESULT_",mod_pattern,"_", outcomeName ,"_", sexs[sex_n], ".tsv"), append=F, quote=F, sep="\t", row.names=F, col.names=T)
		} else {
			write.table(df, paste0("RESULT_",mod_pattern,"_", outcomeName ,"_", sexs[sex_n], ".tsv"), append=T, quote=F, sep="\t", row.names=F, col.names=F)
		}
	
	}

}



