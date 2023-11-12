
### Association between diseases and childlessness using time-varying coxPH model --------------------------

## Model 2: cond_coxvary_sibmatch: (need to use the entire population (born in 1956-1982) with same-sex full sibs since it can manage the censoring) -------------------
# cox HR model stratified by full-sibling families for childlessness in sibling pairs disconcordant on outcomes, using a sibling-match design and considering disease status as time-varing covariable
# To avoid the confounding from reverse causation and the higher prevalence of diseases in older people, 
#    we use a matched pair case-control study design within full-siblings and then apply cox HR model stratified by full-sib families. 
# The sibling matching processes could be: 
# 1) remove individuals with the outcome event prior to the disease exposure. 
#    By removing individuals with the outcome event prior to the disease exposure, we can avoid assigning potential affected individuals to be unaffected; 
# 2) Keep only families having full-siblings discordant on outcome status, that is at least one sibling with children and one being childless; 
# 3) Within each family, randomly select one sibling with children as control, and the childless siblings with closest birth year with the control as case; 
# 4) If the case and the control are also discordant on disease status at the time point when the event occurs to the case that is the age at first birth for case for childless, 
#    then we enroll this full-sibling pair into our matched pair case-control study;




#########################################################
#           Input variables for each scenario    ´      #
#########################################################

######## Only change this part when running for different scenarios ########
# which sex?
sex_n <- 1
sexs <- c("male","female")

# which outcome?
# outcomeName <- "childless"  
outcomeName <- "childless_NoART"  

# which diagnose to use?
# mod_pattern <- "cond_coxvary_sibmatch"  

# which sibs to include? Here, we compare four strategies to select siblings from each family.
sib_pattern <- "withchildClosest"  # 1_case : 1_control regarding disease outcome


#############################################################################




#########################################################
#           Set working environment                     #
#########################################################

# for disease endpoint with too few cases
df_NA <- matrix(NA, ncol=29, nrow=1)
colnames(df_NA) <- c("Estimate","HR","SE","Z","P_val","HR_025","HR_975", "avg_disease","avg_age","avg_age2","beta_age","beta_age2", "vcv_11","vcv_12","vcv_13","vcv_22","vcv_23","vcv_33", "N_Family","N_Sample","N_Record","N_00","N_11","N_01","N_10","model","sex","LRS","Endpoint")


# work directory and function
setwd("/homes/aliu/DSGE_LRS/output/registry_edit/REGRESSION/")
r_dir <- "/homes/aliu/DSGE_LRS/input/r_files/"


library(data.table)
library(dplyr)
options(tibble.width = Inf)
library(survival)

'%!in%' <- function(x,y)!('%in%'(x,y))



#########################################################
#       Read in data and extract for a specific sex     #
#########################################################

## Demographic info, the same for all analyses, born from 1956 to 1982
indexW <- data.frame(get(load(paste0(r_dir, "indexW_5682_fullsib.Rdata")))) %>% filter(index_sex==sex_n) %>% 
      select(index_id, index_sex, parent_id, index_birth_date, endfollowup4550_date, endfollowup4550_age)
nrow(indexW)  # 444,031 for male


outcome_age_dat <- data.frame(get(load(paste0(r_dir, "indexW_5682_everyone_",outcomeName,".Rdata")))) %>% filter(index_id %in% indexW$index_id) %>%
	                     mutate(outcome_age=outcome_age-1)   # change from age at giving birth to age at pregnancy
nrow(outcome_age_dat)  # 317,984


## Date at diagnosis
ry_first_index_1 <- data.frame(get(load(paste0(r_dir,"indexW_5682_everyone_endpoint.Rdata")))) %>% filter(ID %in% indexW$index_id) %>% select(ID, ENDPOINT, EVENT_F_DATE)
ry_first_index_2 <- data.frame(get(load(paste0(r_dir,"indexW_5682_everyone_Qendpoint.Rdata")))) %>% filter(ID %in% indexW$index_id) %>% select(ID, ENDPOINT, EVENT_F_DATE)
ry_first_index_3 <- data.frame(get(load(paste0(r_dir,"indexW_5682_everyone_Diabetes.Rdata")))) %>% filter(ID %in% indexW$index_id) %>% select(ID, ENDPOINT, EVENT_F_DATE)
ry_first_index <- rbind(ry_first_index_1, ry_first_index_2, ry_first_index_3)
nrow(ry_first_index)  # 13,704 for diabetes, 2,242 for malformation and 986,015 for other categories 


## Disease endpoint list (skip endpoints already done)
mod_pattern <- "cond_coxvary_sibmatch"   
disease30_lst <- ry_first_index %>% mutate(ENDPOINT=as.character(ENDPOINT)) %>% group_by(ENDPOINT) %>% count() %>% filter(n>=30)  # 1,253
if (file.exists(paste0("RESULT_",mod_pattern,"_", outcomeName ,"_", sexs[sex_n], ".tsv"))){
	done <- read.table(paste0("RESULT_",mod_pattern,"_", outcomeName ,"_", sexs[sex_n], ".tsv"), header=T, sep="\t")
	disease30_lst <- disease30_lst %>% filter(ENDPOINT %!in% done$Endpoint)
}  # 1,398



###################################################################################
#        start/stop/status for outcome events only (ignore disease exposure)      #
###################################################################################

outcome_window <- indexW %>% left_join(outcome_age_dat,by="index_id") %>% 
      mutate(start=16, 
             stop=ifelse(is.na(outcome_age), endfollowup4550_age, outcome_age), 
             status=ifelse(is.na(outcome_age), 0, 1)) %>% 
      select(index_id, index_sex, index_birth_date, endfollowup4550_age, parent_id, start, stop, status) 


for (disease in disease30_lst$ENDPOINT ) {
	# disease <- "F5_SCHIZPER"
	# disease <- "N14_MALEINFERT"
	print(paste0(sexs[sex_n],": ",which(disease==disease30_lst$ENDPOINT),": ", disease))


	###########################################################################
	#             Add time-varing disease exposure variable                   #
	###########################################################################

	## Select disease endpoint of interest and merge with outcome_window
	outcome_window_endpoint <- ry_first_index %>% mutate(index_id=as.numeric(as.character(ID))) %>% filter(ENDPOINT==disease) %>% select(index_id, EVENT_F_DATE) %>% 
	      right_join(outcome_window, by="index_id") %>% 
	      mutate(age_onset=round(as.numeric((as.Date(as.character(EVENT_F_DATE),"%Y%m%d") - as.Date(as.character(index_birth_date),"%Y%m%d"))/365.25),3)) %>% 
	      select(-EVENT_F_DATE)

	
	## remove diagnose after stop 
	outcome_window_endpoint <- outcome_window_endpoint %>% mutate(age_onset=ifelse(age_onset>stop,NA,age_onset)) %>% filter(stop>start)
	dim(outcome_window_endpoint)   # 443,824      9



	###########################################################################
	#         Sibling-pairs disconcordant on outcome                          #
	###########################################################################

	## Keep only families having full-siblings discordant on outcome status, that is at least one sibling with children and one being childless; 
	parent_outcome1 <- outcome_window_endpoint %>% filter(status==1) %>% select(parent_id) %>% unique()  # 177,327
	parent_outcome0 <- outcome_window_endpoint %>% filter(status==0) %>% select(parent_id) %>% unique()  # 97,080
	p <- outcome_window_endpoint %>% filter(parent_id %in% parent_outcome0$parent_id) %>% filter(parent_id %in% parent_outcome1$parent_id)
	dim(p)  # 189,077      9


	## Within each family, randomly select one sibling with children as control, and the childless siblings with closest birth year with the control as case;
	p_outcome0 <- p %>% filter(status==0) %>% group_by(parent_id) %>% sample_n(1) %>% ungroup()   # 79,122
	p_outcome1 <- p %>% filter(status==1) %>% inner_join(p_outcome0[,c("parent_id","index_birth_date")], by="parent_id") %>% mutate(dif=abs(index_birth_date.y-index_birth_date.x)) %>% 
	                   group_by(parent_id) %>% filter(dif==min(dif)) %>% sample_n(1) %>% 
	                   select(-index_birth_date.y, -dif) %>% rename(index_birth_date="index_birth_date.x") %>% ungroup()

	dim(p_outcome0)  # 79,207     10
	dim(p_outcome1)  # 79,207     10

	## disease status at the time point when the event occurs to the case that is the age at first birth for the individual with children
	p_outcome0 <- p_outcome0 %>% inner_join(p_outcome1[,c("parent_id","stop")], by="parent_id") %>% mutate(stop=ifelse(stop.x>stop.y, stop.y, stop.x), age_onset=ifelse(age_onset>stop, NA, age_onset)) %>% select(-stop.x, -stop.y)

	## combine childless samples and samples with children
	ppp <- rbind(p_outcome0, p_outcome1) %>% mutate(disease=ifelse(is.na(age_onset),0,1)) %>% 
	            mutate(age=2018-as.numeric(substr(as.character(index_birth_date),1,4)), age2=age^2)
	dim(ppp)  # 158,414     12

	ppp_dif <- intersect(ppp[ppp$disease==0,"parent_id"],ppp[ppp$disease==1,"parent_id"])   # 234   1
	pppp <- ppp %>% filter(parent_id %in% ppp_dif$parent_id)   # 468  12, also disconcordant on disease status, for reporting N (e.g N_Record) which reflect the sample size contribute to the estimate of disease status.


	###########################################################################
	#             Add time-varing disease exposure variable                   #
	###########################################################################

	## set age_onset to 16 for disease diganose before age 16
	ppp <- ppp %>% mutate(age_onset=ifelse(age_onset<16, 16, age_onset))

	# (1) No diagnose (|-----|)
	ppp_1 <- ppp %>% filter(is.na(age_onset)) %>% mutate(disease=0)

	# (2) Windows with start<age_onset (|---*..| or |-----|*)
	ppp_2 <- ppp %>% filter(!is.na(age_onset) & start<age_onset) %>% mutate(disease=0) %>% 
	      mutate(stop=ifelse(stop<age_onset, stop, age_onset), status=ifelse(stop<age_onset, status, 0))
	
	# (3) Windows with stop>age_onset (|..*---| or *|-----|)
	ppp_3 <- ppp %>% filter(!is.na(age_onset) & stop>age_onset) %>% mutate(disease=1) %>% 
	      mutate(start=ifelse(start<age_onset, age_onset, start))

	ppp <- bind_rows(ppp_1, ppp_2, ppp_3) 
	dim(ppp)  # 158,646       8
	rm(ppp_1, ppp_2, ppp_3)
	print("This step 2 is done.")



	###########################################################################
	#  Survival analysis considering disease status as time-varing covariate  #
	###########################################################################

	if (ppp %>% filter(disease==1) %>% select(index_id) %>% unique() %>% nrow()>=30){
		ppp$parent_id <- factor(ppp$parent_id)
		m_mod_sib <- coxph(Surv(start,stop,status) ~ age + age2 + cluster(index_id) + disease + strata(parent_id), data=ppp)
		df <- as.data.frame(t(as.data.frame(summary(m_mod_sib)$coeff["disease",])))
		rownames(df) <- "Endpoint"
		df <- df[,c("coef", "exp(coef)", "robust se", "z","Pr(>|z|)")]
		colnames(df) <- c("Estimate", "HR", "SE","Z", "P_val")
		df$HR <- summary(m_mod_sib)$conf.int["disease","exp(coef)"]
		df$HR_025 <- summary(m_mod_sib)$conf.int["disease","lower .95"]
		df$HR_975 <- summary(m_mod_sib)$conf.int["disease","upper .95"]
		
		df$avg_disease <- mean(as.numeric(as.character(ppp$disease)))
		df$avg_age <- mean(ppp$age)
		df$avg_age2 <- mean(ppp$age2)
		
		df$beta_age <- summary(m_mod_sib)$coeff["age","coef"]
		df$beta_age2 <- summary(m_mod_sib)$coeff["age2","coef"]
		
		vcv <- vcov(m_mod_sib)
		df$vcv_11 <- vcv[1,1]
		df$vcv_12 <- vcv[1,2]
		df$vcv_13 <- vcv[1,3]
		df$vcv_22 <- vcv[2,2]
		df$vcv_23 <- vcv[2,3]
		df$vcv_33 <- vcv[3,3]
		
	} else {
		df <- df_NA
	}
	
	df[1,"N_Family"] <- pppp %>% select(parent_id) %>% unique() %>% nrow()  # 2*2 table for large SE results
	df[1,"N_Sample"] <- pppp %>% select(index_id) %>% unique() %>% nrow()  
	df[1,"N_Record"] <- nrow(pppp)
	df[1,"N_00"] <- pppp %>% filter(disease==0 & status==0) %>% nrow()
	df[1,"N_11"] <- pppp %>% filter(disease==1 & status==1) %>% nrow()
	df[1,"N_01"] <- pppp %>% filter(disease==0 & status==1) %>% nrow()
	df[1,"N_10"] <- pppp %>% filter(disease==1 & status==0) %>% nrow()
	df[1,"model"] <- paste0("cox_vary_",mod_pattern)
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


