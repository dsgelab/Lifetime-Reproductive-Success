# LRS results visualization 

#### packages ####
library(shiny)
library(shinyalert)
library(ggplot2)
library(data.table)
library(dplyr)
library(shinyWidgets)
library(purrr)
library(gridExtra)
library(cowplot)
library(RColorBrewer)



###############################################################################################
####                         FUNCTIONS TO DEFINE PLOTS                                     ####
###############################################################################################

## Main plot
# set color key to be the same as the manuscript
nb.cols <- 16
mycolors <- data.frame(matrix(c(colorRampPalette(brewer.pal(8, "Set1"))(nb.cols), c("Blood and immune mechanism", "Circulatory system", "Digestive system","Ear and mastoid process","Endocrine-nutritional-metabolic","Eye and adnexa",
      "Genitourinary system", "Infectious-parasitic", "Mental-behavioural","Musculoskeletal system", "Neoplasms","Nervous system", "Other","Respiratory system", "Skin and subcutaneous tissue", "Congenital anomalies")), byrow=F, ncol=2))
colnames(mycolors) <- c("mycol","grp_text")

mycolors <- mycolors %>% mutate(grp_text=factor(grp_text, levels=c("Infectious-parasitic","Neoplasms", "Blood and immune mechanism","Endocrine-nutritional-metabolic", "Mental-behavioural", 
                                "Nervous system", "Eye and adnexa", "Ear and mastoid process", "Circulatory system", "Respiratory system", "Digestive system", 
                                "Skin and subcutaneous tissue", "Musculoskeletal system", "Genitourinary system", "Congenital anomalies", "Other"))) %>% 
                         arrange(-desc(grp_text))


plot_cap <- expression(paste('*Only disease diagnoses that are significantly associated with childlessness after multiple-testing correction (P<1.5x',10^{-4},') are colored.'))

mhd_index <- function(dat, names_dis, plot_title, plot_ylab){
	### names_dis is either names_dis_m or names_dis_f
	# labels <- names_dis %>% filter(Endpoint %in% dat$Endpoint) %>% arrange(match(Endpoint,dat$Endpoint)) %>% mutate(seq=1:n()) %>% group_by(grp_text) 
	# labels_axis <- labels %>% summarize(xlab_index=median(seq))
	# labels_lines <- labels %>% summarize(xlab_index=max(seq))
	
	hlines <- aggregate(seq ~ grp_text, dat, max)[,"seq"] + 0.5
	xlabs <- 0.5 * aggregate(seq ~ grp_text, dat, min)[,"seq"] + 0.5 * aggregate(seq ~ grp_text, dat, max)[,"seq"]
	xlabs_text <- aggregate(seq ~ grp_text, dat, min)[,"grp_text"]
	
	pp <- ggplot(data=dat, aes(x=seq, y=OR.meta, group=grp_text)) +
	      geom_point(aes(color=grp_text, fill=grp_text, alpha=al), shape=25, size=2) + 
	      labs(y=plot_ylab, x="", title=plot_title, subtitle="(Click on the plot to select a disease)", caption=plot_cap) + 
	      geom_hline(yintercept=1, linetype="dotted", color="grey", size=0.6) +
	      #  geom_vline(xintercept=labels_lines$xlab_index, linetype="dotted", color="grey", size=0.6) + 
	      #  scale_x_continuous(breaks=labels_axis$xlab_index, label=labels_axis$grp_text) + 
	      geom_vline(xintercept=hlines, linetype="dotted", color="grey", size=0.6) + 
	      scale_x_continuous(breaks=xlabs, label=xlabs_text) +
	      scale_y_log10(limits=c(0.4,25),breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25),labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25))+ 
	      scale_color_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
	      scale_fill_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
	      theme_classic()
	      pp <- pp + theme(axis.text.x=element_text(hjust=0, vjust=0, size=8, angle=-60, face="bold"), axis.title.y=element_text(size=8, face="bold"), plot.title=element_text(size=12, face="bold"), plot.caption.position="plot", plot.caption=element_text(hjust=0, vjust=6, size=9, face="italic"), plot.margin=unit(c(0,0,0.1,0.01), "null"))
	      # pp <- pp + theme(axis.text.x=element_text(hjust=0, vjust=0, size=8, angle=-60, face="bold"), axis.title.y=element_text(size=8, face="bold"), plot.title=element_text(size=12, face="bold"))
	      pp <- pp + theme(legend.position="none")
	return(pp)
}


# To compare measures
# get ylim values to have same y-axis scale for 
# 3 plots next to each other for both sexes, if both are plotted
get_ylims <- function(sex, index_meta_m_endpoint, index_meta_m_mod_endpoint, index_meta_m_cond_endpoint, med_meta_m_endpoint, index_meta_f_endpoint, index_meta_f_mod_endpoint, index_meta_f_cond_endpoint, med_meta_f_endpoint){ 
  if (sex == "female"){
    # use first only female data
    min025 <- min(
      index_meta_f_endpoint$OR025.meta, 
      index_meta_f_cond_endpoint$OR025.Swedish, 
      index_meta_f_cond_endpoint$OR025.Finnish,
      med_meta_f_endpoint$OR025.meta,
      na.rm = TRUE)
    
    max975 <- max(
      index_meta_f_endpoint$OR975.meta, 
      index_meta_f_cond_endpoint$OR975.Swedish, 
      index_meta_f_cond_endpoint$OR975.Finnish,
      med_meta_f_endpoint$OR975.meta,
      na.rm = TRUE)
    
    if(nrow(index_meta_m_cond_endpoint) != 0){
      # don't test for significance because plots are shown 
      # even if results of main analysis are not significant
      # if(index_meta_m_cond_endpoint$al == 0.8){
      # get values using both sexes
      min025_m <- min(
        index_meta_m_endpoint$OR025.meta,
        index_meta_m_cond_endpoint$OR025.Swedish,
        index_meta_m_cond_endpoint$OR025.Finnish,
        med_meta_m_endpoint$OR025.meta,
        na.rm = TRUE)
      
      max975_m <- max(
        index_meta_m_endpoint$OR975.meta,
        index_meta_m_cond_endpoint$OR975.Swedish,
        index_meta_m_cond_endpoint$OR975.Finnish,
        med_meta_m_endpoint$OR975.meta,
        na.rm = TRUE)
      
      min025 <- min(min025, min025_m)
      max975 <- max(max975, max975_m)
      #}
    }
  }
  
  if (sex == "male"){
    min025 <-   min025 <- min(
      index_meta_m_endpoint$OR025.meta,
      index_meta_m_cond_endpoint$OR025.Swedish,
      index_meta_m_cond_endpoint$OR025.Finnish,
      med_meta_m_endpoint$OR025.meta,
      na.rm = TRUE)
    
    max975 <- max(
      index_meta_m_endpoint$OR975.meta,
      index_meta_m_cond_endpoint$OR975.Swedish,
      index_meta_m_cond_endpoint$OR975.Finnish,
      med_meta_m_endpoint$OR975.meta,
      na.rm = TRUE)
    
    if(nrow(index_meta_f_cond_endpoint) != 0){
      # don't test for significance because plots are shown 
      # even if results of main analysis are not significant
      # if(index_meta_f_cond_endpoint$al == 0.8){
      # get values using both sexes
      min025_f <- min(
        index_meta_f_endpoint$OR025.meta, 
        index_meta_f_cond_endpoint$OR025.Swedish, 
        index_meta_f_cond_endpoint$OR025.Finnish,
        med_meta_f_endpoint$OR025.meta,
        na.rm = TRUE)
      
      max975_f <- max(
        index_meta_f_endpoint$OR975.meta, 
        index_meta_f_cond_endpoint$OR975.Swedish, 
        index_meta_f_cond_endpoint$OR975.Finnish,
        med_meta_f_endpoint$OR975.meta,
        na.rm = TRUE)
      
      min025 <- min(min025, min025_f)
      max975 <- max(max975, max975_f)
      #}
    }
  }
  
  # ylim_min
  ylim_min <- round(min025 - 0.05, digits = 1)
  
  # ylim_max
  ylim_max <- round(max975 + 0.05, digits = 1)
  
  return(c(ylim_min, ylim_max))
}



comp_measures <- function(dat, ylims, title_sex){
	# save data used for plot. needed for hovering
	data_comp_meas_f <<- dat

	pp <- ggplot(data=dat, aes(x=LRS, y=OR.meta, ymin=OR025.meta, ymax=OR975.meta)) +
		# geom_pointrange() +
		geom_pointrange(shape=25, size=0.35, aes(fill=grp_text, color=grp_text)) +
		geom_hline(yintercept=1, linetype="dotted", color="grey", size=0.3) +
		scale_color_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
		scale_fill_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
		theme_classic() + 
		scale_x_discrete("",labels=c("childless"="Childless\n", "spouseless"="Partnerless", "childless_WithSpouse"="Childless\n(with partner)")) + 
		theme(axis.text.x=element_text(size=8, face="bold"), axis.title.y=element_text(size=8, face="bold"), plot.title=element_text(size=9, face="bold")) + 
		labs(y="Odds Ratio", x="", title="Relationship of disease diagnoses\nwith different reproductive outcomes") + 
		# ylim(ylims) + 
		scale_y_log10(limits=c(0.4,25),breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25),labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25))+ 
		theme(legend.position="none")
	
	return(pp)
}



# To compare models for has_child 
comp_models <- function(dat, ylims, title_sex){
	# save data used for plot. needed for hovering
	data_comp_meas_f <<- dat
	
	max_OR <- max(dat[,"OR975.meta"])
	if (max_OR<25) {
		pp <- ggplot(data=dat, aes(x=model, y=OR.meta, ymin=OR025.meta, ymax=OR975.meta)) +
			# geom_pointrange() +
			geom_pointrange(shape=25, size=0.35, aes(fill=grp_text, color=grp_text)) + 
			geom_hline(yintercept = 1, linetype="dotted", color = "grey", size=0.3) +
			scale_color_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
			scale_fill_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
			theme_classic() +
			scale_x_discrete("",labels=c("Main analysis"="Main\nanalysis", "Alive by age 45"="Alive\nby age 45", "Alive by age 50"="Alive\nby age 50", 
			                             "Time-dependent Cox PH"="Stratified\nCox PH", "Population-based"="Population\nbased")) + 
			theme(axis.text.x=element_text(size=8, face="bold"), axis.title=element_text(size=8, face="bold"), plot.title=element_text(size=9, face="bold")) + 
			labs(y="Odds Ratio of being childless", x="", title="Relationship of disease diagnoses\nwith childlessness using different models") + 
			# ylim(ylims)
			scale_y_log10(limits=c(0.4,25),breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25),labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25)) + 
			theme(legend.position="none")
	} 
	
	if (max_OR>=25 & max_OR<50) {
		pp <- ggplot(data=dat, aes(x=model, y=OR.meta, ymin=OR025.meta, ymax=OR975.meta)) +
			# geom_pointrange() +
			geom_pointrange(shape=25, size=0.35, aes(fill=grp_text, color=grp_text)) + 
			geom_hline(yintercept = 1, linetype="dotted", color = "grey", size=0.3) +
			scale_color_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
			scale_fill_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
			theme_classic() +
			scale_x_discrete("",labels=c("Main analysis"="Main\nanalysis", "Alive by age 45"="Alive\nby age 45", "Alive by age 50"="Alive\nby age 50", 
			                             "Time-dependent Cox PH"="Stratified\nCox PH", "Population-based"="Population\nbased")) + 
			theme(axis.text.x=element_text(size=8, face="bold"), axis.title=element_text(size=8, face="bold"), plot.title=element_text(size=9, face="bold")) + 
			labs(y="Odds Ratio of being childless", x="", title="Relationship of disease diagnoses\nwith childlessness using different models") + 
			# ylim(ylims)
			scale_y_log10(limits=c(0.4,50),breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25,35,50),labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25,35,50)) + 
			theme(legend.position="none")
	} 
	
	if (max_OR>=50) {
	pp <- ggplot(data=dat, aes(x=model, y=OR.meta, ymin=OR025.meta, ymax=OR975.meta)) +
		# geom_pointrange() +
		geom_pointrange(shape=25, size=0.35, aes(fill=grp_text, color=grp_text)) + 
		geom_hline(yintercept = 1, linetype="dotted", color = "grey", size=0.3) +
		scale_color_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
		scale_fill_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
		theme_classic() +
		scale_x_discrete("",labels=c("Main analysis"="Main\nanalysis", "Alive by age 45"="Alive\nby age 45", "Alive by age 50"="Alive\nby age 50", 
		                             "Time-dependent Cox PH"="Stratified\nCox PH", "Population-based"="Population\nbased")) + 
		theme(axis.text.x=element_text(size=8, face="bold"), axis.title=element_text(size=8, face="bold"), plot.title=element_text(size=9, face="bold")) + 
		labs(y="Odds Ratio of being childless", x="", title="Relationship of disease diagnoses\nwith childlessness using different models") + 
		# ylim(ylims)
		scale_y_log10(limits=c(0.4,75),breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25,35,50,75),labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25,35,50,75)) + 
		theme(legend.position="none")
	}
	return(pp)
}


# To compare Swedish and Finnish
comp_fin_swe <- function(dat, ylims){
	# return plot only data from both countries are available
	if (is.na(dat$OR.Swedish)){
		return(as.character("Sweden"))
	} else if (is.na(dat$OR.Finnish)){
		return("Finland")
	} else{
	
	tt <- data.frame(OR=c(dat$OR.Swedish,dat$OR.Finnish),
	                 OR025=c(dat$OR025.Swedish,dat$OR025.Finnish),
	                 OR975=c(dat$OR975.Swedish,dat$OR975.Finnish),
	                 P_val=c(dat$P_val.Swedish, dat$P_val.Finnish),
	                 group=c("Sweden","          Finland"),
	                 grp_text=c(dat$grp_text,dat$grp_text))
	
	pp <- ggplot(data=tt, aes(x=group, y=OR, ymin=OR025, ymax=OR975)) +
	      geom_pointrange(shape=25, size=0.35, aes(fill=grp_text, color=grp_text)) +
	      geom_hline(yintercept=1, linetype="dotted", color="grey", size=0.3) + 
	      scale_x_discrete("",labels=c("Finland"="Finland\n\n","Sweden"="Sweden\n\n")) + 
	      theme_classic() + 
	      theme(axis.text.x=element_text(size=8, face="bold"), axis.title.y=element_text(size=8, face="bold"), plot.title=element_text(size=9, face="bold")) + 
	      labs(y="Odds Ratio of being childless", x="", title="Relationship of disease diagnoses\nwith childlessness in Finland and Sweden") + 
	      # ylim(ylims)
	      scale_y_log10(limits=c(0.4,25),breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25),labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25)) + 
	      scale_color_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
	      scale_fill_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) + 
	      theme(legend.position="none")
	
	return(pp)
	}
}



# Age-specific effect
age_onset<- function(dat){
	pp <- ggplot(data=dat, aes(x=Age_onset_grp, y=OR.meta, ymin=OR025.meta, ymax=OR975.meta, group=1)) +
		# geom_pointrange() +
		geom_pointrange(shape=25, size=0.35, aes(fill=grp_text, color=grp_text)) + 
		geom_line(aes(color=grp_text), linetype="dashed", size=0.3) + 
		scale_color_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
		scale_fill_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
		geom_hline(yintercept = 1, linetype="dotted", color = "grey", size=0.3) +
		theme_classic() + 
		theme(axis.text.x=element_text(size=8, face="bold"), axis.title=element_text(size=8, face="bold"), plot.title=element_text(size=9, face="bold")) + 
		labs(y="Odds Ratio of being childless", x="Age of onset", title="Age-of-onset-dependent effects") + 
		theme(plot.title = element_text(size=10)) + 
		theme(legend.position="none")
	return(pp)
}



## change to mediation analysis (total effect, indirect effect, direct effect)
compare_mediation <- function(dat, ylims, title_sex) {
	dat$EffectType <- factor(dat$EffectType, levels=c("Total","Indirect","Direct"))
	data_comp_med_f <<- dat    # save data used for plot. needed for hovering
	
	pp <- ggplot(data=dat, aes(x=EffectType, y=OR.meta, ymin=OR025.meta, ymax=OR975.meta)) +
	      geom_pointrange(shape=25, size=0.35, aes(fill=grp_text, color=grp_text)) +
	      geom_hline(yintercept=1, linetype="dotted", color = "grey", size=0.3) + 
	      theme_classic() + 
	      labs(y="Odds Ratio of being childless", x="", title="Mediation effect by singlehood") + 
	      # ylim(ylims)
	      scale_y_log10(limits=c(0.4,25),breaks=c(0.5,1,2,3,4,6,8,10,12,15,20,25),labels=c(0.5,1,2,3,4,6,8,10,12,15,20,25)) + 
	      theme(axis.text.x=element_text(size=8, face="bold"), axis.title.y=element_text(size=8, face="bold"), plot.title=element_text(size=9, face="bold")) + 
	      scale_color_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) +
	      scale_fill_manual(breaks=mycolors[,"grp_text"], values=mycolors[,"mycol"]) + 
	      theme(legend.position="none")
	
	return(pp)
}



get_style <- function(hover, left_lim, top_lim){
  # point position
  left_px <- hover$coords_css$x
  top_px <- hover$coords_css$y
  
  # to make sure that box doesn't go out of the plot (too far right or down)
  left_px <- min(left_px, left_lim)
  top_px <- min(top_px, top_lim)
  
  # create style property for tooltip
  # background color is set so tooltip is a bit transparent
  # z-index is set so we are sure are tooltip will be on top
  
  style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                  "left:", left_px, "px; top:", top_px, "px;", "pointer-events: none;",
                  "font-size:10px;", "color: black;", 
                  # "margin-bottom: 0px; border-bottom: 0px;",  #has no effect
                  "padding-top: 3px; padding-bottom: 0px; padding-left: 3px; padding-right: 3px;"
  )
  return(style)
}



#### FUNCTIONS TO DEFINE TEXT ####
endpoint_description <- function(dat){
	return(c(paste("Name: ",dat[1,"LONGNAME"]), # select first row from two/four
		paste("ICD8: ",dat[1, "ICD8"]),
		paste("ICD9: ",dat[1, "ICD9"]),
		paste("ICD10: ",dat[1, "ICD10"])
	))
}



get_url <- function(endpoint){
  #endpoint <- input$choicepan
  # replace characters that are not goog for URL
  # reserved or unwise characters
  # NOTE, replacement not done!
  
  # NOTE, \ is not replaced
  # endpoint <- gsub(pattern = "[\]", replacement =  "%22", endpoint, fixed = TRUE) #not working
  # endpoint <- gsub(pattern = ";", replacement =  "%22", endpoint,  fixed = TRUE)
  # endpoint <- gsub(pattern = "/", replacement =  "%22", endpoint,  fixed = TRUE)
  # endpoint <- gsub(pattern = "?", replacement =  "%22", endpoint,  fixed = TRUE)
  # endpoint <- gsub(pattern = ":", replacement =  "%22", endpoint,  fixed = TRUE)
  # endpoint <- gsub(pattern = "@", replacement =  "%22", endpoint,  fixed = TRUE)
  # endpoint <- gsub(pattern = "&", replacement =  "%22", endpoint,  fixed = TRUE)
  # endpoint <- gsub(pattern = "=", replacement =  "%22", endpoint,  fixed = TRUE)
  # #endpoint <- gsub(pattern = "+", replacement =  "%22", endpoint,  fixed = TRUE)
  # endpoint <- gsub(pattern = "$", replacement =  "%22", endpoint,  fixed = TRUE)
  # endpoint <- gsub(pattern = ",", replacement =  "%22", endpoint,  fixed = TRUE)
  # endpoint <- gsub(pattern = "%", replacement =  "%22", endpoint,  fixed = TRUE)
  # endpoint <- gsub(pattern = "#", replacement =  "%22", endpoint,  fixed = TRUE)
  # endpoint <- gsub(pattern = "{", replacement =  "%22", endpoint,  fixed = TRUE)
  # endpoint <- gsub(pattern = "}", replacement =  "%22", endpoint,  fixed = TRUE)
  # endpoint <- gsub(pattern = "|", replacement =  "%22", endpoint,  fixed = TRUE)
  # endpoint <- gsub(pattern = "^", replacement =  "%22", endpoint,  fixed = TRUE)
  # endpoint <- gsub(pattern = "[", replacement =  "%22", endpoint,  fixed = TRUE)
  # endpoint <- gsub(pattern = "]", replacement =  "%22", endpoint,  fixed = TRUE)
  # endpoint <- gsub(pattern = "`", replacement =  "%22", endpoint,  fixed = TRUE)
  
  # + is a reserved character -> replace space by + AFTER replacing reserved characters
  endpoint <- gsub(pattern = " ", replacement =  "+", endpoint)
  
  # output$link <- renderUI({
  #   tagList(paste0("http://www.google.com/search?q=", endpoint))
  # })
  url_href <- paste0("http://www.google.com/search?q=", endpoint)
  
  url <- a("(Search disease in Google)", href= url_href, target= "_blank")
  
  return(url)
}




###############################################################################################
####                                   LOAD DATA                                           ####
###############################################################################################

# index_meta 
load("index_meta_m_out.Rdata")
load("index_meta_f_out.Rdata")

# moddel comparions for has_child
load("index_meta_m_mod.Rdata")
load("index_meta_f_mod.Rdata")

# cond results for LRS on having a child for the main plot
load("m_cond_haschild.Rdata")
load("f_cond_haschild.Rdata")

# mediation analysis 
load("med_meta_m.Rdata")
load("med_meta_f.Rdata")

# index_age_onset
load("index_age_onset_m.Rdata") 
load("index_age_onset_f.Rdata")

# results of conditional anlysis, for population compare
load("index_meta_m_pop.Rdata")
load("index_meta_f_pop.Rdata")

# sex specific disease names
load("names_dis_m.Rdata")
load("names_dis_f.Rdata")

# all diseases, grouped by disease type
load("grouped_list.Rdata")

# prevalence
# NOTE, a typo fixed in one long name afterwards -> remember not to use long names for matching
load("prevalence.Rdata")



###############################################################################################
####                                       UI                                              ####
###############################################################################################

ui <- fluidPage(
	
	# Include custom CSS to prevent horizontal scrollbar to appear from splitLayout
	tags$head(
		tags$style(HTML('.shiny-split-layout>div {overflow: hidden;}')),
	),
	
	theme="simplex.min.css",
	tags$style(
		type="text/css",
		"label {font-size: 12px;}",
		".recalculating {opacity: 1.0;}"
	),
	
	
	## zero line: Endpoint information --------------------------------------------------
	# tags$h2("The relationship of major diseases with childlessness"),
	# tags$h4("   –– A sibling matched case-control and population register study in Finland and Sweden", style="white-space: pre-wrap"),
	tags$h2("Evidence from Finland and Sweden on the relationship between early-life diseases and lifetime childlessness in men and women"),
	tags$h4("https://doi.org/10.1038/s41562-023-01763", style="white-space: pre-wrap"),
	hr(),

	fluidRow(
		column(
			width = 4, 
			pickerInput(inputId="choicepan",
			label="Choose a disease",
			selected="F5_SCHZPHR",
			choices=grouped_list,   # list of disease from both sexes
			options = list(`live-search`=TRUE),
			multiple = F)
		),
#		column(
#			width = 3, radioButtons(inputId = "selected_plot",
#			label = "Switch main plot",
#			choiceNames = list("Impact of disease on number of children", "Impact of disease on having children"),
#			choiceValues = list("number", "having"))
#		),
		column(
			width = 8, 
			uiOutput("endpoint_info")
		)
	),
	
	fluidRow(
		column(
			width = 12, 
			tableOutput("n_table")
		)
	),
	
	
	## first line plots --------------------------------------------------
	hr(),
	fluidRow(
		column(
			width = 6, 
			tags$h3("Men")
		),
		column(
			width = 6, 
			tags$h3("Women")
		)
	),
	
	fluidRow(
		column(
			width = 6, 
			style = 'padding:2px',
			plotOutput("m_distPlot", height = "400px",click = "m_plot_click")
		),
		column(
			width = 6, 
			style = 'padding:2px',
			plotOutput("f_distPlot", height = "400px",click = "f_plot_click")
		)
	),
	
	fluidRow(
		column(
			width = 6,
			uiOutput("m_plot_info")
		),
		column(
			width = 6,
			uiOutput("f_plot_info")
		)
	),
	
	fluidRow(
		column(
			width = 6, 
			uiOutput("m_missing_plot_info")
		),
		column(
			width = 6, 
			uiOutput("f_missing_plot_info")
		)
	),
	
	
	## second line plots --------------------------------------------------
	splitLayout(    # use splitLayout to divide page to two
		# second line plots for males 
		fluidRow(
			style='padding-right:25px',   # align with main plots
			tabsetPanel(
				id = "m_conditional_plots",
				type = "hidden",
				tabPanel("two_plots",
					fluidRow(
							# NOTE! fluidRow has 12-unit wide grid even though
							# space is limited to half page by splitLayout 
							column(
								width = 6, 
								style='padding-left:25px; padding-right:0px',
								plotOutput(
									"m_plot_comp_m_2p", 
									height ="220px",
									hover = hoverOpts("hover_m_plot_comp_m_2p", delay = 100)
								),
								uiOutput("m_plot_comp_m_2p_info")
							),
							column(
								width = 6, 
								style='padding-left:1px; padding-right:4px',
								plotOutput(
									"m_plot_comp_mod_2p", 
									height ="225px",
									hover = hoverOpts("hover_m_plot_comp_mod_2p", delay = 100)
								),
								uiOutput("m_plot_comp_mod_2p_info")
							)
						)
					),
				tabPanel("three_plots",
					fluidRow(
						column(
							width = 4,
							style='padding-left:25px;padding-right:0px',
							plotOutput(
								"m_plot_comp_m_3p", 
								height ="220px", 
								hover = hoverOpts("hover_m_plot_comp_m_3p", delay = 100)
							),
							uiOutput("m_plot_comp_m_3p_info")
						),
						column(
							width = 4, 
							style='padding-left:1px;padding-right:0px',
							plotOutput(
								"m_plot_comp_mod_3p", 
								height ="225px",
								hover = hoverOpts("hover_m_plot_comp_mod_3p", delay = 100)
							),
							uiOutput("m_plot_comp_mod_3p_info")
						),
						column(
							width = 4, 
							style='padding-left:1px;padding-right:4px',
							plotOutput(
								"m_plot_comp_c_3p", 
								height ="220px",
								hover = hoverOpts("hover_m_plot_comp_c_3p",
								delay = 100)
							),
							uiOutput("m_plot_comp_c_3p_info")
						)
					)
				)
				
			)
		),
		
		# second line plots for females
		fluidRow(
			# style='padding-left:20px',   # align with main plots
			tabsetPanel(
				id = "f_conditional_plots",
				type = "hidden",
				tabPanel("two_plots",
					fluidRow(
						column(
							width = 6, 
							style='padding-left:25px; padding-right:0px;',
							plotOutput(
								"f_plot_comp_m_2p", 
								height ="220px",
								hover = hoverOpts("hover_f_plot_comp_m_2p", delay = 100)
							),
							uiOutput("f_plot_comp_m_2p_info")
						),
						column(
							width = 6, 
							style='padding-left:1px; padding-right:4px',
							plotOutput(
								"f_plot_comp_mod_2p", 
								height ="225px",
								hover = hoverOpts("hover_f_plot_comp_mod_2p", delay = 100)
							),
							uiOutput("f_plot_comp_mod_2p_info")
						)
					)
				),
				tabPanel("three_plots",
					fluidRow(
						column(
							width = 4, 
							style='padding-left:25px;padding-right:0px',
							plotOutput("f_plot_comp_m_3p", 
								height ="220px",
								hover = hoverOpts("hover_f_plot_comp_m_3p", delay = 100)
							),
							uiOutput("f_plot_comp_m_3p_info")
						),
						column(
							width = 4, 
							style='padding-left:1px;padding-right:0px', 
							plotOutput(
								"f_plot_comp_mod_3p", 
								height ="225px",
								hover = hoverOpts("hover_f_plot_comp_mod_3p", delay = 100)
							),
							uiOutput("f_plot_comp_mod_3p_info")
						),
						column(
							width = 4,
							style='padding-left:1px;padding-right:4px',
							plotOutput(
								"f_plot_comp_c_3p", 
								height ="220px",
								hover = hoverOpts("hover_f_plot_comp_c_3p", delay = 100)
							),
							uiOutput("f_plot_comp_c_3p_info")
						)
					)
				)
			)
		)
	),
	
	
	## third line plots --------------------------------------------------
	# age onset males, malf_pret, age onset females
	fluidRow(
		column(
			width = 3, 
			style = 'padding:2px',
			plotOutput(
				"m_plot_comp2", 
				height = "200px",
				hover = hoverOpts(id = "hover_m_plot_comp2", delay = 100)
			),
			uiOutput("hover_m_plot_comp2_info")
		),
		column(
			width = 3, 
			style = 'padding:2px',
			plotOutput(
				"m_plot_comp_med", 
				height = "200px",
				hover = hoverOpts(id = "hover_m_plot_comp_med", delay = 100)
			),
			uiOutput("m_plot_comp_med_info")
		),
		column(
			width = 3,
			style = 'padding:2px',
			plotOutput(
				"f_plot_comp2", 
				height = "200px", 
				hover = hoverOpts(id = "hover_f_plot_comp2", delay = 100)
			),
			uiOutput("hover_f_plot_comp2_info")
		),
		column(
			width = 3, 
			style = 'padding:2px',
			plotOutput(
				"f_plot_comp_med", 
				height = "200px",
				hover = hoverOpts(id = "hover_f_plot_comp_med", delay = 100)
			),
			uiOutput("f_plot_comp_med_info")
		)
	),
	
	fluidRow(
		column(
			width = 3, 
			uiOutput("m_age_onset_info")
		),
		column(
			width = 3, 
			uiOutput("m_med_info")
		),
		column(
			width = 3,
			uiOutput("f_age_onset_info")
		),
		column(
			width = 3, 
			uiOutput("f_med_info")
		)
	)
	
)




###############################################################################################
####                                     SERVER                                            ####
###############################################################################################

server <- function(session, input, output) {
	
	## ALLOWS TO CLICK ON THE PLOT AND GET AN UPDATED ENDPOINT, # update Endpoint input only if point in the plot is clicked, not empty space
	observeEvent(input$m_plot_click, {   # reset text underneath the other plot
		output$f_plot_info <- renderUI({tags$p("") })
		
		# selected_point is a row from the data that was used for making the plot, used data is saved in reactiveValues based on which main plot was selected
		selected_point = nearPoints(
			rv_plot_data$m_plot_data, 
			yvar="OR.meta", 
			xvar="seq", 
			input$m_plot_click, 
			addDist = FALSE,
			maxpoints=1
		)
		if (nrow(selected_point)==1){
			# Update selected endpoint  
			updatePickerInput(session, "choicepan", selected = selected_point$Endpoint)
			# reset text underneath the plot
			output$m_plot_info <- renderUI({tags$p("") })   
		} else if(nrow(selected_point)==0){
			# update text box underneath the plot 
			output$m_plot_info <- renderUI({tags$p("Nothing was selected. Please select a point from a plot.")})
		}
	})
	
	observeEvent(input$f_plot_click, {   # reset text underneath the other plot
		output$m_plot_info <- renderUI({tags$p("") })
		
		# used data is saved in reactiveValues based on which main plot was selected
		selected_point <- nearPoints(
			rv_plot_data$f_plot_data,
			yvar="OR.meta", 
			xvar="seq", 
			input$f_plot_click, 
			addDist = FALSE,
			maxpoints=1
		)
		if (nrow(selected_point)==1) {
			# Update selected endpoint 
			updatePickerInput(session, "choicepan", selected = selected_point$Endpoint)
			# reset text box underneath the plot
			output$f_plot_info <- renderUI({tags$p("") })
		} else if (nrow(selected_point)==0){
			# update text box underneath the plot 
			output$f_plot_info <- renderUI({tags$p("Nothing was selected. Please select a point from a plot.")})
		}
	})
	
	
	### DEFINE INPUT DATASETS ###
	# Index meta, gee
	# index_meta_m_endpoint <- reactive({index_meta_m_gee %>% filter(Endpoint %in% input$choicepan)})
	# index_meta_f_endpoint <- reactive({index_meta_f_gee %>% filter(Endpoint %in% input$choicepan)})
	index_meta_m_out_endpoint <- reactive({index_meta_m_out %>% filter(Endpoint %in% input$choicepan)})
	index_meta_f_out_endpoint <- reactive({index_meta_f_out %>% filter(Endpoint %in% input$choicepan)})

	
	index_meta_m_mod_endpoint <- reactive({index_meta_m_mod %>% filter(Endpoint %in% input$choicepan)})
	index_meta_f_mod_endpoint <- reactive({index_meta_f_mod %>% filter(Endpoint %in% input$choicepan)})
	
	# Index meta, gee, LRS == "has_child_Age4550"
	m_cond_haschild_endpoint <- reactive({m_cond_haschild %>% filter(Endpoint%in% input$choicepan)})
	f_cond_haschild_endpoint <- reactive({f_cond_haschild %>% filter(Endpoint%in% input$choicepan)})
	
	# Index meta, conditional 
	# index_meta_m_cond_endpoint <- reactive({index_meta_m_cond %>% filter(Endpoint %in% input$choicepan)})
	# index_meta_f_cond_endpoint <- reactive({index_meta_f_cond %>% filter(Endpoint %in% input$choicepan)})

	index_meta_m_pop_endpoint <- reactive({index_meta_m_pop %>% filter(Endpoint %in% input$choicepan)})
	index_meta_f_pop_endpoint <- reactive({index_meta_f_pop %>% filter(Endpoint %in% input$choicepan)})

	# Mediation meta
	med_meta_m_endpoint <- reactive({med_meta_m %>% filter(Endpoint %in% input$choicepan)})
	med_meta_f_endpoint <- reactive({med_meta_f %>% filter(Endpoint %in% input$choicepan)})
	
	# Age at onset
	index_age_onset_m_endpoint <- reactive({index_age_onset_m %>% filter(Endpoint %in% input$choicepan)})
	index_age_onset_f_endpoint <- reactive({index_age_onset_f %>% filter(Endpoint %in% input$choicepan)})
	
#	# Age at onset stats
#	index_age_onset_stats_m_endpoint <- reactive({index_age_onset_stats_m %>% filter(Endpoint %in% input$choicepan)})
#	index_age_onset_stats_f_endpoint <- reactive({index_age_onset_stats_f %>% filter(Endpoint %in% input$choicepan)})
	
	# prevalence
	prevalence_endpoint <- reactive({prevalence %>% filter(Endpoint %in% input$choicepan)})
	
	
	#### PLOTS ####
	## Main plots - conditional analysis ##
	# reactive values for saving used data
	rv_plot_data <- reactiveValues(f_plot_data = NULL, m_plot_data = NULL)
	
	main_plot_m <- mhd_index(dat=m_cond_haschild, names_dis_m, "Relationship of disease diagnoses with childlessness", "Odds ratio of being childless by age 50")
	main_plot_f <- mhd_index(dat=f_cond_haschild, names_dis_f, "Relationship of disease diagnoses with childlessness", "Odds ratio of being childless by age 45")
	
	
	output$m_distPlot <- renderPlot({main_plot_m + geom_pointrange(aes(x=seq, y=OR.meta, ymin=OR025.meta, ymax=OR975.meta), data=m_cond_haschild_endpoint(), shape=20, size=1, col="red")})
	output$f_distPlot <- renderPlot({main_plot_f + geom_pointrange(aes(x=seq, y=OR.meta, ymin=OR025.meta, ymax=OR975.meta), data=f_cond_haschild_endpoint(), shape=20, size=1, col="red")}) 
	rv_plot_data$f_plot_data <- f_cond_haschild
	rv_plot_data$m_plot_data <- m_cond_haschild
	
	
	## Second line plots ##
	# make reactive values. 
	# output plots change always when endpoint is changed, so no need to reset these plots
	rv <- reactiveValues(f_plot_comp_m = NULL, f_plot_comp_mod = NULL, f_plot_comp_c = NULL, f_plot_comp_med = NULL, 
	                     m_plot_comp_m = NULL, m_plot_comp_mod = NULL, m_plot_comp_c = NULL, m_plot_comp_med = NULL, ylims = NULL)
	
	
	#### observe change of selected endpoint ####
	observeEvent(input$choicepan, {
		
		#### TEXT ABOVE MAIN PLOT WITH ENDPOINT INFO ####
		# select data for endpoint description, 
		# use conditional analysis results data because that has results for non-significant results too
		# use index_meta_m_pop_endpoint if endpoint is available only for men
		# -> description available for all endpoints
		info_vector <- endpoint_description(prevalence_endpoint())
		
		# use long name of endpoint 
		url <- get_url(prevalence_endpoint()[1,"LONGNAME"]) 
		output$endpoint_info <- renderUI({
			tagList(
				tags$p(info_vector[1]),
				tags$p(url),   # link to Google search
				tags$p(info_vector[2]),
				tags$p(info_vector[3]),
				tags$p(info_vector[4])
			)
		})
		
		
		# prevalence and number of cases table 
		# table_data <- as.data.frame(prevalence_endpoint()[, c("Prevalence","N","N_case","N_Age_0_16","N_Age_17_24","N_Age_25_29","N_Age_30_34","N_Age_35_39","N_Age_40_44","N_Age_45_50")])
		table_data <- as.data.frame(prevalence_endpoint()[, c("Prevalence","N","N_Case","N_Age_0_15","N_Age_16_20","N_Age_21_25","N_Age_26_30","N_Age_31_35","N_Age_36_40","N_Age_41_45","N_Age_46_50")])
		# table_data <- as.data.frame(prevalence_endpoint())
		
		# prevent small prevalences to be rounded to zero 
		table_data$Prevalence <- as.character(table_data$Prevalence)
		row.names(table_data) <- prevalence_endpoint()$table_rows
		output$n_table <- renderTable(table_data, rownames = TRUE)
		
		#### reset text boxes when selected endpoint changes ####
		output$f_missing_plot_info <- renderUI({tags$p("")})
		output$m_missing_plot_info <- renderUI({tags$p("")})
		
		##### reactive second line PLOTS #####
		# show plots only if results are available for that sex
		
		#### Second line plots for females #### ------------------------
		# reset ylims every time when endpoint is changed
		rv$ylims <- NULL
		
		if (!(input$choicepan %in% names_dis_f$Endpoint)){
			# results not for woman, don't make female plots
			output$f_missing_plot_info <- renderUI({tags$p("Selected endpoint is not available for women.")})
			output$f_plot_comp_m_2p <- NULL
			output$f_plot_comp_mod_2p <- NULL
			output$f_plot_comp_m_3p <- NULL
			output$f_plot_comp_mod_3p <- NULL
			output$f_plot_comp_c_3p <- NULL
			output$f_plot_comp2 <- NULL
			output$f_plot_comp_med <- NULL
			output$f_age_onset_info <- NULL
			output$f_med_info <- NULL
		} else{
			## make female plots ##
			# get limit values for y-axes 
			rv$ylims <- get_ylims(sex="female", index_meta_m_out_endpoint(), index_meta_m_mod_endpoint(), index_meta_m_pop_endpoint(), med_meta_m_endpoint(), 
			                                    index_meta_f_out_endpoint(), index_meta_f_mod_endpoint(), index_meta_f_pop_endpoint(), med_meta_f_endpoint())
			                                    
			# Plots comparison between different measures
			rv$f_plot_comp_m <- comp_measures(index_meta_f_out_endpoint(), rv$ylims, title_sex = "women")
			
			# Plots comparison between different models for has_child 
			rv$f_plot_comp_mod <- comp_models(index_meta_f_mod_endpoint(), rv$ylims, title_sex = "women")
			
			# Plots comparing Finnish and Swedish
			rv$f_plot_comp_c <- comp_fin_swe(index_meta_f_pop_endpoint(), rv$ylims)
			
			# Render 
			# if not Swedish results, rv$f_plot_comp_c is NULL, render only two plots and give message
			if (!(is.character(rv$f_plot_comp_c))){
				# All 3 plots are available
				updateTabsetPanel(inputId = "f_conditional_plots", selected="three_plots")
				output$f_plot_comp_m_3p <- renderPlot(rv$f_plot_comp_m)
				output$f_plot_comp_mod_3p <- renderPlot(rv$f_plot_comp_mod)
				output$f_plot_comp_c_3p <- renderPlot(rv$f_plot_comp_c)
				
				## pop-up info boxes ##
				# get the data that was used to make the plot -> needed for hovering 
				f_data_comp_m <- rv$f_plot_comp_m$data
				output$f_plot_comp_m_3p_info <- renderUI({
					hover <- input$hover_f_plot_comp_m_3p
					point <- nearPoints(f_data_comp_m, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
					if(nrow(point) == 0) return(NULL)
					wellPanel(
						style = get_style(hover, 170, 80),
						p(HTML(paste0("Odds ratio: ", round(point$OR.meta, 2),"<br/>", "95% CI: [", round(point$OR025.meta, 2), ", ", round(point$OR975.meta, 2), "]", "<br/>", "P-value: ", point$P_val.meta)))
					)
				})
				
				f_data_comp_mod <- rv$f_plot_comp_mod$data 
				output$f_plot_comp_mod_3p_info <- renderUI({
					hover <- input$hover_f_plot_comp_mod_3p
					point <- nearPoints(f_data_comp_mod, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
					if(nrow(point) == 0) return(NULL)
					wellPanel(
						style = get_style(hover, 170, 80),
						p(HTML(paste0("Odds ratio: ", round(point$OR.meta, 2),"<br/>", "95% CI: [", round(point$OR025.meta, 2), ", ", round(point$OR975.meta, 2), "]", "<br/>", "P-value: ", point$P_val.meta)))
					)
				})
				
				# get data. here availability of the data is already tested
				f_data_comp_c <- rv$f_plot_comp_c$data
				output$f_plot_comp_c_3p_info <- renderUI({
					hover <- input$hover_f_plot_comp_c_3p
					point <- nearPoints(f_data_comp_c, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
					if(nrow(point) == 0) return(NULL)
					wellPanel(
						style = get_style(hover, 170, 80), 
						p(HTML(paste0("Odds ratio: ", round(point$OR, 2),"<br/>", "95% CI: [", round(point$OR025, 2),  ", ", round(point$OR975, 2), "]", "<br/>","P-value: ", point$P_val)))
					)
				})
				
			} else {
				updateTabsetPanel(inputId = "f_conditional_plots", selected = "two_plots")
				output$f_plot_comp_m_2p <- renderPlot(rv$f_plot_comp_m)
				output$f_plot_comp_mod_2p <- renderPlot(rv$f_plot_comp_mod)
				output$f_missing_plot_info <- renderUI({tags$p(paste("There are no results from ", rv$f_plot_comp_c, ".", sep = ""))})
				
				# pop-up info boxes
				# get the data that was used to make the plot -> needed for hovering 
				f_data_comp_m <- rv$f_plot_comp_m$data
				output$f_plot_comp_m_2p_info <- renderUI({
					hover <- input$hover_f_plot_comp_m_2p
					point <- nearPoints(f_data_comp_m, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
					if(nrow(point) == 0) return(NULL)
					wellPanel(
						style = get_style(hover, 252, 80),
						p(HTML(paste0("Odds ratio: ", round(point$OR.meta, 2),"<br/>","95% CI: [", round(point$OR025.meta, 2),  ", ", round(point$OR975.meta, 2), "]", "<br/>","P-value: ", point$P_val.meta)))
					)
				})
				
				f_data_comp_mod <- rv$f_plot_comp_mod$data
				output$f_plot_comp_mod_2p_info <- renderUI({
					hover <- input$hover_f_plot_comp_mod_2p
					point <- nearPoints(f_data_comp_mod, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
					if(nrow(point) == 0) return(NULL)
					wellPanel(
						style = get_style(hover, 252, 80),
						p(HTML(paste0("Odds ratio: ", round(point$OR.meta, 2),"<br/>","95% CI: [", round(point$OR025.meta, 2),  ", ", round(point$OR975.meta, 2), "]", "<br/>","P-value: ", point$P_val.meta)))
					)
				})
				
			}
			
			### Third line plots for females ###  ------------------------
			## Age at onset ##
			if (nrow(index_age_onset_f_endpoint())==0){
				output$f_plot_comp2 <- NULL
				output$f_age_onset_info <- renderUI({ "There are no age of onset results for this endpoint." })
			} else {
				# output$f_age_onset_info <- renderUI({ tags$p(paste0("Age at onset effect is: ",index_age_onset_stats_f_endpoint()$Polynomial)) })
				output$f_age_onset_info <- renderUI({tags$p("") })
				rv$f_plot_comp2 <- age_onset(index_age_onset_f_endpoint()) 
				output$f_plot_comp2 <- renderPlot(rv$f_plot_comp2)
				
				# pop-up info boxes
				age_onset_data_f <- rv$f_plot_comp2$data   # get used data
				output$hover_f_plot_comp2_info <- renderUI({
					hover <- input$hover_f_plot_comp2
					point <- nearPoints(age_onset_data_f, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
					if(nrow(point) == 0) return(NULL)
					wellPanel(
						style = get_style(hover, 300, 80),
						# p(HTML(paste0("Odds ratio: ", round(point$OR, 2),"<br/>", "95% CI: [", round(point$OR025, 2),  ", ", round(point$OR975, 2), "]", "<br/>","P-value: ", point$P_val)))
						p(HTML(paste0("Age group: ", point$Age_onset_grp, "<br/>", "<br/>","Odds ratio: ", round(point$OR.meta, 2),"<br/>", "95% CI: [", round(point$OR025.meta, 2),  ", ", round(point$OR975.meta, 2), "]", "<br/>","P-value: ", point$P_val.meta)))
					)
				})
			}
			
			## Mediation analysis 
			if (nrow(med_meta_f_endpoint())==0){ 
				output$f_plot_comp_med <- NULL
				output$f_med_info <- renderUI({tags$p("There are no mediation analyses for this endpoint.") })
			} else{
				rv$f_plot_comp_med <- compare_mediation(med_meta_f_endpoint(), rv$ylims, title_sex="women") 
				output$f_plot_comp_med <- renderPlot(rv$f_plot_comp_med) 
				output$f_med_info <- renderUI({tags$p("") })
				
				# pop-up info boxes
				f_data_comp_med <- rv$f_plot_comp_med$data   # get used data
				output$f_plot_comp_med_info <- renderUI({
					hover <- input$hover_f_plot_comp_med
					point <- nearPoints(f_data_comp_med, hover, threshold=5, maxpoints=1, addDist=TRUE)
					if(nrow(point) == 0) return(NULL)
					wellPanel(
						style = get_style(hover, 300, 80),
						# p(HTML(paste0("Odds ratio: ", round(point$OR, 2),"<br/>", "95% CI: [", round(point$OR025, 2),  ", ", round(point$OR975, 2), "]", "<br/>","P-value: ", point$P_val)))
						p(HTML(paste0("Odds ratio: ", round(point$OR.meta, 2),"<br/>", "95% CI: [", round(point$OR025.meta, 2),  ", ", round(point$OR975.meta, 2), "]", "<br/>","P-value: ", point$P_val.meta)))
					)
				})
			}

			
			
			
		}
		
		
		### Second line plots for males ###  ------------------------
		if (!(input$choicepan %in% names_dis_m$Endpoint)){
			# results not for males, don't make male plots
			output$m_missing_plot_info <- renderUI({tags$p("Selected endpoint is not available for men.") })
			output$m_plot_comp_m_2p <- NULL
			output$m_plot_comp_mod_2p <- NULL
			output$m_plot_comp_m_3p <- NULL
			output$m_plot_comp_mod_3p <- NULL
			output$m_plot_comp_c_3p <- NULL
			output$m_plot_comp2 <- NULL
			output$m_plot_comp_med <- NULL
			output$m_age_onset_info <- NULL
			output$m_med_info <- NULL
		} else{
			## make male plots ##
			# get limit values for y-axes, if not already calculated when doing female plots
			if (is.null(rv$ylims)){
				rv$ylims <- get_ylims(sex="male", index_meta_m_out_endpoint(), index_meta_m_mod_endpoint(), index_meta_m_pop_endpoint(), med_meta_m_endpoint(), 
				                                  index_meta_f_out_endpoint(), index_meta_f_mod_endpoint(), index_meta_f_pop_endpoint(), med_meta_f_endpoint())
			}
			
			# Plots comparison between different measures
			# rv$m_plot_comp_m <- comp_measures(index_meta_m_endpoint(), rv$ylims, title_sex = "men")
			rv$m_plot_comp_m <- comp_measures(index_meta_m_out_endpoint(), rv$ylims, title_sex = "men")

			# Plots comparison between different models for has_child
			rv$m_plot_comp_mod <- comp_models(index_meta_m_mod_endpoint(), rv$ylims, title_sex = "men")
			
			# Plots comparing Finnish and Swedish
			rv$m_plot_comp_c <- comp_fin_swe(index_meta_m_pop_endpoint(), rv$ylims)
			
			# Render
			# if not Swedish results, rv$m_plot_comp_c is NULL, render only two plots and give message
			if (!(is.character(rv$m_plot_comp_c))){
				updateTabsetPanel(inputId = "m_conditional_plots", selected = "three_plots")
				output$m_plot_comp_m_3p <- renderPlot(rv$m_plot_comp_m)
				output$m_plot_comp_mod_3p <- renderPlot(rv$m_plot_comp_mod)
				output$m_plot_comp_c_3p <- renderPlot(rv$m_plot_comp_c)
				
				## pop-up info boxes ##
				m_data_comp_m <- rv$m_plot_comp_m$data    # get used data
				output$m_plot_comp_m_3p_info <- renderUI({
					hover <- input$hover_m_plot_comp_m_3p
					point <- nearPoints(m_data_comp_m, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
					if(nrow(point) == 0) return(NULL)
					wellPanel(
						style = get_style(hover, 170, 80),
						p(HTML(paste0("Odds ratio: ", round(point$OR.meta, 2),"<br/>","95% CI: [", round(point$OR025.meta, 2),  ", ", round(point$OR975.meta, 2), "]", "<br/>","P-value: ", point$P_val.meta)))
					)
				})
				
				m_data_comp_mod <- rv$m_plot_comp_mod$data
				output$m_plot_comp_mod_3p_info <- renderUI({
					hover <- input$hover_m_plot_comp_mod_3p
					point <- nearPoints(m_data_comp_mod, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
					if(nrow(point) == 0) return(NULL)
					wellPanel(
						style = get_style(hover, 170, 80),
						p(HTML(paste0("Odds ratio: ", round(point$OR.meta, 2),"<br/>","95% CI: [", round(point$OR025.meta, 2),  ", ", round(point$OR975.meta, 2), "]", "<br/>","P-value: ", point$P_val.meta)))
					)
				})
				
				m_data_comp_c <- rv$m_plot_comp_c$data
				output$m_plot_comp_c_3p_info <- renderUI({
					hover <- input$hover_m_plot_comp_c_3p
					point <- nearPoints(m_data_comp_c, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
					if(nrow(point) == 0) return(NULL)
					wellPanel(
						style = get_style(hover, 170, 80),
						p(HTML(paste0("Odds ratio: ", round(point$OR, 2),"<br/>","95% CI: [", round(point$OR025, 2),  ", ", round(point$OR975, 2), "]", "<br/>","P-value: ", point$P_val)))
					)
				})
				
			} else {
				output$m_missing_plot_info <- renderUI({tags$p(paste("There are no results from ", rv$m_plot_comp_c, ".", sep = "")) })
				updateTabsetPanel(inputId = "m_conditional_plots", selected = "two_plots")
				output$m_plot_comp_m_2p <- renderPlot(rv$m_plot_comp_m)
				output$m_plot_comp_mod_2p <- renderPlot(rv$m_plot_comp_mod)
				
				## pop-up info boxes ##
				m_data_comp_m <- rv$m_plot_comp_m$data    # get used data
				output$m_plot_comp_m_2p_info <- renderUI({
					hover <- input$hover_m_plot_comp_m_2p
					point <- nearPoints(m_data_comp_m, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
					if(nrow(point) == 0) return(NULL)
					wellPanel(
						style = get_style(hover, 252, 80),
						p(HTML(paste0("Odds ratio: ", round(point$OR.meta, 2),"<br/>","95% CI: [", round(point$OR025.meta, 2),  ", ", round(point$OR975.meta, 2), "]", "<br/>","P-value: ", point$P_val.meta)))
					)
				})
				
				m_data_comp_mod <- rv$m_plot_comp_mod$data
				output$m_plot_comp_mod_2p_info <- renderUI({
					hover <- input$hover_m_plot_comp_mod_2p
					point <- nearPoints(m_data_comp_mod, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
					if(nrow(point) == 0) return(NULL)
					wellPanel(
						style = get_style(hover, 252, 80),
						p(HTML(paste0("Odds ratio: ", round(point$OR.meta, 2),"<br/>","95% CI: [", round(point$OR025.meta, 2),  ", ", round(point$OR975.meta, 2), "]", "<br/>","P-value: ", point$P_val.meta)))
					)
				})
			}
			
			
			### Third line plots for males ### ------------------------
			## Age at onset plot ##
			if (nrow(index_age_onset_m_endpoint())==0){
				output$m_plot_comp2 <- NULL
				output$m_age_onset_info <- renderUI({ tags$p("There are no age of onset results for this endpoint.") })
			} else {
				# output$m_age_onset_info <- renderUI({ tags$p(paste0("Age at onset effect is: ",index_age_onset_stats_m_endpoint()$Polynomial)) })
				output$m_age_onset_info <- renderUI({tags$p("") })
				rv$m_plot_comp2 <- age_onset(index_age_onset_m_endpoint()) 
				output$m_plot_comp2 <- renderPlot(rv$m_plot_comp2)
				
				# pop-up info boxes
				age_onset_data_m <- rv$m_plot_comp2$data   # get used data
				output$hover_m_plot_comp2_info <- renderUI({
					hover <- input$hover_m_plot_comp2
					point <- nearPoints(age_onset_data_m, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
					if(nrow(point) == 0) return(NULL)
					wellPanel(
						style = get_style(hover, 300, 80),
						p(HTML(paste0("Age group: ", point$Age_onset_grp, "<br/>", "<br/>","Odds ratio: ", round(point$OR.meta, 2),"<br/>", "95% CI: [", round(point$OR025.meta, 2),  ", ", round(point$OR975.meta, 2), "]", "<br/>","P-value: ", point$P_val.meta)))
					)
				})
			}
			
			## Mediation analysis 
			if (nrow(med_meta_m_endpoint())==0){ 
				output$m_plot_comp_med <- NULL
				output$m_med_info <- renderUI({tags$p("There are no mediation analyses for this endpoint.") })
			} else{
				rv$m_plot_comp_med <- compare_mediation(med_meta_m_endpoint(), rv$ylims, title_sex="men") 
				output$m_plot_comp_med <- renderPlot(rv$m_plot_comp_med) 
				output$m_med_info <- renderUI({tags$p("") })
				
				# pop-up info boxes
				m_data_comp_med <- rv$m_plot_comp_med$data   # get used data
				output$m_plot_comp_med_info <- renderUI({
					hover <- input$hover_m_plot_comp_med
					point <- nearPoints(m_data_comp_med, hover, threshold=5, maxpoints=1, addDist=TRUE)
					if(nrow(point) == 0) return(NULL)
					wellPanel(
						style = get_style(hover, 300, 80),
						# p(HTML(paste0("Odds ratio: ", round(point$OR, 2),"<br/>", "95% CI: [", round(point$OR025, 2),  ", ", round(point$OR975, 2), "]", "<br/>","P-value: ", point$P_val)))
						p(HTML(paste0("Odds ratio: ", round(point$OR.meta, 2),"<br/>", "95% CI: [", round(point$OR025.meta, 2),  ", ", round(point$OR975.meta, 2), "]", "<br/>","P-value: ", point$P_val.meta)))
					)
				})
			}

			
		}
		
	})   # End of observeEvent for choicepan
}

shinyApp(ui = ui, server = server)

