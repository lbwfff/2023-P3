#####################################################################
#先批量把单因素做了，然后再多因素
#结果：ukb-b-6353
library(TwoSampleMR)
library(MVMR)
library(tidyr)
library(tibble)
library(dplyr)

source('mvmr.R')
load('available_outcomes.RData')
filao<-ao[ao$population=='European',]

id_exposure<-c('ukb-b-5192','ukb-a-6','ukb-b-4943','ukb-b-969',
                'ukb-b-4094','ukb-b-6811' )  #感觉驾驶和其他几个因子联系有点过于弱了，所以没加到方案里，暂时先不用


{
exposure_dat <- extract_instruments(id_exposure,
                                    p1 = 5e-08)

outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ukb-b-6353")

dat <- harmonise_data(
  exposure_dat = exposure_dat,
  outcome_dat = outcome_dat)

res <- mr(dat)

result_2smr <- res %>%
  split_outcome() %>%
  separate(outcome, "outcome", sep="[(]") %>% 
  mutate(outcome=stringr::str_trim(outcome))%>% 
  generate_odds_ratios() %>% 
  select(-id.exposure, -id.outcome) %>% 
  tidy_pvals()

f <-function(x) unlist(strsplit(x['exposure'],'[:]'))[2]
result_2smr$id<-apply(result_2smr,1,f)
cache<-filao[match(result_2smr$id,filao$id),]
result_2smr$trait<-cache$trait

} 

sigmr<-function(name,inid,outid){
  
  path<-paste0('./',name,'/')
  dir.create(path)
  
  library(TwoSampleMR)
  
  exposure_dat <- extract_instruments(inid,p1 = 5e-08)
  write.table(exposure_dat, 
              file = paste0(path,"01.exposure_dat_input.xls"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
  
  outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=outid)
  write.table(outcome_dat, 
              file = paste0(path,"02.outcome_dat_input.xls"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
  
  dat <- harmonise_data(
    exposure_dat = exposure_dat,
    outcome_dat = outcome_dat)
  write.table(dat, 
              file = paste0(path,"03.Harmonize_effects.xls"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
  
  df_MR_hetero <- mr_heterogeneity(dat)
  write.table(df_MR_hetero, 
              file = paste0(path,"05.MR_heterogeneity.xls"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
  
  df_MR_pleio <- mr_pleiotropy_test(dat)
  write.table(df_MR_pleio, 
              file = paste0(path,"06.Horizontal_pleiotropy.xls"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
  
  
  
  if (df_MR_hetero$Q_pval[nrow(df_MR_hetero)]<0.05) {
    res=mr(dat, method_list = c("mr_ivw_mre",
                                "mr_egger_regression","mr_weighted_median", 
                                'mr_simple_mode','mr_weighted_mode'))} else{
                                  res=mr(dat, method_list = c("mr_ivw_fe",
                                                              "mr_egger_regression","mr_weighted_median", 
                                                              'mr_simple_mode','mr_weighted_mode'
                                                              ))}
  
  write.table(res, file = paste0(path,"04.MR_results.xls"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
  
  OR <-generate_odds_ratios(res)
  OR
  write.table(OR, 
              file = paste0(path,"04.MR_results_OR.xls"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
  
  res_loo <- mr_leaveoneout(dat)
  pdf(file = paste0(path,"07.Leave-one-out.pdf"), height = 12, width = 6)
  print(mr_leaveoneout_plot(res_loo))
  dev.off()
  
  pdf(file = paste0(path,"08.Scatter.pdf"), height = 8, width = 8)
  print(mr_scatter_plot(res, dat))
  dev.off()
  
  res_single <- mr_singlesnp(dat)
  write.table(df_MR_pleio, 
              file = paste0(path,"09.singlesnp_res.xls"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
  
  pdf(file = paste0(path,"10.forest.pdf"), height = 12, width = 6)
  print(mr_forest_plot(res_single))
  dev.off()
  
  pdf(file = paste0(path,"11.funnel.pdf"), height = 8, width = 8)
  print(mr_funnel_plot(res_single))
  dev.off()
}

ana<-result_2smr[!duplicated(result_2smr$id),]

for (i in 1:nrow(ana)){
  sigmr(ana$id[i],ana$id[i],'ukb-b-6353')
}

cache<-as.data.frame(array(NA,c(1,6)))
colnames(cache)<-c('trait','id','method','or','or_lci95','or_uci95')

for (i in 1:nrow(ana)){
  
  path<-paste0('./',ana$id[i],'/04.MR_results_OR.xls')
  or<-read.table(path,sep = '\t',header = T)
  
  cache2<-data.frame(trait=c(ana$trait[i]),
                     id=c(ana$id[i]),
                     method=c('MR Egger','Weighted median','Inverse variance weighted'),
                     or=c(or$or),
                     or_lci95=c(or$or_lci95),
                     or_uci95=c(or$or_uci95))
  
  cache<-rbind(cache,cache2)
  
}

plot<-cache[-1,]

library(ggplot2)
library(cowplot)
library(MetBrewer)

pdf("sig_3way.pdf",width = 6,height = 6)
ggplot(plot, aes(y=trait, x=or,  colour=method)) + 
  geom_errorbarh(aes(xmin=or_lci95, xmax=or_uci95), height=.3) +
  geom_point(size=2)+
  # xlim(0,2)+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 5,type ='discrete'))+
  geom_vline(xintercept=1, linetype=3) +
  theme_minimal_hgrid(10, rel_small = 1) +
  facet_wrap(~method, ncol=1)+
  labs(color = "",y = "", x = "Odds Ratio",
       title= paste0("Univariate MR results for ",
                     'Myopia' ,", 95% CI") )+
  theme(legend.position = "none")
dev.off()
#唯一的问题在于，看电视被算成了保护因素，就很奇怪

##############################################################
#然后是多因素的部分
id_exposure<-c('ukb-a-6','ukb-b-5192','ukb-b-969','ukb-b-4943')
exposure_dat <- mv_extract_exposures(id_exposure,find_proxies = FALSE,
                                     force_server = FALSE)
outcome_dat <- extract_outcome_data(exposure_dat$SNP, 'ukb-b-6353')

mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)

res <- mv_multiple(mvdat,plots=T)

result_2smr <- res$result %>%  #调整一下结果格式计算OR值
  split_outcome() %>%
  separate(outcome, "outcome", sep="[(]") %>% 
  mutate(outcome=stringr::str_trim(outcome))%>% 
  generate_odds_ratios() %>% 
  select(-id.exposure, -id.outcome) %>% 
  tidy_pvals()

write.csv(result_2smr,file = 'mvmr.csv')
result_2smr$type<-c('TwoSampleMR')

f <-function(x) unlist(strsplit(x['exposure'],'[|]'))[1]
result_2smr$name<-apply(result_2smr,1,f)
f <-function(x) substr(x['name'],1,nchar(x['name'])-1)
result_2smr$name<-apply(result_2smr,1,f)

pdf("multi_1way.pdf",width = 6,height = 2)
ggplot(result_2smr, aes(y=name, x=or, label=outcome, colour=type)) +  #shape=exposure, 
  geom_errorbarh(aes(xmin=or_lci95, xmax=or_uci95), height=.3) +
  geom_point(size=2)+
  # xlim(0,2)+
  scale_color_manual(values=met.brewer("VanGogh1", 3,type ='discrete'))+
  # scale_shape_manual(values = c(19,20,17)) +
  geom_vline(xintercept=1, linetype=3) +
  theme_minimal_hgrid(10, rel_small = 1) +
  # facet_wrap(~outcome, ncol=1)+
  labs(color = "",y = "", x = "Odds Ratio",
       title= paste0("Multivariate MR results for ",
                     'Myopia' ,", 95% CI") )+
  theme(legend.position = "none")
dev.off()
#

#然后在这里我想加一张单变量和多变量一起的图

sig<-plot[plot$method=='Inverse variance weighted',]
sig<-sig[sig$trait %in% result_2smr$name,]

sigandmul<-data.frame(name=c(sig$trait,result_2smr$name),
                      or=c(sig$or,result_2smr$or),
                      or_lci95=c(sig$or_lci95,result_2smr$or_lci95),
                      or_uci95=c(sig$or_uci95,result_2smr$or_uci95),
                      type=c(rep('Univariate',nrow(sig)),rep('Multivariate',nrow(result_2smr))))
sigandmul$type<-factor(sigandmul$type,levels = c('Univariate','Multivariate'))

pdf("sig_and_mul.pdf",width = 6,height = 4)
ggplot(sigandmul, aes(y=name, x=or, label=type, colour=type)) +  #shape=exposure, 
  geom_errorbarh(aes(xmin=or_lci95, xmax=or_uci95), height=.3) +
  geom_point(size=2)+
  # xlim(0,2)+
  scale_color_manual(values=met.brewer("VanGogh1", 3,type ='discrete'))+
  # scale_shape_manual(values = c(19,20,17)) +
  geom_vline(xintercept=1, linetype=3) +
  theme_minimal_hgrid(10, rel_small = 1) +
  facet_wrap(~type, ncol=1)+
  labs(color = "",y = "", x = "Odds Ratio",
       title= paste0("MR results for ",
                     'Myopia' ,", 95% CI") )+
  theme(legend.position = "none")
dev.off()

#然后是MVMR
mvmr_input <- make_mvmr_input(exposure_dat, outcome.id.mrbase= 'ukb-b-6353') #这个input好像也是用twosamplemr的function做的，只不过把结果调整为了MVMR的input形式

mvmr_out <- format_mvmr(BXGs = mvmr_input$XGs %>% select(contains("beta")),  
                        BYG = mvmr_input$YG$beta.outcome,                     
                        seBXGs = mvmr_input$XGs %>% select(contains("se")),  
                        seBYG = mvmr_input$YG$se.outcome,                     
                        RSID = mvmr_input$XGs$SNP)   

mvmr_res <- ivw_mvmr(r_input=mvmr_out) #多变量分析

result_mvmr <-
  mvmr_res %>% 
  tidy_mvmr_output() %>% 
  mutate(exposure = mvmr_input$exposures,
         outcome = 'Breast cancer') %>% 
  select(exposure, outcome, everything()) %>% 
  tidy_pvals() 

result_mvmr$type<-c('MVMR')

f <-function(x) unlist(strsplit(x['exposure'],'[|]'))[1]
result_mvmr$name<-apply(result_mvmr,1,f)
f <-function(x) substr(x['name'],1,nchar(x['name'])-1)
result_mvmr$name<-apply(result_mvmr,1,f)

merge<-data.frame(exposure=c(result_2smr$name,result_mvmr$name),
                  OR=c(result_2smr$or,result_mvmr$or),outcome=c('myopia'),
                  orlow=c(result_2smr$or_lci95,result_mvmr$or_lci95),
                  orup=c(result_2smr$or_uci95,result_mvmr$or_uci95),
                  type=c(result_2smr$type,result_mvmr$type))  #把twosampleMR的结果和MVMR整合到一起，两个output格式改得还不一样需要手动调一下

merge$type<-factor(merge$type,levels = c('TwoSampleMR','MVMR'))

pdf("multi_twoway.pdf",width = 6,height = 4)
ggplot(merge, aes(y=exposure, x=OR, label=outcome, colour=type)) +  #shape=exposure, 
  geom_errorbarh(aes(xmin=orlow, xmax=orup), height=.3) +
  geom_point(size=2)+
  # xlim(0.88,1.14)+
  scale_color_manual(values=met.brewer("VanGogh1", 3,type ='discrete'))+
  # scale_shape_manual(values = c(19,20,17)) +
  geom_vline(xintercept=1, linetype=3) +
  theme_minimal_hgrid(10, rel_small = 1) +
  facet_wrap(~type, ncol=1)+
  labs(color = "",y = "", x = "Odds Radio",
       title= paste0("Multivariate MR results for ",
                     'Myopia' ,", 95% CI") )+
  theme(legend.position = "none")
dev.off()

