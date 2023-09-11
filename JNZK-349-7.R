
########################################################
#主要是对单因素的出图的美化，主要是改的话九个药物都需要改

##1、lisinopril

name<-'lisinopril'
path<-paste0('./',name,'/')
dir.create(path)

library(TwoSampleMR)
library(phenoscanner)

exposure_dat <- extract_instruments('ukb-b-12095',p1 = 5e-08)
write.table(exposure_dat, 
            file = paste0(path,"01.exposure_dat_input.xls"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

# res <- phenoscanner(snpquery=c(exposure_dat$SNP))

outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes='ukb-b-964')
write.table(outcome_dat, 
            file = paste0(path,"02.outcome_dat_input.xls"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

dat <- harmonise_data(
  exposure_dat = exposure_dat,
  outcome_dat = outcome_dat)

# run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05)
# 
# dat$EAF2 <- (1 - dat$eaf.exposure)
# dat$MAF <- pmin(dat$eaf.exposure, dat$EAF2)
# PVEfx <- function(BETA, MAF, SE, N){
#   pve <- (2*(BETA^2)*MAF*(1 - MAF))/((2*(BETA^2)*MAF*(1 - MAF)) + ((SE^2)*2*N*MAF*(1 - MAF)))
#   return(pve) 
# }
# dat$PVE <- mapply(PVEfx, dat$beta.exposure, dat$MAF, dat$se.exposure, N = dat$samplesize.exposure)
# dat$FSTAT <- ((dat$samplesize.exposure - 1 - 1)/1)*(dat$PVE/(1 - dat$PVE))
# dat<-dat[dat$FSTAT>10,]

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

p<-mr_leaveoneout_plot(res_loo)[[1]]
p<-p+scale_color_brewer(palette = 'Set2')+
  xlab('MR leave−one−out sensitivity analysis')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'none',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))
  
pdf(file = paste0(path,"07.Leave-one-out.pdf"), height = 8, width = 6)
print(p)
dev.off()

p<-mr_scatter_plot(res, dat)[[1]]
p[["labels"]][["y"]]<-unlist(strsplit(p[["labels"]][["y"]],'[|]'))[1]
p[["labels"]][["x"]]<-unlist(strsplit(p[["labels"]][["x"]],'[|]'))[1]
p[["layers"]][[4]][["data"]]<-p[["layers"]][[4]][["data"]][c(1,3,4),]
p<-p+scale_color_brewer(palette = 'Set1')+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'top',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 12))+ #我不太明白，这个截距到底是值重要还是视觉感受重要？
  theme(aspect.ratio=1)+
  scale_x_continuous(expand = c(0,0),limits = c(0,NA)) 

pdf(file = paste0(path,"08.Scatter.pdf"), height = 8, width = 8)
print(p)
dev.off()

res_single <- mr_singlesnp(dat)
write.table(df_MR_pleio, 
            file = paste0(path,"09.singlesnp_res.xls"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

p<-mr_forest_plot(res_single)[[1]]
p[["data"]]<-p[["data"]][-nrow(p[["data"]]),]  #改完虚线移位了
p[["layers"]][[4]]<-NULL

p<-p+scale_color_brewer(palette = 'Set2')+
  geom_hline(yintercept='', linetype='dotted', col = 'black')+
  xlab('MR effect size')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'none',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))

pdf(file = paste0(path,"10.forest.pdf"), height = 8, width = 7)
print(p)
dev.off()


p<-mr_funnel_plot(res_single)[[1]]
p[["layers"]][[2]][["data"]]<-p[["layers"]][[2]][["data"]][1,]
p<-p+scale_color_brewer(palette = 'Set1')+
  xlab('MR leave−one−out sensitivity analysis')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'top',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))+
  theme(aspect.ratio=1)

pdf(file = paste0(path,"11.funnel.pdf"), height = 8, width = 8)
print(p)
dev.off()

################################################################

##2、bendroflumethiazide

name<-'bendroflumethiazide'
path<-paste0('./',name,'/')
dir.create(path)

library(TwoSampleMR)

exposure_dat <- extract_instruments('ukb-b-18799',p1 = 5e-08)
write.table(exposure_dat, 
            file = paste0(path,"01.exposure_dat_input.xls"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes='ukb-b-964')
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

p<-mr_leaveoneout_plot(res_loo)[[1]]
p<-p+scale_color_brewer(palette = 'Set2')+
  xlab('MR leave−one−out sensitivity analysis')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'none',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))

pdf(file = paste0(path,"07.Leave-one-out.pdf"), height = 8, width = 6)
print(p)
dev.off()

p<-mr_scatter_plot(res, dat)[[1]]
p[["labels"]][["y"]]<-unlist(strsplit(p[["labels"]][["y"]],'[|]'))[1]
p[["labels"]][["x"]]<-unlist(strsplit(p[["labels"]][["x"]],'[|]'))[1]
p[["layers"]][[4]][["data"]]<-p[["layers"]][[4]][["data"]][c(1,3,4),]
p<-p+scale_color_brewer(palette = 'Set1')+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'top',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 12))+ #我不太明白，这个截距到底是值重要还是视觉感受重要？
  theme(aspect.ratio=1)+
  scale_x_continuous(expand = c(0,0),limits = c(0,NA)) 

pdf(file = paste0(path,"08.Scatter.pdf"), height = 8, width = 8)
print(p)
dev.off()

res_single <- mr_singlesnp(dat)
write.table(df_MR_pleio, 
            file = paste0(path,"09.singlesnp_res.xls"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

p<-mr_forest_plot(res_single)[[1]]
p[["data"]]<-p[["data"]][-nrow(p[["data"]]),]  #改完虚线移位了
p[["layers"]][[4]]<-NULL

p<-p+scale_color_brewer(palette = 'Set2')+
  geom_hline(yintercept='', linetype='dotted', col = 'black')+
  xlab('MR effect size')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'none',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))

pdf(file = paste0(path,"10.forest.pdf"), height = 8, width = 7)
print(p)
dev.off()

p<-mr_funnel_plot(res_single)[[1]]
p[["layers"]][[2]][["data"]]<-p[["layers"]][[2]][["data"]][1,]
p<-p+scale_color_brewer(palette = 'Set1')+
  xlab('MR leave−one−out sensitivity analysis')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'top',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))+
  theme(aspect.ratio=1)

pdf(file = paste0(path,"11.funnel.pdf"), height = 8, width = 8)
print(p)
dev.off()

#########################################################

##2、gliclazide

name<-'gliclazide'
path<-paste0('./',name,'/')
dir.create(path)

library(TwoSampleMR)

exposure_dat <- extract_instruments('ukb-b-8602',p1 = 5e-08)
write.table(exposure_dat, 
            file = paste0(path,"01.exposure_dat_input.xls"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes='ukb-b-964')
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

p<-mr_leaveoneout_plot(res_loo)[[1]]
p<-p+scale_color_brewer(palette = 'Set2')+
  xlab('MR leave−one−out sensitivity analysis')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'none',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))

pdf(file = paste0(path,"07.Leave-one-out.pdf"), height = 8, width = 6)
print(p)
dev.off()

p<-mr_scatter_plot(res, dat)[[1]]
p[["labels"]][["y"]]<-unlist(strsplit(p[["labels"]][["y"]],'[|]'))[1]
p[["labels"]][["x"]]<-unlist(strsplit(p[["labels"]][["x"]],'[|]'))[1]
p[["layers"]][[4]][["data"]]<-p[["layers"]][[4]][["data"]][c(1,3,4),]
p<-p+scale_color_brewer(palette = 'Set1')+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'top',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 12))+ #我不太明白，这个截距到底是值重要还是视觉感受重要？
  theme(aspect.ratio=1)+
  scale_x_continuous(expand = c(0,0),limits = c(0,NA)) 

pdf(file = paste0(path,"08.Scatter.pdf"), height = 8, width = 8)
print(p)
dev.off()

res_single <- mr_singlesnp(dat)
write.table(df_MR_pleio, 
            file = paste0(path,"09.singlesnp_res.xls"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

p<-mr_forest_plot(res_single)[[1]]
p[["data"]]<-p[["data"]][-nrow(p[["data"]]),]  #改完虚线移位了
p[["layers"]][[4]]<-NULL

p<-p+scale_color_brewer(palette = 'Set2')+
  geom_hline(yintercept='', linetype='dotted', col = 'black')+
  xlab('MR effect size')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'none',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))

pdf(file = paste0(path,"10.forest.pdf"), height = 8, width = 7)
print(p)
dev.off()

p<-mr_funnel_plot(res_single)[[1]]
p[["layers"]][[2]][["data"]]<-p[["layers"]][[2]][["data"]][1,]
p<-p+scale_color_brewer(palette = 'Set1')+
  xlab('MR leave−one−out sensitivity analysis')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'top',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))+
  theme(aspect.ratio=1)

pdf(file = paste0(path,"11.funnel.pdf"), height = 8, width = 8)
print(p)
dev.off()


########################################################
##3、prednisolone

name<-'prednisolone'
path<-paste0('./',name,'/')
dir.create(path)

library(TwoSampleMR)

exposure_dat <- extract_instruments('ukb-b-8241',p1 = 5e-08)
write.table(exposure_dat, 
            file = paste0(path,"01.exposure_dat_input.xls"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes='ukb-b-964')
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

p<-mr_leaveoneout_plot(res_loo)[[1]]
p<-p+scale_color_brewer(palette = 'Set2')+
  xlab('MR leave−one−out sensitivity analysis')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'none',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))

pdf(file = paste0(path,"07.Leave-one-out.pdf"), height = 8, width = 6)
print(p)
dev.off()

p<-mr_scatter_plot(res, dat)[[1]]
p[["labels"]][["y"]]<-unlist(strsplit(p[["labels"]][["y"]],'[|]'))[1]
p[["labels"]][["x"]]<-unlist(strsplit(p[["labels"]][["x"]],'[|]'))[1]
p[["layers"]][[4]][["data"]]<-p[["layers"]][[4]][["data"]][c(1,3,4),]
p<-p+scale_color_brewer(palette = 'Set1')+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'top',
        legend.text = element_text(size = 12))+ #我不太明白，这个截距到底是值重要还是视觉感受重要？
  theme(aspect.ratio=1)+
  scale_x_continuous(expand = c(0,0),limits = c(0,NA))

pdf(file = paste0(path,"08.Scatter.pdf"), height = 8, width = 8)
print(p)
dev.off()

res_single <- mr_singlesnp(dat)
write.table(df_MR_pleio, 
            file = paste0(path,"09.singlesnp_res.xls"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

p<-mr_forest_plot(res_single)[[1]]
p[["data"]]<-p[["data"]][-nrow(p[["data"]]),]  #改完虚线移位了
p[["layers"]][[4]]<-NULL

p<-p+scale_color_brewer(palette = 'Set2')+
  geom_hline(yintercept='', linetype='dotted', col = 'black')+
  xlab('MR effect size')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'none',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))

pdf(file = paste0(path,"10.forest.pdf"), height = 8, width = 7)
print(p)
dev.off()

p<-mr_funnel_plot(res_single)[[1]]
p[["layers"]][[2]][["data"]]<-p[["layers"]][[2]][["data"]][1,]
p<-p+scale_color_brewer(palette = 'Set1')+
  xlab('MR leave−one−out sensitivity analysis')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'top',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))+
  theme(aspect.ratio=1)

pdf(file = paste0(path,"11.funnel.pdf"), height = 8, width = 8)
print(p)
dev.off()

#################################################################
##4、ramipril

name<-'ramipril'
path<-paste0('./',name,'/')
dir.create(path)

library(TwoSampleMR)

exposure_dat <- extract_instruments('ukb-b-11895',p1 = 5e-08)
write.table(exposure_dat, 
            file = paste0(path,"01.exposure_dat_input.xls"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes='ukb-b-964')
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

p<-mr_leaveoneout_plot(res_loo)[[1]]
p<-p+scale_color_brewer(palette = 'Set2')+
  xlab('MR leave−one−out sensitivity analysis')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'none',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))

pdf(file = paste0(path,"07.Leave-one-out.pdf"), height = 8, width = 6)
print(p)
dev.off()

p<-mr_scatter_plot(res, dat)[[1]]
p[["labels"]][["y"]]<-unlist(strsplit(p[["labels"]][["y"]],'[|]'))[1]
p[["labels"]][["x"]]<-unlist(strsplit(p[["labels"]][["x"]],'[|]'))[1]
p[["layers"]][[4]][["data"]]<-p[["layers"]][[4]][["data"]][c(1,3,4),]
p<-p+scale_color_brewer(palette = 'Set1')+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'top',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 12))+ #我不太明白，这个截距到底是值重要还是视觉感受重要？
  theme(aspect.ratio=1)+
  scale_x_continuous(expand = c(0,0),limits = c(0,NA)) 

pdf(file = paste0(path,"08.Scatter.pdf"), height = 8, width = 8)
print(p)
dev.off()

res_single <- mr_singlesnp(dat)
write.table(df_MR_pleio, 
            file = paste0(path,"09.singlesnp_res.xls"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

p<-mr_forest_plot(res_single)[[1]]
p[["data"]]<-p[["data"]][-nrow(p[["data"]]),]  #改完虚线移位了
p[["layers"]][[4]]<-NULL

p<-p+scale_color_brewer(palette = 'Set2')+
  geom_hline(yintercept='', linetype='dotted', col = 'black')+
  xlab('MR effect size')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'none',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))

pdf(file = paste0(path,"10.forest.pdf"), height = 8, width = 7)
print(p)
dev.off()

p<-mr_funnel_plot(res_single)[[1]]
p[["layers"]][[2]][["data"]]<-p[["layers"]][[2]][["data"]][1,]
p<-p+scale_color_brewer(palette = 'Set1')+
  xlab('MR leave−one−out sensitivity analysis')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'top',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))+
  theme(aspect.ratio=1)

pdf(file = paste0(path,"11.funnel.pdf"), height = 8, width = 8)
print(p)
dev.off()

##############################################################
##4、ramipril

name<-'atenolol'
path<-paste0('./',name,'/')
dir.create(path)

library(TwoSampleMR)

exposure_dat <- extract_instruments('ukb-b-11632',p1 = 5e-08)
write.table(exposure_dat, 
            file = paste0(path,"01.exposure_dat_input.xls"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes='ukb-b-964')
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

p<-mr_leaveoneout_plot(res_loo)[[1]]
p<-p+scale_color_brewer(palette = 'Set2')+
  xlab('MR leave−one−out sensitivity analysis')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'none',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))

pdf(file = paste0(path,"07.Leave-one-out.pdf"), height = 8, width = 6)
print(p)
dev.off()

p<-mr_scatter_plot(res, dat)[[1]]
p[["labels"]][["y"]]<-unlist(strsplit(p[["labels"]][["y"]],'[|]'))[1]
p[["labels"]][["x"]]<-unlist(strsplit(p[["labels"]][["x"]],'[|]'))[1]
p[["layers"]][[4]][["data"]]<-p[["layers"]][[4]][["data"]][c(1,3,4),]
p<-p+scale_color_brewer(palette = 'Set1')+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'top',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 12))+ #我不太明白，这个截距到底是值重要还是视觉感受重要？
  theme(aspect.ratio=1)+
  scale_x_continuous(expand = c(0,0),limits = c(0,NA)) 

pdf(file = paste0(path,"08.Scatter.pdf"), height = 8, width = 8)
print(p)
dev.off()

res_single <- mr_singlesnp(dat)
write.table(df_MR_pleio, 
            file = paste0(path,"09.singlesnp_res.xls"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

p<-mr_forest_plot(res_single)[[1]]
p[["data"]]<-p[["data"]][-nrow(p[["data"]]),]  #改完虚线移位了
p[["layers"]][[4]]<-NULL

p<-p+scale_color_brewer(palette = 'Set2')+
  geom_hline(yintercept='', linetype='dotted', col = 'black')+
  xlab('MR effect size')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'none',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))

pdf(file = paste0(path,"10.forest.pdf"), height = 8, width = 7)
print(p)
dev.off()

p<-mr_funnel_plot(res_single)[[1]]
p[["layers"]][[2]][["data"]]<-p[["layers"]][[2]][["data"]][1,]
p<-p+scale_color_brewer(palette = 'Set1')+
  xlab('MR leave−one−out sensitivity analysis')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'top',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))+
  theme(aspect.ratio=1)

pdf(file = paste0(path,"11.funnel.pdf"), height = 8, width = 8)
print(p)
dev.off()

##############################################################
##4、paracetamol

name<-'paracetamol'
path<-paste0('./',name,'/')
dir.create(path)

library(TwoSampleMR)

exposure_dat <- extract_instruments('ukb-a-194',p1 = 5e-08)
write.table(exposure_dat, 
            file = paste0(path,"01.exposure_dat_input.xls"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes='ukb-b-964')
write.table(outcome_dat, 
            file = paste0(path,"02.outcome_dat_input.xls"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

dat <- harmonise_data(
  exposure_dat = exposure_dat,
  outcome_dat = outcome_dat)

dat<-dat[-1,]

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

p<-mr_leaveoneout_plot(res_loo)[[1]]
p<-p+scale_color_brewer(palette = 'Set2')+
  xlab('MR leave−one−out sensitivity analysis')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'none',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))

pdf(file = paste0(path,"07.Leave-one-out.pdf"), height = 8, width = 6)
print(p)
dev.off()

p<-mr_scatter_plot(res, dat)[[1]]
p[["labels"]][["y"]]<-unlist(strsplit(p[["labels"]][["y"]],'[|]'))[1]
p[["labels"]][["x"]]<-unlist(strsplit(p[["labels"]][["x"]],'[|]'))[1]
p[["layers"]][[4]][["data"]]<-p[["layers"]][[4]][["data"]][c(1,3,4),]
p<-p+scale_color_brewer(palette = 'Set1')+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'top',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 12))+ #我不太明白，这个截距到底是值重要还是视觉感受重要？
  theme(aspect.ratio=1)+
  scale_x_continuous(expand = c(0,0),limits = c(0,NA)) 

pdf(file = paste0(path,"08.Scatter.pdf"), height = 8, width = 8)
print(p)
dev.off()

res_single <- mr_singlesnp(dat)
write.table(df_MR_pleio, 
            file = paste0(path,"09.singlesnp_res.xls"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

p<-mr_forest_plot(res_single)[[1]]
p[["data"]]<-p[["data"]][-nrow(p[["data"]]),]  #改完虚线移位了
p[["layers"]][[4]]<-NULL

p<-p+scale_color_brewer(palette = 'Set2')+
  geom_hline(yintercept='', linetype='dotted', col = 'black')+
  xlab('MR effect size')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'none',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))

pdf(file = paste0(path,"10.forest.pdf"), height = 8, width = 7)
print(p)
dev.off()

p<-mr_funnel_plot(res_single)[[1]]
p[["layers"]][[2]][["data"]]<-p[["layers"]][[2]][["data"]][1,]
p<-p+scale_color_brewer(palette = 'Set1')+
  xlab('MR leave−one−out sensitivity analysis')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'top',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))+
  theme(aspect.ratio=1)

pdf(file = paste0(path,"11.funnel.pdf"), height = 8, width = 8)
print(p)
dev.off()


######################################
##4、atorvastatin

name<-'atorvastatin'
path<-paste0('./',name,'/')
dir.create(path)

library(TwoSampleMR)

exposure_dat <- extract_instruments('ukb-b-10008',p1 = 5e-08)
write.table(exposure_dat, 
            file = paste0(path,"01.exposure_dat_input.xls"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes='ukb-b-964')
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

p<-mr_leaveoneout_plot(res_loo)[[1]]
p<-p+scale_color_brewer(palette = 'Set2')+
  xlab('MR leave−one−out sensitivity analysis')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'none',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))

pdf(file = paste0(path,"07.Leave-one-out.pdf"), height = 8, width = 6)
print(p)
dev.off()

p<-mr_scatter_plot(res, dat)[[1]]
p[["labels"]][["y"]]<-unlist(strsplit(p[["labels"]][["y"]],'[|]'))[1]
p[["labels"]][["x"]]<-unlist(strsplit(p[["labels"]][["x"]],'[|]'))[1]
p[["layers"]][[4]][["data"]]<-p[["layers"]][[4]][["data"]][c(1,3,4),]
p<-p+scale_color_brewer(palette = 'Set1')+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'top',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 12))+ #我不太明白，这个截距到底是值重要还是视觉感受重要？
  theme(aspect.ratio=1)+
  scale_x_continuous(expand = c(0,0),limits = c(0,NA)) 

pdf(file = paste0(path,"08.Scatter.pdf"), height = 8, width = 8)
print(p)
dev.off()

res_single <- mr_singlesnp(dat)
write.table(df_MR_pleio, 
            file = paste0(path,"09.singlesnp_res.xls"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

p<-mr_forest_plot(res_single)[[1]]
p[["data"]]<-p[["data"]][-nrow(p[["data"]]),]  #改完虚线移位了
p[["layers"]][[4]]<-NULL

p<-p+scale_color_brewer(palette = 'Set2')+
  geom_hline(yintercept='', linetype='dotted', col = 'black')+
  xlab('MR effect size')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'none',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))

pdf(file = paste0(path,"10.forest.pdf"), height = 8, width = 7)
print(p)
dev.off()

p<-mr_funnel_plot(res_single)[[1]]
p[["layers"]][[2]][["data"]]<-p[["layers"]][[2]][["data"]][1,]
p<-p+scale_color_brewer(palette = 'Set1')+
  xlab('MR leave−one−out sensitivity analysis')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'top',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))+
  theme(aspect.ratio=1)

pdf(file = paste0(path,"11.funnel.pdf"), height = 8, width = 8)
print(p)
dev.off()


##########################################
##4、atorvastatin

name<-'amlodipine'
path<-paste0('./',name,'/')
dir.create(path)

library(TwoSampleMR)

exposure_dat <- extract_instruments('ukb-b-9207',p1 = 5e-08)
write.table(exposure_dat, 
            file = paste0(path,"01.exposure_dat_input.xls"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes='ukb-b-964')
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

p<-mr_leaveoneout_plot(res_loo)[[1]]
p<-p+scale_color_brewer(palette = 'Set2')+
  xlab('MR leave−one−out sensitivity analysis')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'none',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))

pdf(file = paste0(path,"07.Leave-one-out.pdf"), height = 8, width = 6)
print(p)
dev.off()

p<-mr_scatter_plot(res, dat)[[1]]
p[["labels"]][["y"]]<-unlist(strsplit(p[["labels"]][["y"]],'[|]'))[1]
p[["labels"]][["x"]]<-unlist(strsplit(p[["labels"]][["x"]],'[|]'))[1]
p[["layers"]][[4]][["data"]]<-p[["layers"]][[4]][["data"]][c(1,3,4),]
p<-p+scale_color_brewer(palette = 'Set1')+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'top',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 12))+ #我不太明白，这个截距到底是值重要还是视觉感受重要？
  theme(aspect.ratio=1)+
  scale_x_continuous(expand = c(0,0),limits = c(0,NA)) 

pdf(file = paste0(path,"08.Scatter.pdf"), height = 8, width = 8)
print(p)
dev.off()

res_single <- mr_singlesnp(dat)
write.table(df_MR_pleio, 
            file = paste0(path,"09.singlesnp_res.xls"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

p<-mr_forest_plot(res_single)[[1]]
p[["data"]]<-p[["data"]][-nrow(p[["data"]]),]  #改完虚线移位了
p[["layers"]][[4]]<-NULL

p<-p+scale_color_brewer(palette = 'Set2')+
  geom_hline(yintercept='', linetype='dotted', col = 'black')+
  xlab('MR effect size')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'none',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))

pdf(file = paste0(path,"10.forest.pdf"), height = 8, width = 7)
print(p)
dev.off()

p<-mr_funnel_plot(res_single)[[1]]
p[["layers"]][[2]][["data"]]<-p[["layers"]][[2]][["data"]][1,]
p<-p+scale_color_brewer(palette = 'Set1')+
  xlab('MR leave−one−out sensitivity analysis')+
  theme_classic(base_size = 16)+
  theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
        axis.ticks = element_line(size = 1),
        legend.background = element_blank(),
        legend.position = 'top',
        # legend.title = element_blank(),
        # legend.position = c(0.1,0.9),
        legend.text = element_text(size = 18))+
  theme(aspect.ratio=1)

pdf(file = paste0(path,"11.funnel.pdf"), height = 8, width = 8)
print(p)
dev.off()

