##################################################
#这比自己写的循环快多了

library(TwoSampleMR)
library(doParallel)

rm(list = ls())
exp_data <- readRDS("./exp_data_2023.4.3_p5e8_idALL.rds")
load('./available_outcomes.RData')
length(unique(exp_data$id.exposure))
length(intersect(ao$id,exp_data$id.exposure))
# 总共26279个暴露信息

# 去掉eqtl的
delete_eqtl <- stringr::str_detect(
  string = exp_data$id.exposure,
  pattern = "eqtl-a", negate = TRUE )
exp_data <- subset(exp_data,delete_eqtl)
length(unique(exp_data$id.exposure))
# 有10090个暴露信息

# 筛选人群, 你研究的人群确定下
table(ao$population) 
# id2 <- ao

id2 <- subset(ao,ao$population=="European")
length(unique(id2$id))

# # 测试选取10000个id部分
# 38065
# id2 <- id2[35001:38065,]  # 实际运行注解掉

dplyr::count(id2,population)
nrow(exp_data)

exp_data <- subset(exp_data,
                   exp_data$id.exposure %in% id2$id)
nrow(exp_data)

start1=Sys.time()

# IEU在线IEU的
# data(gwas_catalog)
outcome=extract_outcome_data(
  snps = unique(exp_data$SNP),
  outcomes = c("ukb-b-6353"), # 近视
  proxies = TRUE)
# save(outcome,file = 'outcomefinn-b-H7_MYOPIA.rda')
end1=Sys.time();end1-start1

exp_data_list=split(exp_data,list(exp_data$id.exposure))

har_loop <- function(exp_data=exp_data_list){
  BBB=TwoSampleMR::harmonise_data(
    exposure_dat = exp_data, outcome_dat = outcome)
  return(BBB)
}

start1=Sys.time()
dat_list <- list()
for(x in 1:length(exp_data_list)){
  dat_list[[names(exp_data_list)[x]]] <- har_loop(exp_data=exp_data_list[[x]])
}
end1=Sys.time();end1-start1

dat <- do.call(rbind,dat_list)
dat <- subset(dat,mr_keep)

dat <- split(dat,list(dat$id.exposure,dat$id.outcome))
length(dat)
names(dat) <- paste0("A",1:length(dat))
deleteSNP=names(dat)[sapply(dat,nrow)<3]
length(deleteSNP)

# [1]-去掉list里面小于3个SNP的
for (deleteSNPid in deleteSNP) {
  dat[[deleteSNPid]] <- NULL
}
length(dat)
dat2 <- do.call(rbind,dat)
length(unique(dat2$id.exposure)) 

# 并行
library(doParallel)
#
choose_MR <- function(dat1=dat){ 
  res_hete <- mr_heterogeneity(dat1) 
  if(nrow(res_hete)==1 & !grepl('Invers',res_hete$method[1])){next}
  if(nrow(res_hete)==0 ){next}
  if (res_hete$Q_pval[nrow(res_hete)]<0.05) {
    res=TwoSampleMR::mr(dat1, method_list = c("mr_egger_regression",
                                              "mr_weighted_median", "mr_ivw_mre"))
  } else{
    res=TwoSampleMR::mr(dat1, method_list = c("mr_egger_regression",
                                              "mr_weighted_median", "mr_ivw_fe"))
  }
  return(res)
}

start1=Sys.time()
res_list <- list()
for(x in 1:length(dat)){
  res_list[[names(dat)[x]]] <- choose_MR(dat1=dat[[x]])
}
end1=Sys.time();end1-start1
res <- do.call(rbind,res_list)

write.csv(res,'res_6353.csv')

# 结果提取
res$pval=round(res$pval,3)

res_ALL <- split(res,list(res$id.exposure))
#
judge_1 <- function(mr_res=res2) {
  mr_res$b_direction <- as.numeric(sign(mr_res$b))
  mr_res$b_direction=ifelse(abs(sum(mr_res$b_direction))==3 ,
                            NA,"Inconsistent direction")
  mr_res$p_no <- NA
  mr_res[mr_res$method=="MR Egger","p_no"] <- ifelse(
    mr_res[mr_res$method=="MR Egger","pval"]<0.05," ",
    "MR Egger")
  mr_res[mr_res$method=="Weighted median","p_no"] <- ifelse(
    mr_res[mr_res$method=="Weighted median","pval"]<0.05," ",
    "Weighted median")
  mr_res[grep(x = mr_res$method,pattern = "Inverse variance"),"p_no"] <- ifelse(
    mr_res[grep(x = mr_res$method,pattern = "Inverse variance"),"pval"]<0.05,
    " ","IVW")
  mr_res$p_no <- paste(mr_res$p_no,collapse = " ")
  mr_res$p_no=trimws(mr_res$p_no,which = c("both"))
  return(mr_res)
}

res_ALL=purrr::map(.x =res_ALL,.f = ~judge_1(.x) )
res_ALL2 <- do.call(rbind,res_ALL)
res_ALL3 <- subset(res_ALL2,
                   is.na(res_ALL2$b_direction) )
bool=stringr::str_detect(string =res_ALL3$p_no,
                         pattern = "IVW",negate = TRUE )
res_ALL4 <- subset(res_ALL3,bool)

res_ALL4_1=subset(res_ALL4,select = exposure)
res_ALL4_1 <- unique(res_ALL4_1)
res_ALL4_1[1,1]
res_ALL4_2 <- tidyr::separate(
  data = res_ALL4_1,col = exposure,sep = "\\|",
  into = c("exposure","delete")) %>%
  dplyr::select(-delete)

# 导出
library(openxlsx)

wb <- createWorkbook("My name here")
## Add a worksheets
addWorksheet(wb, "sheet1", gridLines = FALSE)
addWorksheet(wb, "sheet2", gridLines = FALSE)
## write data to worksheet 1
writeData(wb,x = res_ALL4,sheet = "sheet1",
          rowNames = FALSE)
writeData(wb,x = res_ALL4_2,sheet = "sheet2",
          rowNames = FALSE)
## style for body
bodyStyle <- createStyle(border = "TopBottom",
                         bgFill ="#e3e9f4",  
                         fgFill = "#e3e9f4")
a=seq(2,nrow(res_ALL4)+1,6)
b=seq(3,nrow(res_ALL4)+1,6)
c=seq(4,nrow(res_ALL4)+1,6)
d=sort(c(a,b,c))
d
addStyle(wb, sheet = 1, bodyStyle, 
         rows = d,
         cols = 1:11, 
         gridExpand = TRUE)
setColWidths(wb, 1, cols = 1, widths = 21) ## set column width for row names column
## Not run: 
saveWorkbook(wb, "MR_6353.xlsx", 
             overwrite = TRUE)


