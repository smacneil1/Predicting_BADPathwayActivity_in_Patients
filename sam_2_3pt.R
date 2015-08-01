
library(sva)
library(ASSIGN)
library(utils)
library(devtools)
install_github("wevanjohnson/ASSIGN",ref="adapt_gene_only")
library(ASSIGN)

#Reading in the signature datasets...

source("~/Dropbox/bild_signatures//bild_signatures/Rmarkdowns_scripts//Key_ASSIGN_functions.Rmd")
#setwd("~/Documents/ThesisWork/GitRepos/bild_signature_validation_old_repo/Datasets")
setwd("/Users/mumtahenarahman/Dropbox/bild_signature/Datasets/")
#expr<-as.matrix(read.table("/Users/mumtahenarahman/Dropbox/Datasets/GFP18_AKT_BAD_HER2_IGF1R_RAF_ERK.tpmlog",sep='\t',row.names=1,header=1))

expr<-as.matrix(read.table("~/Documents/ThesisWork/U01Project_copy/GFP18_AKT_BAD_HER2_IGF1R_RAF_ERK.tpmlog",sep='\t',row.names=1,header=1))

control<-subset(expr, select=GFP.1:GFP.12)
her2<-subset(expr, select=HER2.1:HER2.6)
akt<-subset(expr,select=AKT.1:AKT.6)
bad<-subset(expr,select=BAD.1:BAD.6)
igf1r<-subset(expr,select=IGF1R.1:IGF1R.6)
raf<-subset(expr,select=RAF.1:RAF.6)
erk<-subset(expr,select=ERK.1:ERK.6)
expr_all<-cbind(control,akt,bad,her2,igf1r,raf,erk)
dim(expr_all)
expr_all_f <-expr_all[apply(expr_all[,1:47]==0,1,mean) < 0.85,]


pt2<-read.table("~/Documents/ThesisWork/GitRepos/Predicting_BADPathwayActivity_in_Patients/Patient2_RNA_Combined.tpm",header=1,row.names=1)
dim(pt2)
pt3<-read.table("~/Documents/ThesisWork/GitRepos/Predicting_BADPathwayActivity_in_Patients/Patient3_RNA_Combined.tpm",header=1,row.names=1)
head(pt3)
colnames(pt2)<-paste(colnames(pt2),"pt2",sep='_')
colnames(pt3)<-paste(colnames(pt3),"pt3",sep='_')
pt2_3<-merge_drop(pt2,pt3)
head(pt2_3)


# is this just egfr??

expr_all_f_pt_2_3<-merge_drop(expr_all_f,pt2_3,by=0)

#gfp_egfr_pt_2_3<-merge_drop(gfp_egfr_multi_f,pt2_3,by=0)
sub<-c(12,6,6,5,6,6,6,6,11)
#sub<-c(6,6,12,6,6,5,6,6,6,17)
pcaplot(mat = expr_all_f_pt_2_3,sub = sub,scale = T)

#pcaplot(mat = gfp_egfr_pt_2_3,sub = sub,scale = T)
#bat1<-as.matrix(cbind(c(colnames(gfp_egfr_pt_2_3)),as.numeric(c(rep(1,12),rep(2,47),rep(3,6),rep(4,11)))))
bat1<-as.matrix(cbind(c(colnames(expr_all_f_pt_2_3)),as.numeric(c(rep(1,47),rep(2,6),rep(3,11)))))

combat_expr1<-ComBat(dat=expr_all_f_pt_2_3, batch=bat1[,2], mod=NULL)

#combat_expr1<-ComBat(dat=gfp_egfr_pt_2_3, batch=bat1[,2], mod=NULL, numCovs=NULL)
pcaplot(combat_expr1,sub)
c_gfp<-subset(combat_expr1, select=GFP.1:GFP.12)
c_akt<-subset(combat_expr1, select=AKT.1:AKT.6)
c_bad<-subset(combat_expr1, select=BAD.1:BAD.6)
c_her2<-subset(combat_expr1, select=HER2.1:HER2.6)
c_igf1r<-subset(combat_expr1, select=IGF1R.1:IGF1R.6)
c_raf<-subset(combat_expr1, select=RAF.1:RAF.6)
c_erk<-subset(combat_expr1, select=ERK.1:ERK.6)
train_egfr<-combat_expr1[,1:12]
c_test<-combat_expr1[,60:ncol(combat_expr1)]
colnames(c_test)

#######getting the genelist#######
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/akt_75_gene_list/adapB_single/output.rda")
akt_genelist<-output.data$processed.data$diffGeneList  
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/bad_200_gene_list/adapB_single/output.rda")
bad_genelist<-output.data$processed.data$diffGeneList  
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/igf1r_75_gene_list/adapB_single/output.rda")
igf1r_genelist<-output.data$processed.data$diffGeneList  
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/erk_250_gene_list/adapB_single/output.rda")
erk_genelist<-output.data$processed.data$diffGeneList  
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/her2_15_gene_list/adapB_single/output.rda")
her2_genelist<-output.data$processed.data$diffGeneList  
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/egfr_25_gene_list/adapB_single/output.rda")
egfr_genelist<-output.data$processed.data$diffGeneList  
load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/raf_100_gene_list/adapB_single/output.rda")
raf_genelist<-output.data$processed.data$diffGeneList  

basedir="~/Desktop/tmp/pt2/may_12"
dir.create(basedir)

trainingLabel<-list(control=list(akt=1:12),akt=13:18)
sub_dir=paste(basedir,"akt",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_akt),test=c_test,trainingLabel1 = trainingLabel,geneList=c(akt_genelist),out_dir_base = sub_dir)

trainingLabel<-list(control=list(bad=1:12),bad=13:18)
sub_dir=paste(basedir,"bad",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_bad),test=c_test,trainingLabel1 = trainingLabel,geneList=c(bad_genelist),out_dir_base = sub_dir)

trainingLabel<-list(control=list(her2=1:12),her2=13:17)
sub_dir=paste(basedir,"her2",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_her2),test=c_test,trainingLabel1 = trainingLabel,geneList=c(her2_genelist),out_dir_base = sub_dir)

trainingLabel<-list(control=list(igf1r=1:12),igf1r=13:18)
sub_dir=paste(basedir,"igf1r",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_igf1r),test=c_test,trainingLabel1 = trainingLabel,geneList=c(igf1r_genelist),out_dir_base = sub_dir)

trainingLabel<-list(control=list(erk=1:12),erk=13:18)
sub_dir=paste(basedir,"erk",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_erk),test=c_test,trainingLabel1 = trainingLabel,geneList=c(erk_genelist),out_dir_base = sub_dir)

trainingLabel<-list(control=list(raf=1:12),raf=13:18)
sub_dir=paste(basedir,"raf",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_raf),test=c_test,trainingLabel1 = trainingLabel,geneList=c(raf_genelist),out_dir_base = sub_dir)

trainingLabel<-list(control=list(egfr=1:6),egfr=7:12)
sub_dir=paste(basedir,"egfr",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = combat_expr1[,1:12],test=c_test,trainingLabel1 = trainingLabel,geneList=c(egfr_genelist),out_dir_base = sub_dir)
#######

trainingLabela<-list(control=list(akt=1:12),akt=13:18)
sub_dir<-paste(basedir,"akt_75_gene_list",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_akt),test=c_test,trainingLabel1 = trainingLabela,g=75,out_dir_base = sub_dir,single = 1)


trainingLabelb<-list(control=list(bad=1:12),bad=13:18)
sub_dir<-paste(basedir,"bad_200_gene_list",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_bad),test=c_test,trainingLabel1 = trainingLabelb,g=200,out_dir_base = sub_dir,single = 1)

trainingLabelh<-list(control=list(her2=1:12),her2=13:17)
sub_dir<-paste(basedir,"her2_15_gene_list",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_her2),test=c_test,trainingLabel1 = trainingLabelh,g=15,out_dir_base = sub_dir,single = 1)

trainingLabeli<-list(control=list(igf1r=1:12),igf1r=13:18)
sub_dir<-paste(basedir,"igf1r_75_gene_list",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_igf1r),test=c_test,trainingLabel1 = trainingLabeli,g=75,out_dir_base = sub_dir,single = 1)

trainingLabelr<-list(control=list(raf=1:12),raf=13:18)
sub_dir<-paste(basedir,"raf_100_gene_list",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_raf),test=c_test,trainingLabel1 = trainingLabelr,g=100,out_dir_base = sub_dir,single = 1)

trainingLabele<-list(control=list(erk=1:12),erk=13:18)
sub_dir<-paste(basedir,"erk_250_gene_list",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_erk),test=c_test,trainingLabel1 = trainingLabele,g=250,out_dir_base = sub_dir,single = 1)

trainingLabele<-list(control=list(egfr=1:12),egfr=13:18)
sub_dir<-paste(basedir,"egfr_25_gene_list",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = train_egfr,test=c_test,trainingLabel1 = trainingLabele,g=25,out_dir_base = sub_dir,single = 1)
#####
setwd("~/Desktop/tmp/pt2/may_12/")
filenames<-system("ls */adap_adap*/pathway_activity_testset*", intern=TRUE)
filenames
data=NULL
for(i in 1:length(filenames))
{
  f<-read.csv(filenames[i], header=1,row.names=1) ###reading in the filess one at a time
  colnames(f)<-paste(filenames[i],colnames(f),sep='/')
  if(i!=1){
    print(i)
    data<-cbind(data,f)
  }
  else{
    data<-f
  }
}


colnames(data)<-gsub(pattern = "/pathway_activity_testset.csv",replacement = "",x = colnames(data))
write.table(data,"pt2_3_7_pathway_adaptive_with_gene_list.txt",col.names=NA,sep='\t',quote=F)
###SLC6A8 is duplicated in the list of upregulated genes; removed the one with lower coefficient
emt_signature_down<-read.table("~/Dropbox/bild_signatures/bild_signatures/EMT_signature_down.txt",header=T,row.names=1,sep='\t')
emt_signature_up<-read.table("~/Dropbox/bild_signatures/bild_signatures/EMT_signature_up.txt",header=T,row.names=1,sep='\t')
pt2<-read.table("~/Desktop/tmp/pt2/pt2_sam/Patient2_RNA_Combined.tpm",header=1,row.names=1)
dim(pt2)
pt3<-read.table("~/Desktop/tmp/pt2/pt2_sam/Patient3_RNA_Combined.tpm",header=1,row.names=1)
head(pt3)
colnames(pt2)<-paste(colnames(pt2),"pt2",sep='_')
colnames(pt3)<-paste(colnames(pt3),"pt3",sep='_')

pt2_3<-merge_drop(pt2,pt3)
pcaplot(pt2_3,sub,scale=F)
basedir="~/Desktop/tmp/pt2/may_12/"
#processed.data <- assign.preprocess(trainingData=NULL,testData=pt2_3,trainingLabel=NULL,geneList=as.list(c(rownames(emt_signature_down),rownames(emt_signature_up))))
#mcmc.chain <- assign.mcmc(Y=processed.data$testData_sub,Bg = processed.data$B_vector,X=processed.data$S_matrix,Delta_prior_p = processed.data$Pi_matrix,iter = 2000, adaptive_B=TRUE,adaptive_S=FALSE, mixture_beta=TRUE)
#mcmc.pos.mean <- assign.summary(test=mcmc.chain, burn_in=1000,iter=2000, adaptive_B=TRUE, adaptive_S=FALSE,mixture_beta=TRUE)
#assign.output(processed.data=processed.data,mcmc.pos.mean.testData=mcmc.pos.mean,trainingData=NULL, testData=pt2_3,trainingLabel=NULL,testLabel=NULL, geneList=as.list(c(rownames(emt_signature_down),rownames(emt_signature_up))),adaptive_B=FALSE, adaptive_S=FALSE,mixture_beta=TRUE, outputDir=paste(basedir,"emt",sep='/'))



assign.wrapper(trainingData=NULL, testData=pt2_3,trainingLabel=NULL, testLabel=NULL,geneList=as.list(c(rownames(emt_signature_down),rownames(emt_signature_up))), n_sigGene=NULL, adaptive_B=TRUE,adaptive_S=TRUE, mixture_beta=TRUE,outputDir= paste(basedir,"emt",sep='/'), iter=2000, burn_in=1000)
