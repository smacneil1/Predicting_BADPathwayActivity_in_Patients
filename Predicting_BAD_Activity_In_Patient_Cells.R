# Predicting Bad Pathway Activity in Patient Samples for Sam Brady's Tumor Het. Project

source("~/Documents/ThesisWork/GitRepos/Predicting_BADPathwayActivity_in_Patients/Key_ASSIGN_functions.Rmd")
library(ASSIGN)
library(sva)

# get this file

expr<-as.matrix(read.table("~/Documents/ThesisWork/U01Project_copy/GFP18_AKT_BAD_HER2_IGF1R.tpmlog",sep='\t',row.names=1,header=1))
dim(expr)
expr_all_f_allsigs <-expr[apply(expr[,1:18]==0,1,mean) < 0.85,]
dim(expr_all_f_allsigs)

BAD<-subset(expr,select=BAD.1:BAD.6)
dim(BAD)
control<-subset(expr, select=GFP.1:GFP.12)
expr_all<-cbind(control,BAD)
dim(expr_all)
expr_all_f <-expr_all[apply(expr_all[,1:18]==0,1,mean) < 0.85,]
dim(expr_all_f)

pt2<-read.table("~/Documents/ThesisWork/GitRepos/Predicting_BADPathwayActivity_in_Patients/Patient2_RNA_Combined.tpm",header=1,row.names=1)
dim(pt2)

pt2_log=log2(pt2+1)
head(pt2_log)
pt3<-read.table("~/Documents/ThesisWork/GitRepos/Predicting_BADPathwayActivity_in_Patients/Patient3_RNA_Combined.tpm",header=1,row.names=1)
head(pt3)
colnames(pt2)<-paste(colnames(pt2),"pt2",sep='_')
colnames(pt3)<-paste(colnames(pt3),"pt3",sep='_')
head(pt2_3)

#gfp_bad_pt_2_log<-merge_drop(expr_all_f,pt2_log,by=0)

gfp_bad_pt_2<-merge_drop(expr_all_f,pt2,by=0)
head(gfp_bad_pt_2)
gfp_bad_pt_3<-merge_drop(expr_all_f,pt3,by=0)
gfp_bad_pt_2_3<-merge_drop(expr_all_f,pt2_3,by=0)
View(gfp_bad_pt_2_3)

#try with all the signatures
allsigs_pt_2<-merge_drop(expr_all_f_allsigs,pt2,by=0)
allsigs_pt_3<-merge_drop(expr_all_f_allsigs,pt3,by=0)
allsigs_pt_2_3<-merge_drop(expr_all_f_allsigs,pt2_3,by=0)
View(allsigs_pt_2)
sub_2<-c(12,6,6)
sub_3<-c(12,6,11)
sub_2_3<-c(12,6,6,11)
sub_all_2<-c(12,6,6,5,6,6)
sub_all_3<-c(12,6,6,5,6,11)
sub_all_2_3<-c(12,6,6,5,6,6,11)

bat_2_3<-as.matrix(cbind(c(colnames(gfp_bad_pt_2_3)),as.numeric(c(rep(1,ncol(expr_all_f)),rep(2,ncol(pt2)),rep(3,ncol(pt3))))))
bat_3<-as.matrix(cbind(c(colnames(gfp_bad_pt_3)),as.numeric(c(rep(1,ncol(expr_all_f)),rep(2,ncol(pt3))))))
bat_2<-as.matrix(cbind(c(colnames(gfp_bad_pt_2)),as.numeric(c(rep(1,ncol(expr_all_f)),rep(2,ncol(pt2))))))
bat_all_2_3<-as.matrix(cbind(c(colnames(allsigs_pt_2_3)),as.numeric(c(rep(1,ncol(expr_all_f_allsigs)),rep(2,ncol(pt2)), rep(3,ncol(pt3))))))
bat_all_2<-as.matrix(cbind(c(colnames(allsigs_pt_2)),as.numeric(c(rep(1,ncol(expr_all_f_allsigs)),rep(2,ncol(pt2))))))
bat_all_3<-as.matrix(cbind(c(colnames(allsigs_pt_3)),as.numeric(c(rep(1,ncol(expr_all_f_allsigs)),rep(2,ncol(pt3))))))

pdf("~/Desktop/Patient2_3_PCA_Plots.pdf")

pcaplot(mat = gfp_bad_pt_2,sub = sub_2,scale = T)
pcaplot(combat_expr_2,sub_2)
pcaplot(mat = gfp_bad_pt_3,sub = sub_3,scale = T)
pcaplot(combat_expr_3,sub_3)
pcaplot(mat = gfp_bad_pt_2_3,sub = sub_2_3,scale = T)
pcaplot(combat_expr_2_3,sub_2_3)

dev.off()

pcaplot(mat = allsigs_pt_2,sub = sub_all_2,scale = T)
pcaplot(mat = allsigs_pt_3,sub = sub_all_3,scale = T)
pcaplot(mat = allsigs_pt_2_3,sub = sub_all_2_3,scale = T)


combat_expr_2<-ComBat(dat=gfp_bad_pt_2, batch=bat_2[,2], mod=NULL)

combat_expr_3<-ComBat(dat=gfp_bad_pt_3, batch=bat_3[,2], mod=NULL)

combat_expr_2_3<-ComBat(dat=gfp_bad_pt_2_3, batch=bat_2_3[,2], mod=NULL)



combat_expr_all_2<-ComBat(dat=allsigs_pt_2_3, batch=bat_all_2[,2], mod=NULL)
combat_expr_all_3<-ComBat(dat=allsigs_pt_2, batch=bat_all_3[,2], mod=NULL)
combat_expr_all_2_3<-ComBat(dat=allsigs_pt_3, batch=bat_all_2_3[,2], mod=NULL)














pcaplot(combat_expr1,sub)


c_bad<-subset(combat_expr1, select=BAD.1:BAD.6)

train_bad<-combat_expr1[,1:12]
c_test<-combat_expr1[,60:ncol(combat_expr1)]
colnames(c_test)

#######getting the genelist#######

load("~/Dropbox/bild_signatures/icbp_15_april_assign_adap_/bad_200_gene_list/adapB_single/output.rda")
bad_genelist<-output.data$processed.data$diffGeneList  

basedir="~/Desktop/tmp/pt2/may_12"
dir.create(basedir)


trainingLabel<-list(control=list(bad=1:12),bad=13:18)
sub_dir=paste(basedir,"bad",sep='/')
dir.create(sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_bad),test=c_test,trainingLabel1 = trainingLabel,geneList=c(bad_genelist),out_dir_base = sub_dir)


trainingLabelb<-list(control=list(bad=1:12),bad=13:18)
sub_dir<-paste(basedir,"bad_200_gene_list",sep='/')
dir.create( sub_dir)
assign_easy_multi(trainingData = cbind(c_gfp,c_bad),test=c_test,trainingLabel1 = trainingLabelb,g=200,out_dir_base = sub_dir,single = 1)


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


