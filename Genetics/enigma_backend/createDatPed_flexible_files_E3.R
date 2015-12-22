#
# * createDatPed.R  create .dat and .ped files for running GWAS on ENIGMA ROIs
# * Derrek Hibar - derrek.hibar@ini.usc.edu
# * Neda Jahanshad - neda.jahanshad@ini.usc.edu
# * ENIGMA 2014.
#

###### Basic INFO

###################################################################
## to run:
# ${Rbin} --no-save --slave --args ${1} ${2} ...   ${8} <  ./.R
# R --no-save --slave --args ${csvFILE} ${localfamFILE} ${combinedROItableFILE} ${ageColumnHeader} ${sexColumnHeader} ${maleIndicator} ${patients} ${related} ${peddatdir} ${ALL_ROIS} ${eName} <  ${run_directory}/createDatPed_flexible_files.R
###################################################################

# 8 INPUTS
######  path to your HM3mds2Rmds.csv file
######  path to your combinedROItable.csv (output from ROI script) file
######  the column header for your age covariate
######  the column header for your sex covariate
######  what is the indicator for males in the sex column (M? 1? 2? ... )
######  does your dataset contain patients (0 for no, 1 for yes)
######  does your dataset contain related individuals (0 for no, 1 for yes)
######  output folder to write the ped and dat files into (doesn't have to exist yet!)
####################################################################

options(stringsAsFactors = FALSE)

versionDate=20151222

# Output files will include .dat and .ped files with all necessary covariates
# standard covariates include age, sex, age^2, age-x-sex, age^2-x-sex, and 4 MDS components
###### sex and interaction terms not included if single-sex study
###### AffectionStatus included if patients are included
###### Any site specific covariates, including imaging site (if multiple) should be included.
###### Each site needs its own dummy covariate column
###### Each patient group needs its own dummy covariate column
####################################################################
cmdargs = commandArgs(trailingOnly=T);
HM3mds2RmdscsvFILE=cmdargs[1]
HM3mergedfamFILE=cmdargs[2]
combinedROItableFILE=cmdargs[3]
ageColumnHeader=cmdargs[4]
sexColumnHeader=cmdargs[5]

maleIndicator=cmdargs[6]
if (is.na(as.numeric(maleIndicator)) == "TRUE"){
    maleIndicator=maleIndicator
    } else {
    maleIndicator=as.numeric(maleIndicator)
}

patients=cmdargs[7] ## have a column where all healthy are marked as 0s and all patients as 1
if (is.na(as.numeric(patients)) == "TRUE"){
    patients=patients
} else {
    patients=as.numeric(patients)
}


related=cmdargs[8]  ## does your dataset contain related individuals (0 for no, 1 for yes)
if (is.na(as.numeric(related)) == "TRUE"){
    related=related
} else {
    related=as.numeric(related)
}
outFolder=cmdargs[9]

####################################################################
ALL_ROIS=cmdargs[10]
ALL_ROIS=as.character(parse(text=ALL_ROIS))

eName=cmdargs[11]

dir.create(outFolder, showWarnings = FALSE)
zz <- file(paste(outFolder,"/","RUN_NOTES.txt",sep=""),"w")
####################################################################
source('ENIGMA_functions.R')

####################################################################
#ALL_ROIS=c("ACR","ALIC","AverageFA","BCC","CC","CGC","CGH","CR","CST","EC","FX","FXST","GCC","IC","IFO","PCR","PLIC","PTR","RLIC","SCC","SCR","SFO","SLF","SS","UNC");

ALL_IDS=c("FID","IID","PID","MID","Sex");

Nids=length(ALL_IDS)
Nrois=length(ALL_ROIS)
####################################################################

InfoFile <- data.frame(read.csv(combinedROItableFILE,colClasses = "character")) #Read in the phenotypes + covariates file

InfoFile$IID = InfoFile[,1] #This just renames a column for easier merging
InfoFile$subjectID = NULL
InfoFile$SubjID = NULL
mds.cluster <- read.csv(HM3mds2RmdscsvFILE,colClasses = "character")  #Read in the MDS components
mds.cluster$SOL <- NULL; #Remove the “SOL” column in the MDS components since this is not a covariate to be included

if(related == 1){
    fam <-read.table(HM3mergedfamFILE,colClasses = "character")
	cat('You have indicated that your sample contains RELATED subjects\n')
} else {
	cat('You have indicated that your sample contains only UNRELATED subjects\n')
}

missing="";
l_missing=0;
for (i in 1:Nrois) {
    columnnames = colnames(InfoFile);
    if (length(InfoFile[,which(columnnames==ALL_ROIS[i])]) == 0) {
        missing=paste(missing,ALL_ROIS[i])
        l_missing=l_missing+1;
    }
}

if (l_missing > 0) {
    stop(paste("ERROR: You are missing the following ROIs:", missing ,". Please re-run latest ROI script!",sep=""))
}

cat('Merging your MDS files with your Phenotype and Covariates Files...')
merged_temp <- merge(mds.cluster, InfoFile, by="IID"); #Merge the MDS and other covariates
cat('Done\n')

### make sure subject names match up!
numsubjects=length(merged_temp[,1])
if (numsubjects==0){
    stop("ERROR: Please make sure your subjectID's in your phenotype csv file are the same as those listed in your HM3mds2R.mds.csv file")
}

### make sure there are no duplicates
dups<-duplicated(merged_temp$IID)
Ldups<-length(which(dups=="TRUE"))  ##**##** depending on R version this (& below) may have to be which(dups,"TRUE")
cat('    There are ',Ldups,' duplicate subjects.\n')
if (Ldups > 0){
    merged_temp<-merged_temp[-which(dups=="TRUE"),]
    print('    Duplicates have been removed.\n')
}

numsubjects = length(merged_temp$IID);

writeLines(paste0('VERSION: ',substr(versionDate,5,6),'/',substr(versionDate,7,8),'/',substr(versionDate,0,4)),con=zz,sep="\n")

### if no related individuals, create dummy paternal and maternal IDs
# otherwise break to make sure these are entered
if (related==0) {
        writeLines(paste('STUDY DESIGN: There are no related individuals.'),con=zz,sep="\n")
        if ( (length(merged_temp[,which(columnnames=="MID")])==0 ) || (length(merged_temp[,which(columnnames=="PID")])==0) ) {
        p=matrix(0,nrow=numsubjects,ncol=1)
        #FID=merged_temp$IID
        PID=data.frame("PID"=p)
        MID=data.frame("MID"=p) 
        zygosity=data.frame("zygosity"=p)}
    }  else if ((related==1) &&  (length(merged_temp[,which(columnnames=="MID")])==0 ) && (length(merged_temp[,which(columnnames=="PID")])==0) ) {
          writeLines(paste('STUDY DESIGN: There are related individuals. We will add in the MID and PID from the local.fam file because you did not put it in your csv file...  '),con=zz,sep="\n")
          matchind = match(merged_temp$IID, fam[,2]);
          merged_temp$MID=fam[matchind,4] ####
          merged_temp$PID=fam[matchind,3] ####
          p=matrix(0,nrow=numsubjects,ncol=1)
          
          #### added July 2015 -- in case no MID or PID are in the plink files, ERROR to create new parental IDs
	        if ( (merged_temp$MID == merged_temp$PID) || ( sum(merged_temp$MID) == 0 ) || ( sum(merged_temp$PID) == 0 ) ) {
	  	    stop("ERROR: You have indicated a related cohort, yet no unique Maternal or Paternal IDs were detected in your data file, or in your  HM3mds2R.mds.csv file, so we can't compile it for you! Please go back to your data csv file and add in two columns MID and PID and enter a unique ID for each parent of the family. You do not have to have phenotype data on them, but this is necessary to correctly model the known family relatedness. ")
	        }
          
    }  else if ((related==1) && (length(merged_temp[,which(columnnames=="zygosity")])==0 )){
          writeLines(paste('STUDY DESIGN: There are related individuals.'),con=zz,sep="\n")
          p=matrix(0,nrow=numsubjects,ncol=1)
          zygosity=data.frame("zygosity"=p)
          print("WARNING: It appears this study is of related individuals but no zygosity column was detected. We are assuming there are no twins in your sample. If this is not the case, please check your combinedROItable.csv file and add in a zygosity column to indicate MZ pairs with odd numbers and DZ pairs with even numbers. An estimated zygosity column will be created from PID, MID, sex, and MDS information. ")
          writeLines(paste('    WARNING: It appears this study is of related individuals but no zygosity column was detected.  If this is not the case, please check your combinedROItable.csv file and add in a zygosity column to indicate MZ pairs with odd numbers and DZ pairs with even numbers. An estimated zygosity column will be created from PID, MID, age, sex, and MDS information.'),con=zz,sep="\n")
    }  else {
        writeLines(paste('STUDY DESIGN: This is a family based study and there are related individuals. PID, MID, and zygosity columns have been identified.'),con=zz,sep="\n")
}

### Find age and sex columns, center, and create new age^2, age-x-sex and age^2-x-sex columns
columnnames = colnames(merged_temp);
age=as.numeric(merged_temp[,which(columnnames==ageColumnHeader)])
age_mean=mean(age)
age_md=median(age)
min_age=min(age)
max_age=max(age)

ageC=(age-age_mean)
ageCsq=ageC*ageC

sex=merged_temp[,which(columnnames==sexColumnHeader)]
males=which(sex==maleIndicator)
sexC=matrix(0,nrow=numsubjects,ncol=1)
sexC[males]<- -0.5
sexC[-males]<- 0.5

StandardSex=data.frame("Sex"=sexC)
StandardSex[which(sexC==-.5),]<- 1
StandardSex[which(sexC==.5),]<- 2

Nm=length(which(sexC==-.5))
Nf=length(which(sexC==.5))

writeLines(paste('Distribution stats for age in the full group:'),con=zz,sep="\n")
writeLines(paste('		mean:',age_mean),con=zz,sep="\n")
writeLines(paste('		median:',age_md),con=zz,sep="\n")
writeLines(paste('		min:',min_age),con=zz,sep="\n")
writeLines(paste('		max:',max_age),con=zz,sep="\n")
writeLines(paste('There are',Nm,'males in this study.'),con=zz,sep="\n")
writeLines(paste('There are',Nf,'females in this study.'),con=zz,sep="\n")

age_sexC=age*sexC
ageCsq_sexC=ageCsq*sexC

### Do not include sex or sex interaction variables if population is all M or all F or age stuff if all the same age
if (sd(sexC) ==0) {
        print("WARNING: It appears this study is of a single sex. If this is not the case, please check your combinedROItable.csv file")
        writeLines(paste('  WARNING: It appears this study is of a single sex. If this is not the case, please check your combinedROItable.csv file.'),con=zz,sep="\n")
        age_sexC=age_sexC*0
        ageCsq_sexC=ageCsq_sexC*0
}

if (sd(ageC) ==0) {
        print("WARNING: It appears this study is of a single age group. If this is not the case, please check your combinedROItable.csv file")
        writeLines(paste('  WARNING: It appears this study is of a single age group. If this is not the case, please check your combinedROItable.csv file.'),con=zz,sep="\n")
        ageCsq=ageCsq*0
        age_sexC=age_sexC*0
        ageCsq_sexC=ageCsq_sexC*0
}

## set columns as variables
attach(merged_temp)

merged_temp_rest=as.data.frame(merged_temp)
columnnames = colnames(merged_temp_rest);
merged_temp_rest=merged_temp_rest[,-which(columnnames==ageColumnHeader)]
columnnames = colnames(merged_temp_rest);
merged_temp_rest=merged_temp_rest[,-which(columnnames==sexColumnHeader)]
columnnames = colnames(merged_temp_rest);
merged_temp_rest=merged_temp_rest[,-which(columnnames=="FID")]
columnnames = colnames(merged_temp_rest);
merged_temp_rest=merged_temp_rest[,-which(columnnames=="IID")]
columnnames = colnames(merged_temp_rest);
if (length(merged_temp_rest[,which(columnnames=="MID")]) > 0) {
    merged_temp_rest=merged_temp_rest[,-which(columnnames=="MID")]
}
columnnames = colnames(merged_temp_rest);
if (length(merged_temp_rest[,which(columnnames=="PID")]) > 0) {
    merged_temp_rest=merged_temp_rest[,-which(columnnames=="PID")]
}
columnnames = colnames(merged_temp_rest);
if (length(merged_temp_rest[,which(columnnames=="zygosity")]) > 0) {
    merged_temp_rest=merged_temp_rest[,-which(columnnames=="zygosity")]
}
columnnames = colnames(merged_temp_rest);

VarNames=names(merged_temp_rest)

FullInfoFile=cbind(FID,IID,PID,MID,StandardSex,sexC,age,ageCsq,age_sexC,ageCsq_sexC,merged_temp_rest)

VarNames=names(FullInfoFile)
numsubjects = length(FullInfoFile$IID);

nVar=dim(FullInfoFile)[2]
Nset=Nids+Nrois
nCov=nVar-Nset ### all the covariates
FullInfoFile=moveMe(FullInfoFile,ALL_ROIS,"after","Sex")

numsubjects = length(FullInfoFile$IID);

if (sum(is.na(FullInfoFile$FID)) > 0 ) {
	FullInfoFile=FullInfoFile[-which(is.na(FullInfoFile$FID)),]
	merged_temp_rest=merged_temp_rest[-which(is.na(FullInfoFile$FID)),]
}

numsubjects = length(FullInfoFile$IID);
VarNames=names(FullInfoFile)
columnnames = colnames(FullInfoFile);

drp=0

if (patients!=0) {
    writeLines(paste('STUDY DESIGN: There are patients in this study.'),con=zz,sep="\n")
    if (patients==1) {
        AffectionStatus=FullInfoFile$AffectionStatus
    } else {
        columnnames = colnames(FullInfoFile)
        AffectionStatus=FullInfoFile[,which(columnnames==patients)]
    }
 }

#Remove covariates with sd = 0 keeping the patient column if it exists
for (l in (Nset+1):length(VarNames)){
      columnnames = colnames(FullInfoFile);
    if (sd(as.numeric(FullInfoFile[,which(columnnames==VarNames[l])]))==0) {
        if (patients==1 && VarNames[l]== "AffectionStatus") {
			next
    	} else if (patients!=0 && VarNames[l] == patients ) {
    		next
		}
    	else {
        print(paste('The standard deviation of column', VarNames[l], 'is zero. Therefore, the column will be removed.'))
        columnnames = colnames(FullInfoFile)
        FullInfoFile=FullInfoFile[,-which(columnnames==VarNames[l])]
		drp=drp+1
		}
    }
}

## if related individuals, create connecting file for MERLIN-offline
if (related==1) {
    detach(merged_temp)
    attach(FullInfoFile)
        
    relatedness<-estimateRelatedness_2gen(FullInfoFile, sex=FullInfoFile$Sex)
    
    zyg=relatedness$zygosity
       
    connectingFile=data.frame(nFID=FID,nIID=IID,nPID=PID,nMID=MID,nSex=Sex,nZyg=zyg)
    if (is.na(as.numeric(PID[1])) == "TRUE"){
        connectingFilePID=data.frame(nFID=FID,nIID=PID,nPID=rep(0,numsubjects),nMID=rep(0,numsubjects),nSex=rep(1,numsubjects),nZyg=rep(0,numsubjects))
        connectingFileMID=data.frame(nFID=FID,nIID=MID,nPID=rep(0,numsubjects),nMID=rep(0,numsubjects),nSex=rep(2,numsubjects),nZyg=rep(0,numsubjects))
    } else {
        connectingFilePID=data.frame(nFID=FID,nIID=PID,nPID=rep("0",numsubjects),nMID=rep("0",numsubjects),nSex=rep(1,numsubjects),nZyg=rep(0,numsubjects))
        connectingFileMID=data.frame(nFID=FID,nIID=MID,nPID=rep("0",numsubjects),nMID=rep("0",numsubjects),nSex=rep(2,numsubjects),nZyg=rep(0,numsubjects))
    }

    names(connectingFile)<-c("nFID","nIID","nPID","nMID","nSex","nZyg")
    names(connectingFilePID)<-c("nFID","nIID","nPID","nMID","nSex","nZyg")
    names(connectingFileMID)<-c("nFID","nIID","nPID","nMID","nSex","nZyg")
    fullconnectingfile=rbind(connectingFile, connectingFilePID, connectingFileMID)
    Finalconnectingfile= fullconnectingfile[!duplicated(fullconnectingfile),]

for (i in 1:(dim(connectingFile))[1]) { 
	if (connectingFile[i,6] %% 2 ==1 ) {
	fam[which(fam[,2]==connectingFile[i,2]),6]= 1
	} else if (connectingFile[i,6]  ==0 ) {
	fam[which(fam[,2]==connectingFile[i,2]),6]= 0
	} else if (connectingFile[i,6]  %% 2 ==0 ) {
	fam[which(fam[,2]==connectingFile[i,2]),6]= 2
} }

	fam[which(fam[,6]==-9),6]=0
   names(fam)<-c("nFID","nIID","nPID","nMID","nSex","nZyg")
	fullFamfile=rbind(fam, connectingFilePID, connectingFileMID)
	fullFamfile= fullFamfile[!duplicated(fullFamfile[,2]),]
}
## if diseases exist, make one .dat file without the AffectionStatus covariates
###### if when removing patients, all healthy individuals have the same value for a covariate, remove that as a covariate too
FullInfoFile_healthy=FullInfoFile
FullInfoFile_disease=FullInfoFile
if (patients!=0) {

    FullInfoFile_healthy <- subset(FullInfoFile,AffectionStatus==0);
    if (nrow(subset(FullInfoFile,AffectionStatus==0)) > 0) {
    VarNames=colnames(FullInfoFile_healthy)
    columnnames = colnames(FullInfoFile_healthy);
    for (l in (Nset+1):length(VarNames)){
        columnnames = colnames(FullInfoFile_healthy);
        if (sd(FullInfoFile_healthy[,which(columnnames==VarNames[l])])==0) {
        print(paste('For healthy individuals only, the standard deviation of column', VarNames[l], 'is zero. Therefore, the column will be removed.'))
        FullInfoFile_healthy=FullInfoFile_healthy[,-which(columnnames==VarNames[l])]
        }
    }
    }
    FullInfoFile_disease <- subset(FullInfoFile,AffectionStatus==1);
    VarNames=colnames(FullInfoFile_disease)
    columnnames = colnames(FullInfoFile_disease);
    for (l in (Nset+1):length(VarNames)){
        columnnames = colnames(FullInfoFile_disease);
        if (sd(FullInfoFile_disease[,which(columnnames==VarNames[l])])==0) {
            print(paste('For disease individuals only, the standard deviation of column', VarNames[l], 'is zero. Therefore, the column will be removed.'))
            FullInfoFile_disease=FullInfoFile_disease[,-which(columnnames==VarNames[l])]
        }
    }
}

#Remove AffectionStatus if sd = 0

if ("AffectionStatus" %in% colnames(FullInfoFile) && patients==1 ){ 
	if (sd(sapply(FullInfoFile["AffectionStatus"], as.numeric))==0) {
		columnnames = colnames(FullInfoFile)
        FullInfoFile["AffectionStatus"] = NULL
		nCov=nCov-1
		drp=drp+1
	}
} else if ( patients!=0 ){ 
	if (sd(sapply(AffectionStatus, as.numeric))==0) {
		columnnames = colnames(FullInfoFile)
        FullInfoFile[,-which(columnnames==patients)] = NULL
		nCov=nCov-1
		drp=drp+1
	}
}

nVar_healthy=dim(FullInfoFile_healthy)[2]
nCov_healthy=nVar_healthy-Nset
nSub_healthy=dim(FullInfoFile_healthy)[1]

nVar_disease=dim(FullInfoFile_disease)[2]
nCov_disease=nVar_disease-Nset
nSub_disease=dim(FullInfoFile_disease)[1]

cat('    There are ',numsubjects,' total subjects\n')
cat('    There are ',nCov,' covariates for all subjects\n')
writeLines(paste('    There are ',numsubjects,' total subjects.'),con=zz,sep="\n")
writeLines(paste('    There are ',nCov,' covariates for all subjects.'),con=zz,sep="\n")
writeLines(paste('     -', cbind(colnames(FullInfoFile)[(Nset+1):nVar])),con=zz,sep="\n")


test=NULL
for (val in (Nids+1):(Nids+Nrois) ) {

A=as.numeric(FullInfoFile[,val])
if ( length(which(is.na(A))) > 0) {
	A=A[-which(is.na(A))]
}
test=rbind(test,cbind(colnames(FullInfoFile)[val],
(mean(as.numeric(A))),
(sd(as.numeric(A))),
(min(as.numeric(A))),
(max(as.numeric(A)))))
}
header=(c('ROI','Mean','SD','Min','Max'))
write.table(test,file=zz,quote=F,col.names=header,row.names=FALSE,sep = "\t");

test=NULL
for (val in (Nset+1):nVar ) {
A=as.numeric(FullInfoFile[,val])
if ( length(which(is.na(A))) > 0) {
	A=A[-which(is.na(A))]
}
test=rbind(test,cbind(colnames(FullInfoFile)[val],
(mean(as.numeric(A))),
(sd(as.numeric(A))),
(min(as.numeric(A))),
(max(as.numeric(A)))))
}
header=(c('Covariate','Mean','SD','Min','Max'))
write.table(test,file=zz,quote=F,col.names=header,row.names=FALSE,sep = "\t");

################ now print out the .dat and the .ped files
#
# Output names have been hard-coded for formatting of follow-up scripts.
#
#Write out ped file

if ( dim(FullInfoFile_healthy)[1] > 0 ){
write.table(FullInfoFile_healthy,paste(outFolder,"/ENIGMA_",eName,"_PEDfile_healthy.ped",sep=""),quote=F,col.names=F,row.names=F);
write.table(FullInfoFile_healthy,paste(outFolder,"/ENIGMA_",eName,"_PEDfile_wColNames_healthy.tbl",sep=""),quote=F,col.names=T,row.names=F);
write.table(colnames(FullInfoFile_healthy),paste(outFolder,"/ENIGMA_",eName,"_PEDfile_healthy.header",sep=""),quote=F,col.names=F,row.names=F);
write.table(cbind(c(rep("T",Nrois-36),rep("S",34),"C","S",rep("C",nCov_healthy)),c(colnames(FullInfoFile_healthy)[(Nids+1):nVar_healthy])),paste(outFolder,"/ENIGMA_",eName,"_DATfile_healthy.dat",sep=""),col.names=F,row.names=F,quote=F);
write.table(cbind(c(rep("S",Nrois-36),rep("T",34),"S","C",rep("C",nCov_healthy)),c(colnames(FullInfoFile_healthy)[(Nids+1):nVar_healthy])),paste(outFolder,"/ENIGMA_",eName,"_DATfile_healthy_wTHICK.dat",sep=""),col.names=F,row.names=F,quote=F);
write.table(cbind(c(rep("T",Nrois),rep("C",nCov_healthy)),c(colnames(FullInfoFile_healthy)[(Nids+1):nVar_healthy])),paste(outFolder,"/ENIGMA_",eName,"_DATfile_healthy_noSA.dat",sep=""),col.names=F,row.names=F,quote=F);
} else {
    cat('    There are no healthy individuals in this group. \n')
    writeLines(paste('    There are no healthy individuals in this group.'),con=zz,sep="\n")
}

#### print multiple dat files if patients
if (patients!=0) {
	if (nSub_healthy > 0) {
		cat('    There are ',nSub_healthy,' healthy subjects\n')
		cat('    There are ',nCov_healthy,' covariates for all healthy subjects\n')
		writeLines(paste('    There are ',nSub_healthy,' healthy subjects.'),con=zz,sep="\n")
		writeLines(paste('    There are ',nCov_healthy,' covariates for all healthy subjects.'),con=zz,sep="\n")
		writeLines(paste('     -', colnames(FullInfoFile_healthy)[(Nset+1):nVar_healthy]),con=zz,sep="\n")
		
		
		write.table(cbind(c(rep("T",Nrois-36),rep("S",34),"C","S",rep("C",nCov-drp)),c(colnames(FullInfoFile)[(Nids+1):(nVar-drp)])),paste(outFolder,"/ENIGMA_",eName,"_DATfile_fullGroup.dat",sep=""),col.names=F,row.names=F,quote=F);
		write.table(cbind(c(rep("S",Nrois-36),rep("T",34),"S","C",rep("C",nCov-drp)),c(colnames(FullInfoFile)[(Nids+1):(nVar-drp)])),paste(outFolder,"/ENIGMA_",eName,"_DATfile_fullGroup_wTHICK.dat",sep=""),col.names=F,row.names=F,quote=F);
		write.table(cbind(c(rep("T",Nrois),rep("C",nCov-drp)),c(colnames(FullInfoFile)[(Nids+1):(nVar-drp)])),paste(outFolder,"/ENIGMA_",eName,"_DATfile_fullGroup_noSA.dat",sep=""),col.names=F,row.names=F,quote=F);
		write.table(FullInfoFile,paste(outFolder,"/ENIGMA_",eName,"_PEDfile_fullGroup.ped",sep=""),quote=F,col.names=F,row.names=F);
		write.table(FullInfoFile,paste(outFolder,"/ENIGMA_",eName,"_PEDfile_wColNames_fullGroup.tbl",sep=""),quote=F,col.names=T,row.names=F);
		write.table(colnames(FullInfoFile),paste(outFolder,"/ENIGMA_",eName,"_PEDfile_fullGroup.header",sep=""),quote=F,col.names=F,row.names=F);
	}
	
	if (nSub_healthy == 0) {
		cat('    There are ',nSub_disease,' patients\n')
		cat('    There are ',nCov_disease,' covariates for all patients\n')
		writeLines(paste('    There are ',nSub_disease,' patients.'),con=zz,sep="\n")
		writeLines(paste('    There are ',nCov_disease,' covariates for all patients.'),con=zz,sep="\n")
		writeLines(paste('     -', colnames(FullInfoFile_disease)[(Nset+1):nVar_disease]),con=zz,sep="\n")
		
		write.table(FullInfoFile_disease,paste(outFolder,"/ENIGMA_",eName,"_PEDfile_patients.ped",sep=""),quote=F,col.names=F,row.names=F);
		write.table(FullInfoFile_disease,paste(outFolder,"/ENIGMA_",eName,"_PEDfile_wColNames_patients.tbl",sep=""),quote=F,col.names=T,row.names=F);
		write.table(colnames(FullInfoFile_disease),paste(outFolder,"/ENIGMA_",eName,"_PEDfile_patients.header",sep=""),quote=F,col.names=F,row.names=F);
		write.table(cbind(c(rep("T",Nrois-36),rep("S",34),"C","S",rep("C",nCov_disease)),c(colnames(FullInfoFile_disease)[(Nids+1):nVar_disease])),paste(outFolder,"/ENIGMA_",eName,"_DATfile_patients.dat",sep=""),col.names=F,row.names=F,quote=F);
		write.table(cbind(c(rep("S",Nrois-36),rep("T",34),"S","C",rep("C",nCov_disease)),c(colnames(FullInfoFile_disease)[(Nids+1):nVar_disease])),paste(outFolder,"/ENIGMA_",eName,"_DATfile_patients_wTHICK.dat",sep=""),col.names=F,row.names=F,quote=F);
		write.table(cbind(c(rep("T",Nrois),rep("C",nCov_disease)),c(colnames(FullInfoFile_disease)[(Nids+1):nVar_disease])),paste(outFolder,"/ENIGMA_",eName,"_DATfile_patients_noSA.dat",sep=""),col.names=F,row.names=F,quote=F);
	}
}
#### print connecting files if related for MERLIN-offline
if (related==1) {
 	write.table(fam,paste(outFolder,"/ENIGMA_",eName,"_connecting.fam",sep=""),col.names=F,row.names=F,quote=F);
 	write.table(fullFamfile,paste(outFolder,"/ENIGMA_",eName,"_connecting.full.fam",sep=""),col.names=F,row.names=F,quote=F);
    write.table(Finalconnectingfile,paste(outFolder,"/ENIGMA_",eName,"_connecting.ped",sep=""),col.names=F,row.names=F,quote=F);
    write.table(cbind("Z","zygosity"),paste(outFolder,"/ENIGMA_",eName,"_connecting.dat",sep=""),quote=F,col.names=F,row.names=F);
}

cat('****** DONE ****** Files have been written to ',outFolder,'\n')
writeLines(paste('****** DONE ****** Files have been written. '),con=zz,sep="\n")
close(zz)

###########################################################
###########################################################
