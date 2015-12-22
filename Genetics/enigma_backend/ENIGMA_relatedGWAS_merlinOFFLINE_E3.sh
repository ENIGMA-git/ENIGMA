#!/bin/bash

# * Neda Jahanshad - neda.jahanshad@ini.usc.edu
# * ENIGMA_DTI 2014.
# */

version=20151222

######## SH INPUTS #############
run_merlin=${1}
genodir=${2}  #give the directory to the imputed output from Mach/minimac
peddatdir=${3}  #give the dir to the ped and dat files just created
samplename=${4}  #give abbreviated name of your sample, no spaces in the name (i.e. TWINS)
merlinout=${5}  #make a folder for the output from mach2qtl
status=${6}  ### this should be H for healthy, D for disease or HD for healthy and disease
Nnodes=${7}  # can split up the processing into N nodes
totalFiles=${8}
eName=${9}
merlinOFFLINE=${10}
SGE_TASK_ID=${11}
mode=${12}
######## END SH INPUTS #########

ped_connect=${peddatdir}/ENIGMA_${eName}_connecting.full.fam  ## needs to match up with mac2merlin files
dat_connect=${peddatdir}/ENIGMA_${eName}_connecting.dat
fam_connect=${peddatdir}/ENIGMA_${eName}_connecting.fam ## needs to be from imputed files!

#ls -1 ${genodir}/chunk*-ready4mach.*.imputed.infer.dat.gz > ${merlinout}/fileList.txt
#totalFiles=`ls ${genodir}/chunk*-ready4mach.*.imputed.infer.dat.gz |wc -w`

if [ "$mode" != "run" ]; then
	echo "Running in MANUAL mode"	
	echo "Commands will be stored in the text files called Step1_Manual_GWAS.txt and Step2_Manual_GZIP.txt"
fi


NchunksPerTask=$((totalFiles/Nnodes+1))
echo $NchunksPerTask
start_pt=$(($((${SGE_TASK_ID}-1))*${NchunksPerTask}+1))
end_pt=$((${SGE_TASK_ID}*${NchunksPerTask}))

if [ "$end_pt" == "$((totalFiles + 1))" ] ###
then
end_pt=$((${totalFiles}))
fi

case $status in

    H)

    for ((i=${start_pt}; i<=${end_pt};i++));
        do
        if [ ${merlinOFFLINE} -eq 0 ]  ##
        then
        fileDose=$(awk -v "line=$i" 'NR == line' ${merlinout}/fileList_${SGE_TASK_ID}.txt )
        filePrefix=`echo ${fileDose%.dose.gz}`
        noPathPrefix=${fileDose##*/}
        noPathPrefix=`echo ${noPathPrefix%.dose.gz}`
        
        if [ "$mode" == "run" ]; then
		${run_merlin}/1KGPminimac2merlin.pl -fam ${fam_connect} -prefix ${filePrefix} -out ${genodir}/${noPathPrefix}
	else
		echo "${run_merlin}/1KGPminimac2merlin.pl -fam ${fam_connect} -prefix ${filePrefix} -out ${genodir}/${noPathPrefix}" >> Step1_Manual_GWAS.txt
	fi
       
        fileDat=${genodir}/${noPathPrefix}.dat.gz
        else
        fileDat=$(awk -v "line=$i" 'NR == line' ${merlinout}/fileList_${SGE_TASK_ID}.txt )
        filePrefix=`echo ${fileDat%.dat.gz}`
        fi
        
        echo "done with minimac2merlin" $filePrefix

        
        fileMap=`echo ${fileDat%.dat.gz}.map.gz`
        fileFreq=`echo ${fileDat%.dat.gz}.freq.gz`
        filePed=`echo ${fileDat%.dat.gz}.ped.gz`

        chr=`basename  ${filePrefix} | awk -F '.' '{print $2}'`
        chunk=$(basename  ${filePrefix} | awk -F '-' '{print $1}')
        
        echo ${chr}
        echo ${chunk}

        outName=${merlinout}/${samplename}_${eName}_healthy_${chr}_${chunk}_${version}.out
		wTHICKOutName=${merlinout}/${samplename}_${eName}_healthy_wTHICK_${chr}_${chunk}_${version}.out
		noSAOutName=${merlinout}/${samplename}_${eName}_healthy_noSA_${chr}_${chunk}_${version}.out
        outTBLName=${merlinout}/${samplename}_${eName}_healthy_${chr}_${chunk}.tbl
		wTHICKOutTBLName=${merlinout}/${samplename}_${eName}_healthy_wTHICK_${chr}_${chunk}.tbl
		noSAOutTBLName=${merlinout}/${samplename}_${eName}_healthy_noSA_${chr}_${chunk}.tbl
        outNamePrefix=${merlinout}/${samplename}_${eName}_healthy_${chr}_${chunk}
		wTHICKOutNamePrefix=${merlinout}/${samplename}_${eName}_healthy_wTHICK_${chr}_${chunk}
		noSAOutNamePrefix=${merlinout}/${samplename}_${eName}_healthy_noSA_${chr}_${chunk}
        datFileName=${peddatdir}/ENIGMA_${eName}_DATfile_healthy.dat
		wTHICKDatFileName=${peddatdir}/ENIGMA_${eName}_DATfile_healthy_wTHICK.dat
		noSADatFileName=${peddatdir}/ENIGMA_${eName}_DATfile_healthy_noSA.dat
        pedFileName=${peddatdir}/ENIGMA_${eName}_PEDfile_healthy.ped


        if [ "$mode" == "run" ]; then
        	echo ${outName}
       		${run_merlin}/merlin-1.1.2/executables/merlin-offline -m ${fileMap} -f ${fileFreq} --pedinfer ${filePed} --datinfer ${fileDat} -p ${ped_connect},${pedFileName} -d ${dat_connect},${datFileName} --useCovariates --tabulate  --prefix ${outNamePrefix} > ${outName}
       		gzip ${outTBLName}
        	gzip ${outName}
			echo ${wTHICKOutName}
       		${run_merlin}/merlin-1.1.2/executables/merlin-offline -m ${fileMap} -f ${fileFreq} --pedinfer ${filePed} --datinfer ${fileDat} -p ${ped_connect},${pedFileName} -d ${dat_connect},${wTHICKDatFileName} --useCovariates --tabulate  --prefix ${wTHICKOutNamePrefix} > ${wTHICKOutName}
       		gzip ${wTHICKOutTBLName}
        	gzip ${wTHICKOutName}
			echo ${noSAOutName}
       		${run_merlin}/merlin-1.1.2/executables/merlin-offline -m ${fileMap} -f ${fileFreq} --pedinfer ${filePed} --datinfer ${fileDat} -p ${ped_connect},${pedFileName} -d ${dat_connect},${noSADatFileName} --useCovariates --tabulate  --prefix ${noSAOutNamePrefix} > ${noSAOutName}
       		gzip ${noSAOutTBLName}
        	gzip ${noSAOutName}
	
	else
		echo "${run_merlin}/merlin-1.1.2/executables/merlin-offline -m ${fileMap} -f ${fileFreq} --pedinfer ${filePed} --datinfer ${fileDat} -p ${ped_connect},${pedFileName} -d ${dat_connect},${datFileName} --useCovariates --tabulate  --prefix ${outNamePrefix} > ${outName}" >> Step1_Manual_GWAS.txt
		echo "gzip ${outTBLName}" >> Step2_Manual_GZIP.txt
		echo "gzip ${outName}" >> Step2_Manual_GZIP.txt
		echo "${run_merlin}/merlin-1.1.2/executables/merlin-offline -m ${fileMap} -f ${fileFreq} --pedinfer ${filePed} --datinfer ${fileDat} -p ${ped_connect},${pedFileName} -d ${dat_connect},${wTHICKDatFileName} --useCovariates --tabulate  --prefix ${wTHICKOutNamePrefix} > ${wTHICKOutName}" >> Step1_Manual_GWAS.txt
		echo "gzip ${wTHICKOutTBLName}" >> Step2_Manual_GZIP.txt
		echo "gzip ${wTHICKOutName}" >> Step2_Manual_GZIP.txt
		echo "${run_merlin}/merlin-1.1.2/executables/merlin-offline -m ${fileMap} -f ${fileFreq} --pedinfer ${filePed} --datinfer ${fileDat} -p ${ped_connect},${pedFileName} -d ${dat_connect},${noSADatFileName} --useCovariates --tabulate  --prefix ${noSAOutNamePrefix} > ${noSAOutName}" >> Step1_Manual_GWAS.txt
		echo "gzip ${noSAOutTBLName}" >> Step2_Manual_GZIP.txt
		echo "gzip ${noSAOutName}" >> Step2_Manual_GZIP.txt
	fi
        
        
        
    done
    ;;

    HD)
# we are actually going to double up here, and run 2x as many runs per node one for healthy only and one for the full group
    for ((i=${start_pt}; i<=${end_pt};i++));
        do
        if [ ${merlinOFFLINE} -eq 0 ]  ##
        then
        fileDose=$(awk -v "line=$i" 'NR == line' ${merlinout}/fileList_${SGE_TASK_ID}.txt )
        filePrefix=`echo ${fileDose%.dose.gz}`
        noPathPrefix=${fileDose##*/}
        noPathPrefix=`echo ${noPathPrefix%.dose.gz}`
        
        if [ "$mode" == "run" ]; then
		${run_merlin}/1KGPminimac2merlin.pl -fam ${fam_connect} -prefix ${filePrefix} -out ${genodir}/${noPathPrefix}
	else
		echo "${run_merlin}/1KGPminimac2merlin.pl -fam ${fam_connect} -prefix ${filePrefix} -out ${genodir}/${noPathPrefix}" >> Step1_Manual_GWAS.txt
	fi
        
        
        fileDat=${genodir}/${noPathPrefix}.dat.gz
        else
        fileDat=$(awk -v "line=$i" 'NR == line' ${merlinout}/fileList_${SGE_TASK_ID}.txt )
        filePrefix=`echo ${fileDat%.dat.gz}`
        fi
	
        fileMap=`echo ${fileDat%.dat.gz}.map.gz`
        fileFreq=`echo ${fileDat%.dat.gz}.freq.gz`
        filePed=`echo ${fileDat%.dat.gz}.ped.gz`

        chr=$(basename  ${filePrefix} | awk -F '.' '{print $2}')
        chunk=$(basename  ${filePrefix} | awk -F '-' '{print $1}')

        ###### run healthy and disease -- full group
        outName=${merlinout}/${samplename}_${eName}_mixedHD_${chr}_${chunk}_${version}.out
		wTHICKOutName=${merlinout}/${samplename}_${eName}_mixedHD_wTHICK_${chr}_${chunk}_${version}.out
        noSAOutName=${merlinout}/${samplename}_${eName}_mixedHD_noSA_${chr}_${chunk}_${version}.out
        outTBLName=${merlinout}/${samplename}_${eName}_mixedHD_${chr}_${chunk}.tbl
		wTHICKOutTBLName=${merlinout}/${samplename}_${eName}_mixedHD_wTHICK_${chr}_${chunk}.tbl
        noSAOutTBLName=${merlinout}/${samplename}_${eName}_mixedHD_noSA_${chr}_${chunk}.tbl
        outNamePrefix=${merlinout}/${samplename}_${eName}_mixedHD_${chr}_${chunk}
		wTHICKOutNamePrefix=${merlinout}/${samplename}_${eName}_mixedHD_wTHICK_${chr}_${chunk}
		noSAOutNamePrefix=${merlinout}/${samplename}_${eName}_mixedHD_noSA_${chr}_${chunk}

        datFileName=${peddatdir}/ENIGMA_${eName}_DATfile_fullGroup.dat
		wTHICKDatFileName=${peddatdir}/ENIGMA_${eName}_DATfile_fullGroup_wTHICK.dat
		noSADatFileName=${peddatdir}/ENIGMA_${eName}_DATfile_fullGroup_noSA.dat
        pedFileName=${peddatdir}/ENIGMA_${eName}_PEDfile_fullGroup.ped

        if [ "$mode" == "run" ]; then
        	echo ${outName}
       		${run_merlin}/merlin-1.1.2/executables/merlin-offline -m ${fileMap} -f ${fileFreq} --pedinfer ${filePed} --datinfer ${fileDat} -p ${ped_connect},${pedFileName} -d ${dat_connect},${datFileName} --useCovariates --tabulate  --prefix ${outNamePrefix} > ${outName}
       		gzip ${outTBLName}
        	gzip ${outName}
			echo ${wTHICKOutName}
       		${run_merlin}/merlin-1.1.2/executables/merlin-offline -m ${fileMap} -f ${fileFreq} --pedinfer ${filePed} --datinfer ${fileDat} -p ${ped_connect},${pedFileName} -d ${dat_connect},${wTHICKDatFileName} --useCovariates --tabulate  --prefix ${wTHICKOutNamePrefix} > ${wTHICKOutName}
       		gzip ${wTHICKOutTBLName}
        	gzip ${wTHICKOutName}
			echo ${noSAOutName}
       		${run_merlin}/merlin-1.1.2/executables/merlin-offline -m ${fileMap} -f ${fileFreq} --pedinfer ${filePed} --datinfer ${fileDat} -p ${ped_connect},${pedFileName} -d ${dat_connect},${noSADatFileName} --useCovariates --tabulate  --prefix ${noSAOutNamePrefix} > ${noSAOutName}
       		gzip ${noSAOutTBLName}
        	gzip ${noSAOutName}
	
	else
		echo "${run_merlin}/merlin-1.1.2/executables/merlin-offline -m ${fileMap} -f ${fileFreq} --pedinfer ${filePed} --datinfer ${fileDat} -p ${ped_connect},${pedFileName} -d ${dat_connect},${datFileName} --useCovariates --tabulate  --prefix ${outNamePrefix} > ${outName}" >> Step1_Manual_GWAS.txt
		echo "gzip ${outTBLName}" >> Step2_Manual_GZIP.txt
		echo "gzip ${outName}" >> Step2_Manual_GZIP.txt
		echo "${run_merlin}/merlin-1.1.2/executables/merlin-offline -m ${fileMap} -f ${fileFreq} --pedinfer ${filePed} --datinfer ${fileDat} -p ${ped_connect},${pedFileName} -d ${dat_connect},${wTHICKDatFileName} --useCovariates --tabulate  --prefix ${wTHICKOutNamePrefix} > ${wTHICKOutName}" >> Step1_Manual_GWAS.txt
		echo "gzip ${wTHICKOutTBLName}" >> Step2_Manual_GZIP.txt
		echo "gzip ${wTHICKOutName}" >> Step2_Manual_GZIP.txt
		echo "${run_merlin}/merlin-1.1.2/executables/merlin-offline -m ${fileMap} -f ${fileFreq} --pedinfer ${filePed} --datinfer ${fileDat} -p ${ped_connect},${pedFileName} -d ${dat_connect},${noSADatFileName} --useCovariates --tabulate  --prefix ${noSAOutNamePrefix} > ${noSAOutName}" >> Step1_Manual_GWAS.txt
		echo "gzip ${noSAOutTBLName}" >> Step2_Manual_GZIP.txt
		echo "gzip ${noSAOutName}" >> Step2_Manual_GZIP.txt
	fi

        ###### run healthy only
        outName=${merlinout}/${samplename}_${eName}_healthy_${chr}_${chunk}_${version}.out
		wTHICKOutName=${merlinout}/${samplename}_${eName}_healthy_wTHICK_${chr}_${chunk}_${version}.out
        noSAOutName=${merlinout}/${samplename}_${eName}_healthy_noSA_${chr}_${chunk}_${version}.out
        outTBLName=${merlinout}/${samplename}_${eName}_healthy_${chr}_${chunk}.tbl
		wTHICKOutTBLName=${merlinout}/${samplename}_${eName}_healthy_wTHICK_${chr}_${chunk}.tbl
        noSAOutTBLName=${merlinout}/${samplename}_${eName}_healthy_noSA_${chr}_${chunk}.tbl
        outNamePrefix=${merlinout}/${samplename}_${eName}_healthy_${chr}_${chunk}
		wTHICKOutNamePrefix=${merlinout}/${samplename}_${eName}_healthy_wTHICK_${chr}_${chunk}
        noSAOutNamePrefix=${merlinout}/${samplename}_${eName}_healthy_noSA_${chr}_${chunk}
        datFileName=${peddatdir}/ENIGMA_${eName}_DATfile_healthy.dat
		wTHICKDatFileName=${peddatdir}/ENIGMA_${eName}_DATfile_healthy_wTHICK.dat
        noSADatFileName=${peddatdir}/ENIGMA_${eName}_DATfile_healthy_noSA.dat
        pedFileName=${peddatdir}/ENIGMA_${eName}_PEDfile_healthy.ped

        if [ "$mode" == "run" ]; then
        	echo ${outName}
       		${run_merlin}/merlin-1.1.2/executables/merlin-offline -m ${fileMap} -f ${fileFreq} --pedinfer ${filePed} --datinfer ${fileDat} -p ${ped_connect},${pedFileName} -d ${dat_connect},${datFileName} --useCovariates --tabulate  --prefix ${outNamePrefix} > ${outName}
       		gzip ${outTBLName}
        	gzip ${outName}
			echo ${wTHICKOutName}
       		${run_merlin}/merlin-1.1.2/executables/merlin-offline -m ${fileMap} -f ${fileFreq} --pedinfer ${filePed} --datinfer ${fileDat} -p ${ped_connect},${pedFileName} -d ${dat_connect},${wTHICKDatFileName} --useCovariates --tabulate  --prefix ${wTHICKOutNamePrefix} > ${wTHICKOutName}
       		gzip ${wTHICKOutTBLName}
        	gzip ${wTHICKOutName}
			echo ${noSAOutName}
       		${run_merlin}/merlin-1.1.2/executables/merlin-offline -m ${fileMap} -f ${fileFreq} --pedinfer ${filePed} --datinfer ${fileDat} -p ${ped_connect},${pedFileName} -d ${dat_connect},${noSADatFileName} --useCovariates --tabulate  --prefix ${noSAOutNamePrefix} > ${noSAOutName}
       		gzip ${noSAOutTBLName}
        	gzip ${noSAOutName}
	
	else
		echo "${run_merlin}/merlin-1.1.2/executables/merlin-offline -m ${fileMap} -f ${fileFreq} --pedinfer ${filePed} --datinfer ${fileDat} -p ${ped_connect},${pedFileName} -d ${dat_connect},${datFileName} --useCovariates --tabulate  --prefix ${outNamePrefix} > ${outName}" >> Step1_Manual_GWAS.txt
		echo "gzip ${outTBLName}" >> Step2_Manual_GZIP.txt
		echo "gzip ${outName}" >> Step2_Manual_GZIP.txt
		echo "${run_merlin}/merlin-1.1.2/executables/merlin-offline -m ${fileMap} -f ${fileFreq} --pedinfer ${filePed} --datinfer ${fileDat} -p ${ped_connect},${pedFileName} -d ${dat_connect},${wTHICKDatFileName} --useCovariates --tabulate  --prefix ${wTHICKOutNamePrefix} > ${wTHICKOutName}" >> Step1_Manual_GWAS.txt
		echo "gzip ${wTHICKOutTBLName}" >> Step2_Manual_GZIP.txt
		echo "gzip ${wTHICKOutName}" >> Step2_Manual_GZIP.txt
		echo "${run_merlin}/merlin-1.1.2/executables/merlin-offline -m ${fileMap} -f ${fileFreq} --pedinfer ${filePed} --datinfer ${fileDat} -p ${ped_connect},${pedFileName} -d ${dat_connect},${noSADatFileName} --useCovariates --tabulate  --prefix ${noSAOutNamePrefix} > ${noSAOutName}" >> Step1_Manual_GWAS.txt
		echo "gzip ${noSAOutTBLName}" >> Step2_Manual_GZIP.txt
		echo "gzip ${noSAOutName}" >> Step2_Manual_GZIP.txt
	fi
	
    done
    ;;
esac





