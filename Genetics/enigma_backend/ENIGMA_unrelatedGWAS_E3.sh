#!/bin/bash

# * Neda Jahanshad - neda.jahanshad@ini.usc.edu
# * ENIGMA_DTI 2014.
# */

version=20151222

######## SH INPUTS #############
run_machdir=${1}  #give the directory to the imputed output from Mach/minimac
machdir=${2}
peddatdir=${3}  #give the dir to the ped and dat files just created
samplename=${4}  #give abbreviated name of your sample, no spaces in the name (i.e. ADNI)
mach2qtlout=${5}  #make a folder for the output from mach2qtl
status=${6}  ### this should be H for healthy, D for disease or HD for healthy and disease
Nnodes=${7}  # can split up the processing into N nodes
totalFiles=${8}
eName=${9}
SGE_TASK_ID=${10}
mode=${11}
######## END SH INPUTS #########


#ls -1 ${machdir}/chunk*-ready4mach.*.imputed.dose.gz > ${mach2qtlout}/fileList.txt
#totalFiles=`ls ${machdir}/chunk*-ready4mach.*.imputed.dose.gz |wc -w`

NchunksPerTask=$((totalFiles/Nnodes))
start_pt=$(($((${SGE_TASK_ID}-1))*${NchunksPerTask}+1))
end_pt=$((${SGE_TASK_ID}*${NchunksPerTask}))

if [ "$end_pt" == "$totalFiles" ]
then
end_pt=$((${totalFiles}))
fi

if [ "$mode" != "run" ]; then
	echo "Running in MANUAL mode"	
	echo "Commands will be stored in the text files called Step1_Manual_GWAS.txt and Step2_Manual_GZIP.txt"
fi

case $status in

    H)

    for ((i=${start_pt}; i<=${end_pt};i++));
        do
        fileDose=$(awk -v "line=$i" 'NR == line' ${mach2qtlout}/fileList_${SGE_TASK_ID}.txt )
        fileInfo=`echo ${fileDose%.dose.gz}.info.gz`
        chr=$(basename  ${fileInfo} | awk -F '.' '{print $2}')
        chunk=$(basename  ${fileInfo} | awk -F '-' '{print $1}')

        outName=${mach2qtlout}/${samplename}_${eName}_healthy_${chr}_${chunk}_${version}.out
		wTHICKOutName=${mach2qtlout}/${samplename}_${eName}_healthy_wTHICK_${chr}_${chunk}_${version}.out
		noSAOutName=${mach2qtlout}/${samplename}_${eName}_healthy_noSA_${chr}_${chunk}_${version}.out
        datFileName=${peddatdir}/ENIGMA_${eName}_DATfile_healthy.dat
		wTHICKDatFileName=${peddatdir}/ENIGMA_${eName}_DATfile_healthy_wTHICK.dat
		noSADatFileName=${peddatdir}/ENIGMA_${eName}_DATfile_healthy_noSA.dat
        pedFileName=${peddatdir}/ENIGMA_${eName}_PEDfile_healthy.ped
		
		if [ "$mode" == "run" ]; then
        	echo "${run_machdir}/executables/mach2qtl --datfile ${datFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${outName}"
        	${run_machdir}/executables/mach2qtl --datfile ${datFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${outName}
        	gzip -f ${outName}
			echo "${run_machdir}/executables/mach2qtl --datfile ${wTHICKDatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${wTHICKOutName}"
        	${run_machdir}/executables/mach2qtl --datfile ${wTHICKDatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${wTHICKOutName}
        	gzip -f ${wTHICKOutName}
			echo "${run_machdir}/executables/mach2qtl --datfile ${noSADatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${noSAOutName}"
			${run_machdir}/executables/mach2qtl --datfile ${noSADatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${noSAOutName}
			gzip -f ${noSAOutName}
		else
			echo "${run_machdir}/executables/mach2qtl --datfile ${datFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${outName}" >> Step1_Manual_GWAS.txt
			echo "gzip -f ${outName}" >> Step2_Manual_GZIP.txt
			echo "${run_machdir}/executables/mach2qtl --datfile ${wTHICKDatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${wTHICKOutName}" >> Step1_Manual_GWAS.txt
			echo "gzip -f ${wTHICKOutName}" >> Step2_Manual_GZIP.txt
			echo "${run_machdir}/executables/mach2qtl --datfile ${noSADatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${noSAOutName}" >> Step1_Manual_GWAS.txt
			echo "gzip -f ${noSAOutName}" >> Step2_Manual_GZIP.txt
		fi
	done
    ;;

    D)
    for ((i=${start_pt}; i<=${end_pt};i++));
        do
        fileDose=$(awk -v "line=$i" 'NR == line' ${mach2qtlout}/fileList_${SGE_TASK_ID}.txt )
        fileInfo=`echo ${fileDose%.dose.gz}.info.gz`
        chr=$(basename  ${fileInfo} | awk -F '.' '{print $2}')
        chunk=$(basename  ${fileInfo} | awk -F '-' '{print $1}')

        outName=${mach2qtlout}/${samplename}_${eName}_disease_${chr}_${chunk}_${version}.out
		wTHICKOutName=${mach2qtlout}/${samplename}_${eName}_disease_wTHICK_${chr}_${chunk}_${version}.out
		noSAOutName=${mach2qtlout}/${samplename}_${eName}_disease_noSA_${chr}_${chunk}_${version}.out
        datFileName=${peddatdir}/ENIGMA_${eName}_DATfile_patients.dat
		wTHICKDatFileName=${peddatdir}/ENIGMA_${eName}_DATfile_patients_wTHICK.dat
		noSADatFileName=${peddatdir}/ENIGMA_${eName}_DATfile_patients_noSA.dat
        pedFileName=${peddatdir}/ENIGMA_${eName}_PEDfile_patients.ped
		if [ "$mode" == "run" ]; then
       	 	echo "${run_machdir}/executables/mach2qtl --datfile ${datFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${outName}"
       	 	${run_machdir}/executables/mach2qtl --datfile ${datFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${outName}
       	 	gzip -f ${outName}
			echo "${run_machdir}/executables/mach2qtl --datfile ${wTHICKDatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${wTHICKOutName}"
       	 	${run_machdir}/executables/mach2qtl --datfile ${wTHICKDatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${wTHICKOutName}
       	 	gzip -f ${wTHICKOutName}
			echo "${run_machdir}/executables/mach2qtl --datfile ${noSADatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${noSAOutName}"
       	 	${run_machdir}/executables/mach2qtl --datfile ${noSADatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${noSAOutName}
       	 	gzip -f ${noSAOutName}
		else
			echo "${run_machdir}/executables/mach2qtl --datfile ${datFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${outName}" >> Step1_Manual_GWAS.txt
			echo "gzip -f ${outName}" >> Step2_Manual_GZIP.txt
			echo "${run_machdir}/executables/mach2qtl --datfile ${wTHICKDatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${wTHICKOutName}" >> Step1_Manual_GWAS.txt
			echo "gzip -f ${wTHICKOutName}" >> Step2_Manual_GZIP.txt
			echo "${run_machdir}/executables/mach2qtl --datfile ${saDatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${noSAOutName}" >> Step1_Manual_GWAS.txt
			echo "gzip -f ${noSAOutName}" >> Step2_Manual_GZIP.txt
		fi
    done
    ;;

    HD)
# we are actually going to double up here, and run 2x as many runs per node one for healthy only and one for the full group
    for ((i=${start_pt}; i<=${end_pt};i++));
        do
        fileDose=$(awk -v "line=$i" 'NR == line' ${mach2qtlout}/fileList_${SGE_TASK_ID}.txt )
        fileInfo=`echo ${fileDose%.dose.gz}.info.gz`
        chr=$(basename  ${fileInfo} | awk -F '.' '{print $2}')
        chunk=$(basename  ${fileInfo} | awk -F '-' '{print $1}')

        ###### run healthy and disease -- full group
        outName=${mach2qtlout}/${samplename}_${eName}_mixedHD_${chr}_${chunk}_${version}.out
		wTHICKOutName=${mach2qtlout}/${samplename}_${eName}_mixedHD_wTHICK_${chr}_${chunk}_${version}.out
		noSAOutName=${mach2qtlout}/${samplename}_${eName}_mixedHD_noSA_${chr}_${chunk}_${version}.out
        datFileName=${peddatdir}/ENIGMA_${eName}_DATfile_fullGroup.dat
		wTHICKDatFileName=${peddatdir}/ENIGMA_${eName}_DATfile_fullGroup_wTHICK.dat
		noSADatFileName=${peddatdir}/ENIGMA_${eName}_DATfile_fullGroup_noSA.dat
        pedFileName=${peddatdir}/ENIGMA_${eName}_PEDfile_fullGroup.ped
		if [ "$mode" == "run" ]; then
        	echo "${run_machdir}/executables/mach2qtl --datfile ${datFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${outName}"
        	${run_machdir}/executables/mach2qtl --datfile ${datFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${outName}
        	gzip -f ${outName}
			echo "${run_machdir}/executables/mach2qtl --datfile ${wTHICKDatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${wTHICKOutName}"
        	${run_machdir}/executables/mach2qtl --datfile ${wTHICKDatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${wTHICKOutName}
        	gzip -f ${wTHICKOutName}
			echo "${run_machdir}/executables/mach2qtl --datfile ${noSADatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${noSAOutName}"
        	${run_machdir}/executables/mach2qtl --datfile ${noSADatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${noSAOutName}
        	gzip -f ${noSAOutName}
		else
			echo "${run_machdir}/executables/mach2qtl --datfile ${datFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${outName}" >> Step1_Manual_GWAS.txt
			echo "gzip -f ${outName}" >> Step2_Manual_GZIP.txt
			echo "${run_machdir}/executables/mach2qtl --datfile ${wTHICKDatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${wTHICKOutName}" >> Step1_Manual_GWAS.txt
			echo "gzip -f ${wTHICKOutName}" >> Step2_Manual_GZIP.txt
			echo "${run_machdir}/executables/mach2qtl --datfile ${noSADatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${noSAOutName}" >> Step1_Manual_GWAS.txt
			echo "gzip -f ${noSAOutName}" >> Step2_Manual_GZIP.txt	
		fi

        ###### run healthy only
        outName=${mach2qtlout}/${samplename}_${eName}_healthy_${chr}_${chunk}_${version}.out
		wTHICKOutName=${mach2qtlout}/${samplename}_${eName}_healthy_wTHICK_${chr}_${chunk}_${version}.out
		noSAOutName=${mach2qtlout}/${samplename}_${eName}_healthy_noSA_${chr}_${chunk}_${version}.out
        datFileName=${peddatdir}/ENIGMA_${eName}_DATfile_healthy.dat
		wTHICKDatFileName=${peddatdir}/ENIGMA_${eName}_DATfile_healthy_wTHICK.dat
		noSADatFileName=${peddatdir}/ENIGMA_${eName}_DATfile_healthy_noSA.dat
        pedFileName=${peddatdir}/ENIGMA_${eName}_PEDfile_healthy.ped
		if [ "$mode" == "run" ]; then
        	echo "${run_machdir}/executables/mach2qtl --datfile ${datFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${outName}"
        	${run_machdir}/executables/mach2qtl --datfile ${datFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${outName}
        	gzip -f ${outName}
			echo "${run_machdir}/executables/mach2qtl --datfile ${wTHICKDatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${wTHICKOutName}"
        	${run_machdir}/executables/mach2qtl --datfile ${wTHICKDatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${wTHICKOutName}
        	gzip -f ${wTHICKOutName}
			echo "${run_machdir}/executables/mach2qtl --datfile ${noSADatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${noSAOutName}"
        	${run_machdir}/executables/mach2qtl --datfile ${noSADatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${noSAOutName}
        	gzip -f ${noSAOutName}
		else
			echo "${run_machdir}/executables/mach2qtl --datfile ${datFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${outName}" >> Step1_Manual_GWAS.txt
			echo "gzip -f ${outName}" >> Step2_Manual_GZIP.txt
			echo "${run_machdir}/executables/mach2qtl --datfile ${wTHICKDatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${wTHICKOutName}" >> Step1_Manual_GWAS.txt
			echo "gzip -f ${wTHICKOutName}" >> Step2_Manual_GZIP.txt
			echo "${run_machdir}/executables/mach2qtl --datfile ${noSADatFileName} --pedfile ${pedFileName} --infofile ${fileInfo} --dosefile ${fileDose} --samplesize > ${noSAOutName}" >> Step1_Manual_GWAS.txt
			echo "gzip -f ${noSAOutName}" >> Step2_Manual_GZIP.txt
		fi
	done
    ;;

esac
