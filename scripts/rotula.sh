# *************************************************************** #
#                                                                 #
# This script verifies the experimental data at 'data' folder and #
# creates files in 'path' folder describing all the observed      #
# data.                                                           #
#                                                                 #
# ***************************************************************
#!/bin/bash

DAT=../data/ # directory with experimental data
SPK=/spikes/01/ # path inside the data folder
DPATH=../path/ # output folder

# Gets the mouses.
MOUSE=`ls ${DAT}`
MOUSE=`echo "${MOUSE}" | grep "ge" | tr -d 'ge' | sort -n`

# Creates the files at output folder.
echo "${MOUSE}" > ${DPATH}mouses.dat
> ${DPATH}HPmousesP1.dat
> ${DPATH}HPmousesP3.dat
> ${DPATH}HPCA1mousesP1.dat
> ${DPATH}HPCA1mousesP3.dat
> ${DPATH}HPDGmousesP1.dat
> ${DPATH}HPDGmousesP3.dat
> ${DPATH}S1mousesP1.dat
> ${DPATH}S1mousesP3.dat
> ${DPATH}V1mousesP1.dat
> ${DPATH}V1mousesP3.dat
> ${DPATH}dataPath.dat

for RAT in ${MOUSE}
do
    echo "${DAT}ge${RAT}/" >> ${DPATH}dataPath.dat
    all=`ls ${DAT}ge${RAT}${SPK}`
    check=`ls ${DAT}ge${RAT} | grep "contacts"`
    if [ ${check} ] # gets min and max time where the mouse was in touch with the geometric objects
    then
	min=`awk '{printf "%d\n", $1}' ${DAT}ge${RAT}/ge${RAT}_contacts.txt | sort -n | head -1`
	max=`awk -F, '{printf "%d\n", $2}' ${DAT}ge${RAT}/ge${RAT}_contacts.txt | sort -n | tail -1`
    fi
	
    # Hippocampus.
    hp=`echo "${all}" | grep "HP" | sed "s/HP_//" | tr -d '.spk' | sed '/i/d'`

    ca1=`echo "${hp}" | grep "CA1" | head -1 | cut -c1-3`
    # Hippocampus Cornu Ammonis 1.
    if [ "$ca1" == "CA1" ]
    then
	# Gets the names of the observed neurons.
	ca1=`echo "${hp}" | grep "CA1" | sed "s/CA1_//"`
	echo "${ca1}" > ${DAT}ge${RAT}/rotulosHPCA1.dat
	count=`echo "${ca1}" | wc -l` # count the number of observed neurons
	inf=0.0
	sup=100000.0
	for ca1rot in ${ca1}
	do
	    # Finds the (maximum) start time.
	    aux=`head -1 ${DAT}ge${RAT}${SPK}HP_CA1_${ca1rot}.spk`
	    if [ `echo "${aux} > ${inf}" | bc` == 1 ]
	    then
		inf=`echo ${aux}`
	    fi
	    # Finds the (minimum) end time.
	    aux=`tail -1 ${DAT}ge${RAT}${SPK}HP_CA1_${ca1rot}.spk`
	    if [ `echo "${aux} < ${sup}" | bc` == 1 ]
	    then
		sup=`echo ${aux}`
	    fi
	done
	if [ ${check} ]
	then
	    echo "${RAT}" "${count}" "${inf}" "${min}" >> ${DPATH}HPCA1mousesP1.dat
	    echo "${RAT}" "${count}" "${max}" "${sup}" >> ${DPATH}HPCA1mousesP3.dat
	    for ca1rot in ${ca1}
	    do
		echo "${ca1rot} ${DAT}ge${RAT}${SPK}HP_CA1_${ca1rot}.spk" >> ${DPATH}HPCA1mousesP1.dat
		echo "${ca1rot} ${DAT}ge${RAT}${SPK}HP_CA1_${ca1rot}.spk" >> ${DPATH}HPCA1mousesP3.dat
	    done
	else # there isn't information about when the mouse got in touch with geometric objects
	    echo "${RAT}" "${count}" "${inf}" "${sup}" >> ${DPATH}HPCA1mousesP1.dat
	    for ca1rot in ${ca1}
	    do
		echo "${ca1rot} ${DAT}ge${RAT}${SPK}HP_CA1_${ca1rot}.spk" >> ${DPATH}HPCA1mousesP1.dat
	    done
	fi

	dg=`echo "${hp}" | grep "DG" | head -1 | cut -c1-2`
        # Hippocampus Dentate Gyrus.
	if [ "$dg" == "DG" ]
	then
            # Gets the names of the observed neurons.
	    dg=`echo "${hp}" | grep "DG" | sed "s/DG_//"`
	    echo "${dg}" > ${DAT}ge${RAT}/rotulosHPDG.dat
	    count=`echo "${dg}" | wc -l` # count the number of observed neurons
	    inf=0.0
	    sup=100000.0
	    for dgrot in ${dg}
	    do
                # Finds the (maximum) start time.
		aux=`head -1 ${DAT}ge${RAT}${SPK}HP_DG_${dgrot}.spk`
		if [ `echo "${aux} > ${inf}" | bc` == 1 ]
		then
		    inf=`echo ${aux}`
		fi
	        # Finds the (minimum) end time.
		aux=`tail -1 ${DAT}ge${RAT}${SPK}HP_DG_${dgrot}.spk`
		if [ `echo "${aux} < ${sup}" | bc` == 1 ]
		then
		    sup=`echo ${aux}`
		fi
	    done
	    if [ ${check} ]
	    then
		echo "${RAT}" "${count}" "${inf}" "${min}" >> ${DPATH}HPDGmousesP1.dat
		echo "${RAT}" "${count}" "${max}" "${sup}" >> ${DPATH}HPDGmousesP3.dat
		for dgrot in ${dg}
		do
		    echo "${dgrot} ${DAT}ge${RAT}${SPK}HP_DG_${dgrot}.spk" >> ${DPATH}HPDGmousesP1.dat
		    echo "${dgrot} ${DAT}ge${RAT}${SPK}HP_DG_${dgrot}.spk" >> ${DPATH}HPDGmousesP3.dat
		done
	    else # there isn't information about when the mouse got in touch with geometric objects
		echo "${RAT}" "${count}" "${inf}" "${sup}" >> ${DPATH}HPDGmousesP1.dat
		for dgrot in ${dg}
		do
		    echo "${dgrot} ${DAT}ge${RAT}${SPK}HP_DG_${dgrot}.spk" >> ${DPATH}HPDGmousesP1.dat
		done
	    fi
	fi
    else
	# Gets the names of the observed neurons.
	echo "${hp}" > ${DAT}ge${RAT}/rotulosHP.dat
	count=`echo "${hp}" | wc -l` # count the number of observed neurons
	inf=0.0
	sup=100000.0
	for hprot in ${hp}
	do
	    # Finds the (maximum) start time.
	    aux=`head -1 ${DAT}ge${RAT}${SPK}HP_${hprot}.spk`
	    if [ `echo "${aux} > ${inf}" | bc` == 1 ]
	    then
		inf=`echo ${aux}`
	    fi
	    # Finds the (minimum) end time.
	    aux=`tail -1 ${DAT}ge${RAT}${SPK}HP_${hprot}.spk`
	    if [ `echo "${aux} < ${sup}" | bc` == 1 ]
	    then
		sup=`echo ${aux}`
	    fi
	done
	if [ ${check} ]
	then
	    echo "${RAT}" "${count}" "${inf}" "${min}" >> ${DPATH}HPmousesP1.dat
	    echo "${RAT}" "${count}" "${max}" "${sup}" >> ${DPATH}HPmousesP3.dat
	    for hprot in ${hp}
	    do
		echo "${hprot} ${DAT}ge${RAT}${SPK}HP_${hprot}.spk" >> ${DPATH}HPmousesP1.dat
		echo "${hprot} ${DAT}ge${RAT}${SPK}HP_${hprot}.spk" >> ${DPATH}HPmousesP3.dat
	    done
	else # there isn't information about when the mouse got in touch with geometric objects
	    echo "${RAT}" "${count}" "${inf}" "${sup}" >> ${DPATH}HPmousesP1.dat
	    for hprot in ${hp}
	    do
		echo "${hprot} ${DAT}ge${RAT}${SPK}HP_${hprot}.spk" >> ${DPATH}HPmousesP1.dat
	    done
	fi
    fi

    # Primary Somatic Sensory Cortex.
    # Gets the names of the observed neurons.
    s1=`echo "${all}" | grep "S1" | sed "s/S1_//" | tr -d '.spk' | sed '/i/d'`
    echo "${s1}" > ${DAT}ge${RAT}/rotulosS1.dat
    count=`echo "${s1}" | wc -l` # count the number of observed neurons
    inf=0.0
    sup=100000.0
    for s1rot in ${s1}
    do
	# Finds the (maximum) start time.
	aux=`head -1 ${DAT}ge${RAT}${SPK}S1_${s1rot}.spk`
	if [ `echo "${aux} > ${inf}" | bc` == 1 ]
	then
	    inf=`echo ${aux}`
	fi
	# Finds the (minimum) end time.
	aux=`tail -1 ${DAT}ge${RAT}${SPK}S1_${s1rot}.spk`
	if [ `echo "${aux} < ${sup}" | bc` == 1 ]
	then
	    sup=`echo ${aux}`
	fi
    done
    if [ ${check} ]
    then
	echo "${RAT}" "${count}" "${inf}" "${min}" >> ${DPATH}S1mousesP1.dat
	echo "${RAT}" "${count}" "${max}" "${sup}" >> ${DPATH}S1mousesP3.dat
	for s1rot in ${s1}
	do
	    echo "${s1rot} ${DAT}ge${RAT}${SPK}S1_${s1rot}.spk" >> ${DPATH}S1mousesP1.dat
	    echo "${s1rot} ${DAT}ge${RAT}${SPK}S1_${s1rot}.spk" >> ${DPATH}S1mousesP3.dat
	done
    else # there isn't information about when the mouse got in touch with geometric objects
	echo "${RAT}" "${count}" "${inf}" "${sup}" >> ${DPATH}S1mousesP1.dat
	for s1rot in ${s1}
	do
	    echo "${s1rot} ${DAT}ge${RAT}${SPK}S1_${s1rot}.spk" >> ${DPATH}S1mousesP1.dat
	done
    fi

    # Primary Visual Cortex.
    # Gets the names of the observed neurons.
    v1=`echo "${all}" | grep "V1" | sed "s/V1_//" | tr -d '.spk' | sed '/i/d'`
    echo "${v1}" > ${DAT}ge${RAT}/rotulosV1.dat
    count=`echo "${v1}" | wc -l` # count the number of observed neurons
    inf=0.0
    sup=100000.0
    for v1rot in ${v1}
    do
	# Finds the (maximum) start time.
	aux=`head -1 ${DAT}ge${RAT}${SPK}V1_${v1rot}.spk`
	if [ `echo "${aux} > ${inf}" | bc` == 1 ]
	then
	    inf=`echo ${aux}`
	fi
	# Finds the (minimum) end time.
	aux=`tail -1 ${DAT}ge${RAT}${SPK}V1_${v1rot}.spk`
	if [ `echo "${aux} < ${sup}" | bc` == 1 ]
	then
	    sup=`echo ${aux}`
	fi
    done
    if [ ${check} ]
    then
	echo "${RAT}" "${count}" "${inf}" "${min}" >> ${DPATH}V1mousesP1.dat
	echo "${RAT}" "${count}" "${max}" "${sup}" >> ${DPATH}V1mousesP3.dat
	for v1rot in ${v1}
	do
	    echo "${v1rot} ${DAT}ge${RAT}${SPK}V1_${v1rot}.spk" >> ${DPATH}V1mousesP1.dat
	    echo "${v1rot} ${DAT}ge${RAT}${SPK}V1_${v1rot}.spk" >> ${DPATH}V1mousesP3.dat
	done
    else # there isn't information about when the mouse got in touch with geometric objects
	echo "${RAT}" "${count}" "${inf}" "${sup}" >> ${DPATH}V1mousesP1.dat
	for v1rot in ${v1}
	do
	    echo "${v1rot} ${DAT}ge${RAT}${SPK}V1_${v1rot}.spk" >> ${DPATH}V1mousesP1.dat
	done
    fi

done

