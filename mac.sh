#!/bin/bash

clear

# Prints the header.
echo -e "/**  *****************************************************  **/"
echo -e "/**       ** Finding the Most Representative Graph **       **/"
echo -e "/**       **    Model for Neuronal Interactions    **       **/"
echo -e "/**       **     via Markov Chain Monte Carlo      **       **/"
echo -e "/**                                                         **/"
echo -e "/**   Author: Pedro Brandimarte Mendonca                    **/"
echo -e "/**                                                         **/"
echo -e "/**  *****************************************************  **/"
echo -e "/**   For info about the code check the 'README' file       **/"
echo -e "/**  *****************************************************  **/"

# Checks if data path directory exists and is readable.
check=path
echo -e "\nChecking data path.\n"
if [ -r ${check} ]
then
    check=`ls path/ | wc -l`

    # Checks if the directory contains 12 files.
    if [ "${check}" != "12" ]
    then
	echo -n "Organizing the data... "
	cd scripts
	./rotula.sh
	cd ../
	echo -e "done\n"
    fi
else
    echo -n "Organizing the data... "
    mkdir path
    cd scripts
    ./rotula.sh
    cd ../
    echo -e "done\n"
fi

# Computes 70% of the total free available in bits.
# Note: if the memory is given in Gb, change "tr -d 'M'" by "tr -d 'G'"
FREEmem=`top -l 1 | grep "PhysMem:" | awk '{print $10}' | tr -d 'M'`
MCmem=`echo "8 * ${FREEmem}" | bc`
echo -e "Available memory: ${MCmem} bits\n"

verif=0
echo -e "Choose the kind of computation:"
while [ ${verif} -eq 0 ]
do
    echo -n "(type 0 for penalty analysis or 1 for best graph) "
    read option
    if [ "${option}" == "0" ]
    then
        # Compiles the program.
	cd src
	echo -e ""
	while [ ${verif} -eq 0 ]
	do
	    echo -n "Do you want to clean old compilations? (y or n) "
	    read option
	    if [ "${option}" == "y" ]
	    then
		verif=1
		echo -e "\nCleaning old compilations:\n"
		make clean
		make graphPenalty
		check=`echo $?`
		if [ ${check} == 2 ]
		then
		    exit -1
		fi
		echo -e "...done\n"
	    elif [ "${option}" == "n" ]
	    then
		verif=1
		echo -e "\nCompiling changes:\n"
		make graphPenalty
		check=`echo $?`
		if [ ${check} == 2 ]
		then
		    exit -1
		fi
		echo -e "...done\n"
	    else
		echo -e "\nWrong option!\n"
	    fi
	done
	cd ../

	# Choosing the mouse
	verif=0
	while [ ${verif} -eq 0 ]
	do
	    echo -n "Choose the mouse you want to study (4, 5, 6, 9, 12, 13 or all): "
	    read opt
	    if [ "${opt}" == "4" ]
	    then
		mouse="4"
		stepsMax=`echo "3 * ${MCmem} / 756" | bc`
		verif=1
	    elif [ "${opt}" == "5" ]
	    then
		mouse="5"
		stepsMax=`echo "3 * ${MCmem} / 520" | bc`
		verif=1
	    elif [ "${opt}" == "6" ]
	    then
		mouse="6"
		stepsMax=`echo "3 * ${MCmem} / 756" | bc`
		verif=1
	    elif [ "${opt}" == "9" ]
	    then
		mouse="9"
		stepsMax=`echo "3 * ${MCmem} / 1722" | bc`
		verif=1
	    elif [ "${opt}" == "12" ]
	    then
		mouse="12"
		stepsMax=`echo "3 * ${MCmem} / 812" | bc`
		verif=1
	    elif [ "${opt}" == "13" ]
	    then
		mouse="13"
		stepsMax=`echo "3 * ${MCmem} / 930" | bc`
		verif=1
	    elif [ "${opt}" == "all" ]
	    then
		mouse="all"
		stepsMax=`echo "3 * ${MCmem} / 1722" | bc`
		verif=1
	    else
		echo -e "\nWrong option!\n"
	    fi
	done

	# Choosing the brain region
	verif=0
	echo -e ""
	while [ ${verif} -eq 0 ]
	do
	    echo -e "Choose the brain region you want to study"
	    echo -e "1. Hippocampus (HP)"
	    echo -e "2. Primary Somatic Sensory Cortex (S1)"
	    echo -e "3. Primary Visual Cortex (V1)"
	    echo -e "4. Hippocampus Dentate Gyrus (HPDG)"
	    echo -e "5. Hippocampus Cornu Ammonis 1 (HPCA1)"
	    echo -e "6. all regions\n"
	    echo -n "Type the option (1, 2, 3, 4, 5 or 6): "
	    read opt
	    if [ "${opt}" == "1" ]
	    then
		region="HP"
		verif=1
	    elif [ "${opt}" == "2" ]
	    then
		region="S1"
		verif=1
	    elif [ "${opt}" == "3" ]
	    then
		region="V1"
		verif=1
	    elif [ "${opt}" == "4" ]
	    then
		region="HPDG"
		verif=1
	    elif [ "${opt}" == "5" ]
	    then
		region="HPCA1"
		verif=1
	    elif [ "${opt}" == "6" ]
	    then
		region="0"
		verif=1
	    else
		echo -e "\nWrong option!\n"
	    fi
	done

        # Ask if the user wants an arbitrary number of MC steps.
	verif=0
	echo -e ""
	while [ ${verif} -eq 0 ]
	do
	    echo -n "Do you want to choose a fixed number of MC steps? (y or n) "
	    read option
	    if [ "${option}" == "y" ]
	    then
		while [ ${verif} -eq 0 ]
		do
		    echo -e ""
		    echo -n "Type the number of MC steps (less than ${stepsMax}): "
		    read mcSteps
		    if [ ${mcSteps} -lt ${stepsMax} ]
		    then
			MCmem=${mcSteps}
			stepsMax=1
			verif=1
		    else
			echo -e "\nWrong option!"
		    fi
		done
	    elif [ "${option}" == "n" ]
	    then
		stepsMax=0
		verif=1
	    else
		echo -e "\nWrong option!\n"
	    fi
	done

        # Checks if output directory exists and is readable.
	check=out/outPenalty
	if [ -r ${check} ]
	then
	    echo -e ""
	    echo -n "Removing old results "
	    if [ "${region}" == "0" ] && [ "${mouse}" == "all" ]
	    then
		rm out/outPenalty/*
		echo -e "...done\n"
	    elif [ "${region}" == "0" ]
	    then
		rm out/outPenalty/*${mouse}*
		echo -e "...done\n"
	    elif [ "${mouse}" == "all" ]
	    then
		rm out/outPenalty/*${region}*
		echo -e "...done\n"
	    else
		rm out/outPenalty/*${mouse}${region}*
		echo -e "...done\n"
	    fi
	else
	    echo -e ""
	    echo -n "Creating directory for output results... "
	    mkdir out/outPenalty
	    echo -e "done\n"
	fi

        # Runs the program.
	cd bin
	echo -e "\nRunning the program...\n"
	echo -n "Start of run: "
	date
	ini=$(date +%s%N) # initial time with nanoseconds accuracy
	if [ "${region}" == "0" ] && [ "${mouse}" == "all" ]
	then
	    for region in "HP" "S1" "V1" "HPDG" "HPCA1"
	    do
		for mouse in "4" "5" "6" "9" "12" "13"
		do
		    ./graphPenalty ../path/ ../out/outPenalty/ ${MCmem} ${stepsMax} ${region} ${mouse}
		done
	    done
	elif [ "${region}" == "0" ]
	then
	    for region in "HP" "S1" "V1" "HPDG" "HPCA1"
	    do
		./graphPenalty ../path/ ../out/outPenalty/ ${MCmem} ${stepsMax} ${region} ${mouse}
	    done
	elif [ "${mouse}" == "all" ]
	then
	    for mouse in "4" "5" "6" "9" "12" "13"
	    do
		./graphPenalty ../path/ ../out/outPenalty/ ${MCmem} ${stepsMax} ${region} ${mouse}
	    done
	else
	    ./graphPenalty ../path/ ../out/outPenalty/ ${MCmem} ${stepsMax} ${region} ${mouse}
	fi
	fim=$(date +%s%N) # final time with nanoseconds accuracy
	echo -e ""
	echo -n "End of run: "
	date
	tempo=`echo "scale = 10; (${fim} - ${ini}) / 60000000000" | bc`
	echo -e "\nRun time: ${tempo} min\n"
	cd ../

    elif [ "${option}" == "1" ]
    then
        # Compiles the program.
	cd src
	echo -e ""
	while [ ${verif} -eq 0 ]
	do
	    echo -n "Do you want to clean old compilations? (y or n) "
	    read option
	    if [ "${option}" == "y" ]
	    then
		verif=1
		echo -e "\nCleaning old compilations:\n"
		make clean
		make bestGraph
		check=`echo $?`
		if [ ${check} == 2 ]
		then
		    exit -1
		fi
		echo -e "...done\n"
	    elif [ "${option}" == "n" ]
	    then
		verif=1
		echo -e "\nCompiling changes:\n"
		make bestGraph
		check=`echo $?`
		if [ ${check} == 2 ]
		then
		    exit -1
		fi
		echo -e "...done\n"
	    else
		echo -e "\nWrong option!\n"
	    fi
	done
	cd ../

	# Choosing the mouse
	verif=0
	while [ ${verif} -eq 0 ]
	do
	    echo -n "Choose the mouse you want to study (4, 5, 6, 9, 12 or 13): "
	    read opt
	    if [ "${opt}" == "4" ]
	    then
		mouse="4"
		stepsMax=`echo "3 * ${MCmem} / 756" | bc`
		verif=1
	    elif [ "${opt}" == "5" ]
	    then
		mouse="5"
		stepsMax=`echo "3 * ${MCmem} / 520" | bc`
		verif=1
	    elif [ "${opt}" == "6" ]
	    then
		mouse="6"
		stepsMax=`echo "3 * ${MCmem} / 756" | bc`
		verif=1
	    elif [ "${opt}" == "9" ]
	    then
		mouse="9"
		stepsMax=`echo "3 * ${MCmem} / 1722" | bc`
		verif=1
	    elif [ "${opt}" == "12" ]
	    then
		mouse="12"
		stepsMax=`echo "3 * ${MCmem} / 812" | bc`
		verif=1
	    elif [ "${opt}" == "13" ]
	    then
		mouse="13"
		stepsMax=`echo "3 * ${MCmem} / 930" | bc`
		verif=1
	    else
		echo -e "\nWrong option!\n"
	    fi
	done

	# Choosing the brain region
	verif=0
	echo -e ""
	while [ ${verif} -eq 0 ]
	do
	    echo -e "Choose the brain region you want to study"
	    echo -e "1. Hippocampus (HP)"
	    echo -e "2. Primary Somatic Sensory Cortex (S1)"
	    echo -e "3. Primary Visual Cortex (V1)"
	    echo -e "4. Hippocampus Dentate Gyrus (HPDG)"
	    echo -e "5. Hippocampus Cornu Ammonis 1 (HPCA1)\n"
	    echo -n "Type the option (1, 2, 3, 4 or 5): "
	    read opt
	    if [ "${opt}" == "1" ]
	    then
		region="HP"
		verif=1
	    elif [ "${opt}" == "2" ]
	    then
		region="S1"
		verif=1
	    elif [ "${opt}" == "3" ]
	    then
		region="V1"
		verif=1
	    elif [ "${opt}" == "4" ]
	    then
		region="HPDG"
		verif=1
	    elif [ "${opt}" == "5" ]
	    then
		region="HPCA1"
		verif=1
	    else
		echo -e "\nWrong option!\n"
	    fi
	done

        # Ask if the user wants an arbitrary number of MC steps.
	verif=0
	echo -e ""
	while [ ${verif} -eq 0 ]
	do
	    echo -n "Do you want to choose a fixed number of MC steps? (y or n) "
	    read option
	    if [ "${option}" == "y" ]
	    then
		while [ ${verif} -eq 0 ]
		do
		    echo -e ""
		    echo -n "Type the number of MC steps (less than ${stepsMax}): "
		    read mcSteps
		    if [ ${mcSteps} -lt ${stepsMax} ]
		    then
			MCmem=${mcSteps}
			stepsMax=1
			verif=1
		    else
			echo -e "\nWrong option!"
		    fi
		done
	    elif [ "${option}" == "n" ]
	    then
		stepsMax=0
		verif=1
	    else
		echo -e "\nWrong option!\n"
	    fi
	done

        # Checks if output directory exists and is readable.
	check=out/outBestGraph
	if [ -r ${check} ]
	then
	    echo -e ""
	    echo -n "Removing old results "
	    rm out/outBestGraph/*${mouse}${region}*
	    echo -e "...done\n"
	else
	    echo -e ""
	    echo -n "Creating directory for output results... "
	    mkdir out/outBestGraph
	    echo -e "done\n"
	fi

	# Choosing the method for probability computation
	verif=0
	while [ ${verif} -eq 0 ]
	do
	    echo -n "Choose the method for probability computation (1, 2 or 3): "
	    read method
	    if [ "${method}" == "1" ]
	    then
		verif=1
	    elif [ "${method}" == "2" ]
	    then
		verif=1
	    elif [ "${method}" == "3" ]
	    then
		verif=1
	    else
		echo -e "\nWrong option!\n"
	    fi
	done

	# Choosing the penalty constant
	echo -e ""
	echo -n "Type the positive penalty constant: "
	read penal

        # Runs the program.
	cd bin
	echo -e "\nRunning the program...\n"
	echo -n "Start of run: "
	date
	ini=$(date +%s%N) # initial time with nanoseconds accuracy
	./bestGraph ../path/ ../out/outBestGraph/ ${MCmem} ${stepsMax} ${region} ${mouse} ${method} ${penal}
	fim=$(date +%s%N) # final time with nanoseconds accuracy
	echo -e "\n"
	echo -n "End of run: "
	date
	tempo=`echo "scale = 10; (${fim} - ${ini}) / 60000000000" | bc`
	echo -e "\nRun time: ${tempo} min\n"
	cd ../
    else
	echo -e "\nWrong option!\n"
    fi
done

