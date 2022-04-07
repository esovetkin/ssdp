#!/bin/bash

# Generate documentation from parseflags
# searches the sources for commentlines following the pattern:
# // PARSEFLAG <KEYWORD> <FUNCTION-NAME> "<arguments>"
#
# it generates a table with keywords and functions in parsedef.h"
FILE="ssdp.man"
if [ -z "$1" ]
then
	echo gen_ssdpdoc.sh needs arguments!
	exit 1
fi
srcdir="$1"
shift
echo ".TH ssdp 1 \"simple sky dome projector\"" >  $FILE
echo ".SH NAME" >>  $FILE
echo "ssdp - simple sky dome projector" >>  $FILE
echo ".SH SYNOPSIS" >>  $FILE
echo ".B ssdp" >>  $FILE
echo "[command1] [command2] [...] [-f script.ssdp] [-i]" >>  $FILE
echo ".SH DESCRIPTION" >>  $FILE
echo ".B ssdp" >>  $FILE
echo "Computes sky domes according to the Perez all weather sky model [1]," >>  $FILE
echo "and can project this sky on a tilted surface. The solar position" >>  $FILE
echo "is computed from longitude, latitude, date/time, air temperature " >> $FILE
echo "and pressure, according to the Solar Position Algorithm (SPA) [2-4]." >>  $FILE
echo "For SPA wee use the freespa package [5]." >>  $FILE
echo "In addition it can process topological data and compute a horizon" >>  $FILE
echo "and thus take into account shading. The ssdp program implements a " >>  $FILE
echo "basic syntax to compute stuff. " >>  $FILE
echo ".SS Options" >>  $FILE
echo ".TP" >>  $FILE
echo ".B [command]" >>  $FILE
echo "ssdp executes any command passed on the commandline" >>  $FILE
echo ".TP" >>  $FILE
echo ".B -f" >>  $FILE
echo "pass a script file with ssdp commands." >>  $FILE
echo ".TP" >>  $FILE
echo ".B -i" >>  $FILE
echo "start an interactive shell." >>  $FILE
echo ".P" >>  $FILE
echo "There is no requirement regarding the order of the arguments." >>  $FILE
echo "This means you can pass commands on the commandline before and " >> $FILE
echo "after executing commands from a file or interactive shell." >>  $FILE
echo ".SS \"General Functionality\"" >>  $FILE
echo "The ssdp program implements a primitive script language to compute stuff." >>  $FILE
echo "The language entails two types of variables, arrays and configuration data." >>  $FILE
echo "The commands in ssdp operate on such variables. Commands have named" >>  $FILE
echo "input and output arguments (i.e. no positional arguments) and the general syntax is:" >>  $FILE
echo ".P" >>  $FILE
echo ".B <command>" >>  $FILE
echo "arg1=.. arg2=.. ..." >>  $FILE
echo ".P" >>  $FILE
echo "The configuration data stores all data needed for the simulations. The arrays are used for" >>  $FILE
echo "general input output. For convenience simple operations on arrays are implemented," >>  $FILE
echo "however, ssdp script is only very basic and does not support tests, loops, functions,. etc." >>  $FILE

# first collect all parsing flags in one file
for s in $@
do
	echo Collecting Docu flags from $s
	awk '/BEGIN_DESCRIPTION/{flag=1;next}/END_DESCRIPTION/{flag=0}flag' "$srcdir/$s" >>  descriptions
done

sect=0
section=""
while read -r line
do 
	m=$(expr match "$line" "SECTION")
	if [ "$m" -gt 0 ]
	then
		sec=$(echo $line |awk '/SECTION/{$1="";print $0}')
		if [ "$sec" != "$section" ]
		then
			section=$sec;
			echo ".SS \"$section Commands\"" >>  $FILE
		fi
	fi
	
	m=$(expr match "$line" "PARSEFLAG")
	if [ "$m" -gt 0 ]
	then
		comm=$(echo $line |awk '/PARSEFLAG/{print $2}')
		#echo .SS $comm >>$FILE
		args=$(echo $line |awk '/PARSEFLAG/{$1=$2=$3="";print $0}'|sed 's/"//g')
		echo .B $comm >>$FILE
		echo $args >>$FILE
		echo .P >>$FILE
		sect=0
	else
		m=$(expr match "$line" "DESCRIPTION")
		if [ "$m" -gt 0 ]
		then
			descr=$(echo $line |awk '/DESCRIPTION/{$1="";print $0}'|sed 's/"//g')
			echo $descr >>$FILE
			echo .P >>$FILE
		fi
		m=$(expr match "$line" "ARGUMENT")
		if [ "$m" -gt 0 ]
		then
			if [ "$sect" -eq 0 ]
			then
				echo .I input: >>$FILE
				sect=1;
			fi			
			arg=$(echo $line |awk '/ARGUMENT/{print $2}'|sed 's/"//g')
			descr=$(echo $line |awk '/ARGUMENT/{$1=$2="";print $0}'|sed 's/"//g')
			echo .TP >>$FILE
			echo .B $arg >>$FILE
			echo $descr >>$FILE
			echo .P >>$FILE
		fi
		m=$(expr match "$line" "OUTPUT")
		if [ "$m" -gt 0 ]
		then
			if [ "$sect" -lt 2 ]
			then
				echo .I output: >>$FILE
				sect=2;
			fi			
			arg=$(echo $line |awk '/OUTPUT/{print $2}'|sed 's/"//g')
			descr=$(echo $line |awk '/OUTPUT/{$1=$2="";print $0}'|sed 's/"//g')
			echo .TP >>$FILE
			echo .B $arg >>$FILE
			echo $descr >>$FILE
			echo .P >>$FILE
		fi
	fi
done < descriptions
rm descriptions
echo ".SH AUTHOR" >>  $FILE
echo "ssdp was developed by Bart E. Pieters at the Forschungszentrum Juelich GmbH" >>  $FILE
echo ".SH BUGS" >>  $FILE
echo "probably many, let me know." >>  $FILE
echo ".SH REFERENCES" >>  $FILE
echo ".TP " >>  $FILE
echo ".B [1] " >>  $FILE
echo "R. Perez, et al. \"All-Weather Model for Sky Luminance Distribution -- Preliminary Configuration and Validation.\" Solar Energy 50.3 (1993): 235-245" >> $FILE
echo ".TP " >>  $FILE
echo ".B [2] " >>  $FILE
echo "I. Reda and A. Andreas, \"Solar position algorithm for solar radiation applications.\" Solar Energy 76.5 (2004): 577-589" >>  $FILE
echo ".TP " >>  $FILE
echo ".B [3] " >>  $FILE
echo "I. Reda and A. Andreas, \"Corrigendum to Solar position algorithm for solar radiation applications.\" Solar Energy 81.6 (2007): 838" >>  $FILE
echo ".TP " >>  $FILE
echo ".B [4] " >>  $FILE
echo "NREL SPA code: http://rredc.nrel.gov/solar/codesandalgorithms/spa/" >>  $FILE
echo ".TP " >>  $FILE
echo ".B [5] " >>  $FILE
echo "freespa: https://jugit.fz-juelich.de/pearl-project/freespa" >>  $FILE
echo ".TP " >>  $FILE
