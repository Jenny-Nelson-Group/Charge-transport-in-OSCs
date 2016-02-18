#!/bin/sh

#Automagic .com CX1 job submitter. A WIP. 
#JMF 2007-09
#Bash is a stranger in an open car...

#Get Options

NCPUS=8
MEM=11800mb   #Simon Burbidge correction - lots of nodes with 12GB physical memory, leaves no overhead for OS
QUEUE="" #default route
TIME="71:58:02" # Two minutes to midnight :^)

function USAGE()
{
 cat << EOF
Jarv's Gaussian com file runner.

USAGE: ./launch_coms.sh [-nm]

OPTIONS:
	-n number of cpus
	-m amount of memory
	-q queue
	-t time

DEFAULTS (+ inspect for formatting):
	NCPUS = ${NCPUS}
	MEM   = ${MEM}
	QUEUE = ${QUEUE}
	TIME = ${TIME}
EOF
}

while getopts ":n:m:q:t:?" Option
do
    case $Option in
        n    ) NCPUS=$OPTARG;;
        m    ) MEM=$OPTARG;;
	q    ) QUEUE=$OPTARG;;
	t    ) TIME="${OPTARG}";;
        ?    ) USAGE
               exit 0;;
        *    ) echo ""
               echo "Unimplemented option chosen."
               USAGE   # DEFAULT
    esac
done

shift $(($OPTIND - 1))
#  Decrements the argument pointer so it points to next argument.
#  $1 now references the first non option item supplied on the command line
#+ if one exists.

#exit if com doesn't exit

for file in "$@"
do

	if [ ! -f ${file} ]
		then
                echo "file $@ does not exist."
		exit 2 

	fi
done

PWD=` pwd `

for COM in $*
do
 WD=${COM%/*} #subdirectory that .com file is in
 FIL=${COM#*/} #strip filetype off
 echo $PWD $WD $FIL

 if [ "${WD}" = "${FIL}" ]
 then
 echo "WD: ${WD} equiv to File: ${FIL}. Resetting WD to null."
  WD=""
 fi
 #filth to prevent error when .com is in current directory

 cat  > ${COM%.*}.sh << EOF
#!/bin/sh
#PBS -l walltime=${TIME}
#PBS -l select=1:ncpus=${NCPUS}:mem=${MEM}

module load gaussian

cp ${PWD}/${WD}/${FIL} ./
cp ${PWD}/${WD}/${FIL%.*}.chk ./

c8609 "${FIL%.*}.chk"   #convert old checkpoints to latest (i.e. for g03 checks run before ~Dec 2009)

pbsexec g09 ${FIL}

cp *.log  ${PWD}/${WD}/ 
cp *.chk  ${PWD}/${WD}/

EOF

 #echo "CAPTURED QSUB COMMAND: "
 qsub -q "${QUEUE}" ${COM%.*}.sh
done