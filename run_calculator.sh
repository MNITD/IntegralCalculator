#!/bin/bash
#script runs IntegralCalculator 

#options
#-c change directory with executable (default current)
#-t set number of times to run (default 1)

#default options args
directory=$PWD
runsNum=1

#other global variables
minTotal=9999999
sumTotal=0
inaccuracy=0


#check option args
#args of option should NOT be another option
checkargs () {
if [[ $OPTARG =~ ^-[c/t]$ ]]
then
echo "Option argument check failed!"
echo "Unknow argument $OPTARG for option $opt!"
exit 1
fi
}

#with ':' after option, because args for option is required
while getopts "c:t:" opt
do
case $opt in

c) checkargs 
directory=$OPTARG
echo "The directory was changed to : $OPTARG"
;;

t) checkargs
runsNum=$OPTARG
echo "The number of times to run was changed to : $OPTARG"
;;

*) 
echo "No reasonable options found!"
;;
esac
done


#other functions

abs() {
   [ $1 -lt 0 ] && echo $((-$1)) || echo $1
}

run () {
    echo $($directory/IntegralCalculator)
}

main () {
    cd $directory
    directory=$PWD
    for i in `seq 0 $[runsNum-1]`;
    do
          outputs[$i]=" $[i+1] - $(run) "
          IFS=': ' read -r -a array <<< "${outputs[$i]}"
          arrayTotal[$i]=${array[3]}
          sumTotal=$[ sumTotal+${arrayTotal[$i]}]
          if  [ ${arrayTotal[$i]} -lt $minTotal ]
          then
              minTotal=${arrayTotal[$i]}
          fi
    done
    
    average=$[sumTotal/runsNum]
    
    for i in `seq 0  $[runsNum-1]`;
    do
        inaccuracy=$(abs $[(arrayTotal[i]-average)*100/average])
        echo "${outputs[$i]} Inaccuracy: $inaccuracy%"
    done
  
    echo "Min Time: $minTotal"
}

main



















