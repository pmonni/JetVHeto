#!/bin/zsh

#run="echo"
run=""


#for file in fixed-order/H125-LHC7-R05*-lin.fxd fixed-order/DY-LHC7-R05*-lin.fxd
#for file in fixed-order/H125-LHC7-R05*-log.fxd fixed-order/DY-LHC7-R05*-log.fxd
#for file in fixed-order/H125-LHC7-R15*-log.fxd
#for file in fixed-order/H125-LHC7-R10*-log.fxd 
#for file in fixed-order/H125-LHC7-R04*-log.fxd 
#for file in fixed-order/H125-LHC8-R05*-log.fxd 
#for file in fixed-order/H125-LHC8-R04*-log.fxd 
#for file in fixed-order/H125-LHC8-R10*-log.fxd
for file in fixed-order/DY-LHC8-R04*-log.fxd 
#for file in fixed-order/DY-LHC8-R05*-log.fxd 
#for file in fixed-order/DY-LHC8-R10*-log.fxd 
do
    #basename=${file:t:r:s/-lin//}
    basename=${file:t:r:s/-log//}
    basename=${basename:t:r:s/-log//}
    
    echo "handling $basename"
    for order in NLL/1 NNLL/2
    do
	for scheme in a b c
	do 
            $run ./jetvheto -in $file -xQ 0.50 -out benchmarks/$basename-xQ050-NLO+${order:h}$scheme.res -scheme $scheme -order ${order:t} -ptmax 500.01 -fixed_order 1
            if [[ $basename =~ "xmur050-xmuf050" ]] ; then
		$run ./jetvheto -in $file -xQ 1.00 -out benchmarks/$basename-xQ100-NLO+${order:h}$scheme.res -scheme $scheme -order ${order:t} -ptmax 500.01 -fixed_order 1
		$run ./jetvheto -in $file -xQ 0.25 -out benchmarks/$basename-xQ025-NLO+${order:h}$scheme.res -scheme $scheme -order ${order:t} -ptmax 500.01 -fixed_order 1
            fi
        # if [[ $basename =~ "xmur100-xmuf100" ]] ; then
        #     $run ./jetvheto -in $file -xQ 1.00 -out fake-benchmarks/$basename-xQ100-NNLO+${order:h}$scheme.res -scheme $scheme -order ${order:t} -ptmax 500.01
        # fi
        # if [[ $basename =~ "xmur025-xmuf025" ]] ; then
        #     $run ./jetvheto -in $file -xQ 0.25 -out fake-benchmarks/$basename-xQ025-NNLO+${order:h}$scheme.res -scheme $scheme -order ${order:t} -ptmax 500.01
        # fi
	done
    done 
done
