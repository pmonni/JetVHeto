#!/bin/zsh

#run="echo"
run=""


# the loop-mass flag is read directly from the fized-order file, therefore loop over all ptjet files
## remember to switch on the rescaling in the expansion.f90 routine
for file in fixed-order-masses/H125-LHC8-R05*-log.fxd
#for file in fixed-order-masses/H125-LHC8-R05*-mb465-mt1735-log.fxd
#for file in fixed-order-masses/H125-LHC8-R05*-mb0-mtinf-log.fxd
#for file in fixed-order-masses/H125-LHC8-R05*-mb0-mt1735-log.fxd

do
    #basename=${file:t:r:s/-lin//}
    basename=${file:t:r:s/-log//}
    basename=${basename:t:r:s/-log//} 
    
    echo "handling $basename"
    for order in NLL/1 NNLL/2
    do
	for scheme in a b c
	do 
            $run ./jetvheto -in $file -xQ 0.50 -out benchmarks-masses/$basename-xQ050-NNLO+${order:h}$scheme.res -scheme $scheme -order ${order:t}  -ptmax 500.01 -fixed_order 2
            if [[ $basename =~ "xmur050-xmuf050" ]] ; then
		$run ./jetvheto -in $file -xQ 1.00 -out benchmarks-masses/$basename-xQ100-NNLO+${order:h}$scheme.res -scheme $scheme -order ${order:t}  -ptmax 500.01 -fixed_order 2
		$run ./jetvheto -in $file -xQ 0.25 -out benchmarks-masses/$basename-xQ025-NNLO+${order:h}$scheme.res -scheme $scheme -order ${order:t}  -ptmax 500.01 -fixed_order 2
            fi
        # if [[ $basename =~ "xmur100-xmuf100" ]] ; then
        #     $run ./jetvheto -in $file -xQ 1.00 -out matched/$basename-xQ100-NNLO+${order:h}$scheme.res -scheme $scheme -order ${order:t} -ptmax 500.01
        # fi
        # if [[ $basename =~ "xmur025-xmuf025" ]] ; then
        #     $run ./jetvheto -in $file -xQ 0.25 -out matched/$basename-xQ025-NNLO+${order:h}$scheme.res -scheme $scheme -order ${order:t} -ptmax 500.01
        # fi
	done
    done 
done


