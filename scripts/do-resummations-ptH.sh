#!/bin/zsh

#run="echo"
run=""



for file in ../fixed-order-masses/H125-LHC8-R05*-mb465-mt1735-log.pth.fxd
#for file in ../fixed-order-masses/H125-LHC8-R05*-mb0-mtinf-log.pth.fxd
#for file in ../fixed-order-masses/H125-LHC8-R05*-mb0-mt1735-log.pth.fxd

do
    #basename=${file:t:r:s/-lin//}
    basename=${file:t:r:s/-log//}
    basename=${basename:t:r:s/-log//} 
    
    echo "handling $basename"
    for order in NLL/1 NNLL/2
    do
	for scheme in a b c
	do 
            $run ../jetvheto -in $file -xQ 0.50 -out ../benchmarks-ptH/$basename-xQ050-NNLO+${order:h}$scheme.res -scheme $scheme -order ${order:t}  -ptmax 500.01 -fixed_order 2 -observable ptB
            if [[ $basename =~ "xmur050-xmuf050" ]] ; then
		$run ../jetvheto -in $file -xQ 1.00 -out ../benchmarks-ptH/$basename-xQ100-NNLO+${order:h}$scheme.res -scheme $scheme -order ${order:t}  -ptmax 500.01 -fixed_order 2 -observable ptB 
		$run ../jetvheto -in $file -xQ 0.25 -out ../benchmarks-ptH/$basename-xQ025-NNLO+${order:h}$scheme.res -scheme $scheme -order ${order:t}  -ptmax 500.01 -fixed_order 2 -observable ptB
            fi
	done
    done 
done
