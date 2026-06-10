#/bin/env bash
# Prepare an ensemble of simulations without DA from DA files

[[ ! -s "ens_list.txt" ]] && echo "Error: ens_list.txt missing" && exit 1
[[ ! -s "antime_list.txt" ]] && echo "Error: antime_list.txt missing" && exit 1

date1=$(head -1 antime_list.txt | awk '{print $1}')
date2=$(tail -2 antime_list.txt | head -1 | awk '{print $1}')

k=0
while read lline; do
    kl=$(printf "%05d" $k)

    simname="ens_sim_$kl"
    fname=$(echo $lline | awk '{print $1}')
    rstname=$(echo $lline | awk '{print $2}')

    [[ ! -s "$fname" ]] || [[ ! -s "$rstname" ]] && echo "Some input file missing." && exit 1

    cat $fname | sed -e "s/NAMESIM/$simname/" | sed -e "s/ITANF/$date1/g" | \
	sed -e "s/ITEND/$date2/g" | sed -e "s/RESTRT/$rstname/g" | \
	sed -e "s/IDTRST/3600/" > $simname.str

    k=$((k+1))
done < ens_list.txt
