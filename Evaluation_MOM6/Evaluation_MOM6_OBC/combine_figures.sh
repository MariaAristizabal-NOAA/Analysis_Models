#!/bin/sh

##################### User input ########################
#prefix1=mom6_sst_f
#prefix2=mom6_sss_f
#prefix_comb=mom6_sst_sss.f

prefix1=mom6_ssh_f
prefix2=mom6_ssv_f
prefix_comb=mom6_ssh_ssv.f

fhours_to_combine="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48"

# direction=- is for figures to be combined vertically 
# direction=+ is for figures to be combined horizontally 
direction=-

##################### End user input ########################

#for fhour in $fhours_to_combine;
for fhour in $(seq 1 120);
do
    #echo "$fhour"
    #echo "${#fhour}"
    echo $fhour
    echo ${#fhour}
    if [ ${#fhour} -eq 1 ];
    then
        echo "yes1"
        convert ${prefix1}${fhour}.png ${prefix2}${fhour}.png ${direction}append ${prefix_comb}00${fhour}.png
    else
        if [ ${#fhour} -eq 2 ];
        then
           echo "yes2"
           convert ${prefix1}${fhour}.png ${prefix2}${fhour}.png ${direction}append ${prefix_comb}0${fhour}.png
        else
            echo "yes3"
            convert ${prefix1}${fhour}.png ${prefix2}${fhour}.png ${direction}append ${prefix_comb}${fhour}.png
        fi
    fi
done

# Make gif video
echo "Combining png figures into gif video"
#convert -delay 60 -loop 0 ${prefix_comb}*.png ${prefix_comb}.gif
convert -delay 30 -loop 0 ${prefix_comb}*.png ${prefix_comb}.gif

