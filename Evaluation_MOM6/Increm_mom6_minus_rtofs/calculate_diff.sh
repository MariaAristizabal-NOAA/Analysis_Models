#!/bin/sh

hycom_tools=/work/noaa/hwrf/save/maristiz/HYCOM-tools

rtofs_restart1=/work/noaa/hwrf/noscrub/maristiz/RTOFS/arch_restart/2020072812/restart_out.a 
rtofs_restart2=/work/noaa/hwrf/noscrub/maristiz/RTOFS/arch_restart/2020072806/restart_out.a 

mom6_restart=/work/noaa/hwrf/noscrub/maristiz/MOM6_analysis_files/hat10_mom6nc2archv_restart_2020_210_12_April19.a

idm=1135
jdm=633
kdm=41
itlrec=247
#itlrec=1
increc=1
numrec=1

# Extract a field (hycom_extract.F)
#hycom_extract - Usage:  hycom_extract fin.a idm jdm kdm itlrec increc numrec fout.a
#                          hycom_extract fin.a idm jdm rec.txt fout.a
#
#                 Outputs the input fields itlrec+(n-1)*increc+(k-1),
#                 for k=1:kdm and n=1:numrec (or n=1:e-o-f if numrec=0).
#                 Or outputs the input fields listed in rec.txt,
#                 which can be in any order (note no kdm).

#  Note that itlrec and increc are w.r.t. idm*jdm sized records, but
#   numrec is w.r.t. idm*jdm*kdm.

#  fin*.a is assumed to contain idm*jdm 32-bit IEEE real values for
#   each array, in standard f77 element order, followed by padding
#   to a multiple of 4096 32-bit words, but otherwise with no control
#   bytes/words, and input values of 2.0**100 indicating a data void.

#${hycom_tools}/bin/hycom_extract ${rtofs_restart} ${idm} ${jdm} ${kdm} ${itlrec} ${increc} ${numrec} rtofs_restart_temp.a 
#${hycom_tools}/bin/hycom_extract ${rtofs_restart2} ${idm} ${jdm} ${kdm} ${itlrec} ${increc} ${numrec} rtofs_restart_temp2.a 

#${hycom_tools}/bin/hycom_extract ${mom6_restart} ${idm} ${jdm} ${kdm} ${itlrec} ${increc} ${numrec} mom6_restart_temp.a 

# Get the differences (hycom_diff_print.F)
#  hycom_diff_print - Usage:  hycom_diff_print f1.a f2.a idm jdm k cfmt [if il jf jl]
#                 prints a list of all values in the (if:il,jf:jl)
#                 sub-array of the k-th (1:idm,1:jdm) arrays in f1.a
#                 and f2.a that differ.
#                 cfmt     - output format of form "(2i5,...)"
#                 if,jf    - first point in sub-array
#                 il,jl    - last  point in sub-array
#                 can have if>il and/or jf>jl for a reversed print order

k=1
${hycom_tools}/bin/hycom_diff_print_mf rtofs_restart_temp.a rtofs_restart_temp2.a ${idm} ${jdm} ${k} "(2i5,E14.7)" > diff_rtofs_vs_rtofs.txt  

#${hycom_tools}/bin/hycom_diff_print_mf rtofs_restart_temp.a mom6_restart_temp.a ${idm} ${jdm} ${k} "(2i5,f10.4,f10.4,f10.4)" > diff_rtofs_vs_mom6.txt  
#${hycom_tools}/bin/hycom_diff_print_mf rtofs_restart_temp.a rtofs_restart_temp2.a ${idm} ${jdm} ${k} "(2i5,f10.4,f10.4,f10.4)" > diff_rtofs_vs_rtofs.txt  

#${hycom_tools}/bin/hycom_diff_print_mf rtofs_restart_temp.a mom6_restart_temp.a ${idm} ${jdm} ${k} "(2i5,f10.4,f10.4,E14.7)" > diff_rtofs_vs_mom6.txt  

#${hycom_tools}/bin/hycom_diff_print_mf rtofs_restart_temp.a rtofs_restart_temp2.a ${idm} ${jdm} ${k} "(2i5,E14.7,E14.7,E14.7)" > diff_rtofs_vs_rtofs.txt  
#${hycom_tools}/bin/hycom_diff_print_mf rtofs_restart_temp.a rtofs_restart_temp2.a ${idm} ${jdm} ${k} "(2i5,f10.4,f10.4,f10.4)" > diff_rtofs_vs_rtofs.txt  
#${hycom_tools}/bin/hycom_diff_print_mf rtofs_restart_temp.a mom6_restart_temp.a ${idm} ${jdm} ${k} "(2i5,f10.4,f10.4,E14.7)" 1 1 1 10 > diff_rtofs_vs_mom6.txt  
#${hycom_tools}/bin/hycom_diff_print_mf rtofs_restart_temp.a mom6_restart_temp.a ${idm} ${jdm} ${k} "(2i5,f10.4,f10.4,E14.7)" > diff_rtofs_vs_mom6.txt  
#${hycom_tools}/bin/hycom_diff_print_mf rtofs_restart_temp.a rtofs_restart_temp2.a ${idm} ${jdm} ${k} "(2i5,f10.4,f10.4,E14.7)" 1 1 1 10 > diff_rtofs_vs_rtofs.txt  
#${hycom_tools}/bin/hycom_diff_print rtofs_restart_temp.a rtofs_restart_temp.a ${idm} ${jdm} ${k} "(2i5,f10.4)" > diff.txt  

# Get a HYCOM a of an ascii file (ascii2hycom.F).
#/bin/ascii2hycom

#  ascii2hycom - Usage:  ascii2hycom  f.txt idm jdm [spval] [i1 j1] fhycom.a
#                        asciif2hycom f.txt idm jdm [spval] [i1 j1] fhycom.a cfmt
#
#  Outputs a HYCOM ".a" copy of an ascii (plain text) file.

#  The input array is (1:idm,1:jdm), output is (i1:idm,j1:jdm)

#  f.txt is assumed to contain idm*jdm numerical values for for
#   each array, with data voids indicated by spval and input via:
#   read(21,*) a(1:idm,1:jdm).
#  for asciif2hycom, input format is in cfmt, e.g. '(6144I1)', via:
#   read(21,cfmt) a(1:idm,1:jdm)
#

#${hycom_tools}/bin/ascii2hycom_mf diff_rtofs_vs_rtofs.txt ${idm} ${jdm} .TRUE. 1 1 diff_temp_rtofs_minus_mom6.a "(2i5,f10.4,f10.4,E14.7)"

cp ${hycom_tools}/bin/ascii2hycom_mf . 

./ascii2hycom_mf diff_rtofs_vs_rtofs.txt ${idm} ${jdm} 0.1267651E+31 1 1 diff_temp_rtofs_minus_mom6.a  
#./ascii2hycom_mf diff_rtofs_vs_rtofs.txt ${idm} ${jdm} '.TRUE.' 1 1 diff_temp_rtofs_minus_mom6.a  
#./ascii2hycom_mf diff_rtofs_vs_rtofs2.txt ${idm} ${jdm} '.TRUE.' diff_temp_rtofs_minus_mom6.a > ascii2hycom_mf.log 
#./ascii2hycom_mf diff_rtofs_vs_mom6.txt ${idm} ${jdm} diff_temp_rtofs_minus_mom6.a  
