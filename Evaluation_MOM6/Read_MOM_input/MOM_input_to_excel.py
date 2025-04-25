# User input
MOM_input_file = '/work/noaa/hurricane/save/maristiz/MOM6-NOAA-EMC-Aug9_2024/ocean_only/NHC_with_MOM_input_from_HAFS/OUTPUT/MOM_parameter_doc.all'

from pandas import DataFrame
import openpyxl

f = open(MOM_input_file,'r')
mom_input = f.read()

flag_name = []
flag_value = []
flag_explanation = []
for l,line in enumerate(mom_input.split('\n')[2:]):
    if line.split('jfdkjksjdkj')[0] == '':
        flag_name.append(' ')
        flag_value.append(' ')
        flag_explanation.append(' ')
    else:
        if line[0] == '!':
            flag_name.append(line[6:-4])
            flag_value.append(' ')
            flag_explanation.append(' ')
        else:
            if line.split('!')[0][0] == ' ':
                flag_name.append(' ')
                flag_value.append(' ')
                flag_explanation.append(line.split('!')[1])
            else:
                flag_name.append(line.split('!')[0].split('=')[0])
                if len(line.split('!'))>1:
                    flag_value.append(line.split('!')[0].split('=')[1])
                    flag_explanation.append(line.split('!')[1])
                else:
                    flag_value.append(' ')
                    flag_explanation.append(' ')
    #print(l, ' ',line, ' ', flag_name)

df = DataFrame({'Flag Name': flag_name, 'Flag Value': flag_value, 'Flag Explanation':flag_explanation})

df.to_excel('MOM_input_all_ePBL_OM4.xlsx', sheet_name='sheet1', index=False)


