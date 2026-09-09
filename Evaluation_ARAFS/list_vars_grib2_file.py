import grib2io
from collections import defaultdict

#grib2_file = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/arafs-input/RRFS/rrfs.t00z.prslev.13km.f000.na.grib2'
#grib2_file = '/gpfs/f6/drsa-hurr1/world-shared/noscrub/Maria.Aristizabal/arafs-input/COMGFSv16/gfs.20260401/00/atmos/gfs.t00z.pgrb2.0p25.f000'
grib2_file = '/gpfs/f6/ar-cpu/scratch/Maria.Aristizabal/scrub/ARAFS_alaska_coupled_ocean_EXP1/2025122000/00E/atm_ic/gfs.t00z.pgrb2ab.0p25.f000_tmp'

#grb = grib2io.open(grib2file,mode='r')

# List of unique variables
with grib2io.open(grib2_file) as grib:
    # Store unique variables in a dictionary with their long descriptions
    unique_vars = {}
    for msg in grib:
        if msg.shortName not in unique_vars:
            unique_vars[msg.shortName] = msg.fullName

print("Unique Variables in File:")
print("-" * 50)
for var, desc in sorted(unique_vars.items()):
    print(f"{var:<12} : {desc}")


# Print all levels
# Use defaultdict to group levels by variable short name
var_levels = defaultdict(list)
var_names = {}

# Open file and scan messages
with grib2io.open(grib2_file) as grib:
    for msg in grib:
        short_name = msg.shortName
        level_str = str(msg.level)

        # Store variable full description if not already captured
        if short_name not in var_names:
            var_names[short_name] = msg.fullName

        # Append level if it hasn't been added yet for this variable
        if level_str not in var_levels[short_name]:
            var_levels[short_name].append(level_str)

# Print formatted output
print(f"{'Variable':<12} | {'Description':<35} | Levels")
print("-" * 100)

for short_name, levels in sorted(var_levels.items()):
    description = var_names[short_name]
    levels_formatted = ", ".join(levels)

    # Truncate long descriptions or level strings for clean output
    if len(levels_formatted) > 45:
        levels_formatted = levels_formatted[:42] + f"... ({len(levels)} levels)"

    print(f"{short_name:<12} | {description[:35]:<35} | {levels_formatted}")

