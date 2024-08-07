import datetime
import fileinput
import os



valid_times = [datetime.datetime.strptime(os.path.basename(d.rstrip()), "diag.%Y-%m-%d_%H.%M.%S.nc") for d in fileinput.input()]
idate = valid_times[0]
for i, vdate in enumerate(valid_times):
    # forecast time in minutes.
    minutes = int((vdate-idate).total_seconds()/60)
    print(f"{i+1:04d} {minutes:05d}")

