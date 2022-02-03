import sys
import multiprocessing
import time
import run_analysis

runs = 5
spawn_delay = 0

config_filename = sys.argv[1]
output_path = sys.argv[2]

children = []
for i in range(1, runs + 1):
    output_filename = output_path + f'_{i}.csv'
    process = multiprocessing.Process(target=run_analysis.main, args=(config_filename, output_filename))
    process.start()
    children.append(process)
    time.sleep(spawn_delay)
for child in children:
    child.join()
