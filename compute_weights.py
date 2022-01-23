import sys
import csv

POP_COL = 8
COMP_COL = 9
COUNTY_COL = 10
VRA_COL = 11

IGNORE_ROWS = 10000

target_score = 100

total_rows = 0
pop_total = 0.0
comp_total = 0.0
county_total = 0.0
vra_total = 0.0

for filename in sys.argv[1:]:
    with open(filename) as infile:
        reader = csv.reader(infile)
        for i, row in enumerate(reader):
            if i < IGNORE_ROWS:
                continue
            total_rows += 1
            pop_total += float(row[POP_COL])
            comp_total += float(row[COMP_COL])
            county_total += float(row[COUNTY_COL])
            vra_total += float(row[VRA_COL])

print('Population weight:\t', target_score * total_rows / pop_total)
print('Compactness weight:\t', target_score * total_rows / comp_total)
print('County weight:\t\t', target_score * total_rows / county_total)
print('VRA weight:\t\t', target_score * total_rows / vra_total)
