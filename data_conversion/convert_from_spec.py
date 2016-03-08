import numpy as np

"""
Takes a file created from a SPEC experiment and converts to a csv file with 3 columns
1: GISAXS image number
2: Time (Epoch)
3: Samx value
This file directly works with gisaxspy methods

input_file: the spec file
output_file: the file you wish to save as
start: first gisaxs image
end: last gisax image
samx_search: the beginning of the line containing the samx position (usually "#P1 ")
samx_col: the column containing the samx value
epoch_col: the column containing the epoch value
"""

def convert_from_spec(input_file, output_file, start, end, samx_search="#P1 ", samx_col=4, epoch_col=1):
    samx = []
    time = []

    record = False
    file = open(input_file)
    for line in file:
        if '#S ' + str(start) + ' ' in line:
            record = True

        if record:
            if samx_search in line:
                samx_value = line.split(' ')[samx_col]
                samx.append(float(samx_value))

            if line[0] != '#' and line != '\n':
                time_value = line.split(' ')[epoch_col]
                time.append(float(time_value))

        if '#S ' + str(end+1) + ' ' in line:
            record = False

    file.close()

    G = list(range(start, end+1)) # GISAXS image number
    data = np.asarray([G, time, samx], dtype=float)


    np.savetxt(output_file, data.transpose(), delimiter=',', fmt='%.2f')

"""
Example usage:
convert_from_spec('specfile', 'output.csv', 32, 547)

The default values correspond to the 2015 experiments
For older experiments change samx_col:
convert_from_spec('specfile', 'output.csv', 32, 547, samx_col=3)
"""