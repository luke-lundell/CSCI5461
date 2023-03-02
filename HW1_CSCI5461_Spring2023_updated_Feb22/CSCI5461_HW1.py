import csv
from bokeh.plotting import figure, show
import statistics

fp_count = open("HW1-GSE62944-count.csv", "r")
fp_clinic = open("HW1-GSE62944-clinical.csv", "r")
fp_DESeq = open("HW1-DESeq2.csv","r")

count = fp_count.readlines()
clinic = fp_clinic.readlines()
DESeq = fp_DESeq.readlines()

# print(count[0])

read_matrix = []
for line in count:
    inner_lst = line.strip().split(",")
    read_matrix.append(inner_lst)

#problem 1

#number of genes
num_genes = len(count)-1
# print(num_genes)

#number of patient samples
num_samples = len(clinic)-1
# print(num_samples)

# def find_short_long():
short = 0
long = 0
for line in clinic:
    patient = line.split(",")
    # print(patient[2])
    id = patient[2]
    id = id[1:-2] #cut out quotes
    if id == "short":
        short+=1
    elif id == "long":
        long+=1


# print(f"short: {short}")
# print(f"long: {long}")



#get the total reads of every patient
pat_read_tots = [0]*num_samples
for i in range(1, len(count)): #skip first line
    gene = count[i].split(",")
    for j in range(1, len(gene)):
        pat_read_tots[j-1] += int(gene[j])                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           

#get all the patient identifiers
pats = count[0].strip().split(",")
pats = pats[1:]

# print(pats)
# print(total_count)


#now for making the graph using bokah
p = figure(x_range=pats, title="Patient Reads", x_axis_label='patientID', y_axis_label='Reads')
p.vbar(x=pats, top=pat_read_tots, width=0.05)
p.xgrid.grid_line_color = None
p.y_range.start = 0

# add a line renderer with legend and line thickness
# p.line(x, y, legend_label="Line 1", line_width=2)


# show the plot
# show(p)

# #normalization:
normalized_matrix = []
nmed = statistics.median(pat_read_tots)
for i in range(1, len(read_matrix)):
    sub_lst = []
    for j in range(len(1, read_matrix[1])):
        scalar = nmed/pat_read_tots[j-1] #this gives us the median/total for a given patient
        sub_lst.append(scalar*int(read_matrix[i][j]))
        # normalized_matrix[i-1][j-1] = read_matrix[i][j]*scalar
    normalized_matrix.append(sub_lst)



    
