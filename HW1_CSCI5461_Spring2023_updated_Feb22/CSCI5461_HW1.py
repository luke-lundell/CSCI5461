import csv
from bokeh.plotting import figure, show
from bokeh.layouts import row
import statistics
import math
import numpy as np
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
import matplotlib.pyplot as plt

fp_count = open("HW1-GSE62944-count.csv", "r")
fp_clinic = open("HW1-GSE62944-clinical.csv", "r")
fp_DESeq = open("HW1-DESeq2.csv","r")

fp_wilcoxon = open("HW1-WilcoxonData.csv", "w")


count = fp_count.readlines()
clinic = fp_clinic.readlines()
DESeq = fp_DESeq.readlines()

#Global Variables

gl_read_matrix = []
for line in count:
    inner_lst = line.strip().split(",")
    gl_read_matrix.append(inner_lst)

gl_genes = []
for index in range(1, len(count)):
    gl_genes.append(gl_read_matrix[index][0])

gl_patients = gl_read_matrix[0][1:]

gl_pvals = []
for index in range(1, len(DESeq)):
    inner_lst = DESeq[index].strip().split(",")
    gl_pvals.append(float(inner_lst[-1]))



#Classes
class qtObj:

    def __init__(self, gene, logval) -> None:
        self.qtTuple = (gene, logval)
        self.gene = gene
        self.val = logval

    def getName(self):
        return self.qtTuple[0]
    def getLogVal(self):
        return self.qtTuple[1]

    def setLogVal(self, val):
        self.qtTuple = (self.gene, val)

    def __str__(self):
        return f'Name: {self.qtTuple[0]}; Value: {self.qtTuple[1]}'
    def __repr__(self):
        return f'Name: {self.qtTuple[0]}; Value: {self.qtTuple[1]}'

#problem 1

#number of genes
num_genes = len(count)-1
# print(num_genes)

#number of patient samples
num_samples = len(clinic)-1
# print(num_samples)

#find out how many short and long patients there are
def find_short_long(clinic):
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
    return (short, long)


# print(f"short: {short}")
# print(f"long: {long}")



#get the total reads of every patient from a csv format
def getTotalsCSV(count):
    pat_read_tots = [0]*num_samples
    for i in range(1, len(count)): #skip first line
        gene = count[i].split(",")
        for j in range(1, len(gene)):
            pat_read_tots[j-1] += int(gene[j])   
    return pat_read_tots   

#get the total reads from a matrix in list format
def getTotalsLST(matrix):
    pat_read_tots = [0]*num_samples
    for i in range(0, len(matrix)): #skip first line
        # gene = count[i].split(",")
        for j in range(0, len(matrix[0])):
            pat_read_tots[j-1] += int(matrix[i][j])   
    return pat_read_tots                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    

#get all the patient identifiers
# pats = count[0].strip().split(",")
# pats = pats[1:]

# print(pats)
# print(total_count)


#now for making the graph using bokah
def graphTotals(sampleTotalMatrix):
    pats = count[0].strip().split(",")
    pats = pats[1:]
    p = figure(x_range=pats, title="Patient Reads", x_axis_label='patientID', y_axis_label='Reads')
    p.vbar(x=pats, top=sampleTotalMatrix, width=0.05)
    p.xgrid.grid_line_color = None
    p.y_range.start = 0
    show(p)


# #normalization based on the median
def normalizeMatrix(read_matrix, totals):
    normalized_matrix = []
    nmed = statistics.median(totals)
    for i in range(1, len(read_matrix)):
        sub_lst = []
        for j in range(1, len(read_matrix[1])):
            scalar = nmed/totals[j-1] #this gives us the median/total for a given patient
            sub_lst.append(scalar*int(read_matrix[i][j]))
            # normalized_matrix[i-1][j-1] = read_matrix[i][j]*scalar
        normalized_matrix.append(sub_lst)
    return normalized_matrix

# test case for first value:
# print("original value: ", read_matrix[1][1])
# print("nmed: ", nmed)
# print("column total: ", pat_read_tots[0])
# print("normalized value: ", normalized_matrix[0][0])

#graph the normalized data
def graphNormalization(normalized_matrix):
    pats = count[0].strip().split(",")
    pats = pats[1:]
    p = figure(x_range=pats, title="Normalized Patient Reads", x_axis_label='patientID', y_axis_label='Reads')
    p.vbar(x=pats, top=normalized_matrix, width=0.05)
    p.xgrid.grid_line_color = None
    p.y_range.start = 0
    show(p)

#normalizes the data using a log transformation
def logNormalize(normalized_matrix):
    log_normalized = []
    for i in range(0, len(normalized_matrix)):
        sub_lst = []
        for j in range(0, len(normalized_matrix[0])):
            sub_lst.append(math.log10(normalized_matrix[i][j]+1))
            # normalized_matrix[i-1][j-1] = read_matrix[i][j]*scalar
        log_normalized.append(sub_lst)
    return log_normalized

#creates and populates bins to be used in the histrogram for the entire dataset
def arrangeDatasetHistogram(log_matrix):
    numbins = 100
    binlst = [0]*numbins
    binMultiples = [0]*numbins
    max = 0
    for i in range(len(log_matrix)):
        for j in range(len(log_matrix[0])):
            if max < log_matrix[i][j]:
                max = log_matrix[i][j]

    binWidth = max/numbins 

    for x in range(len(binMultiples)):
        binMultiples[x] = x*binWidth

    for i in range(len(log_matrix)):
        for j in range(len(log_matrix[0])):
            slot = int(log_matrix[i][j]//binWidth)
            if slot == len(binlst):
                binlst[slot-1] += 1
            else:
                binlst[slot] += 1
    return (binlst, binMultiples)

#creates and populates bins to be used in the histrogram for a single column
def arrangeColumnHistogram(logColumn):
    numbins = 100
    binlst = [0]*numbins
    binMultiples = [0]*numbins
    max = 0
    for i in range(len(logColumn)):
        if max < logColumn[i]:
            max = logColumn[i]

    binWidth = max/numbins 

    for x in range(len(binMultiples)):
        binMultiples[x] = x*binWidth

    for i in range(len(logColumn)):
        slot = int(logColumn[i]//binWidth)
        if slot == len(binlst):
            binlst[slot-1] +=1
        else:
            binlst[slot] += 1
    return (binlst, binMultiples)

#Creates histogram for the log matrix
def createHistogram(binlist, binMultiples):
    # pats = count[0].strip().split(",")
    # pats = pats[1:]
    p = figure( title="Dataset Log Data", x_axis_label='Log data', y_axis_label='Frequency')
    p.vbar(x=binMultiples, top=binlist, width=0.1)
    p.xgrid.grid_line_color = None
    p.y_range.start = 0
    show(p)

#generates a page of five graphs representing logarithmic data from the fist five patients
def firstFive(log_matrix):
    (bl1, bm1) = arrangeColumnHistogram(column(log_matrix, 0))
    sample1 = figure(title="Log sample 1", x_axis_label="expression levels", y_axis_label="frequency")
    sample1.vbar(x=bm1, top=bl1, width=0.05)
    sample1.xgrid.grid_line_color = None
    sample1.y_range.start = 0
    # show(sample1)

    (bl2, bm2) = arrangeColumnHistogram(column(log_matrix, 1))
    sample2 = figure(title="Log sample 2", x_axis_label="expression levels", y_axis_label="frequency")
    sample2.vbar(x=bm2, top=bl2, width=0.05)
    sample2.xgrid.grid_line_color = None
    sample2.y_range.start = 0
    # show(sample2)

    (bl3, bm3) = arrangeColumnHistogram(column(log_matrix, 2))
    sample3 = figure(title="Log sample 3", x_axis_label="expression levels", y_axis_label="frequency")
    sample3.vbar(x=bm3, top=bl3, width=0.05)
    sample3.xgrid.grid_line_color = None
    sample3.y_range.start = 0
    # show(sample3)

    (bl4, bm4) = arrangeColumnHistogram(column(log_matrix, 3))
    sample4 = figure(title="Log sample 4", x_axis_label="expression levels", y_axis_label="frequency")
    sample4.vbar(x=bm4, top=bl4, width=0.05)
    sample4.xgrid.grid_line_color = None
    sample4.y_range.start = 0
    # show(sample4)

    (bl5, bm5) = arrangeColumnHistogram(column(log_matrix, 5))
    sample5 = figure(title="Log sample 5", x_axis_label="expression levels", y_axis_label="frequency")
    sample5.vbar(x=bm5, top=bl5, width=0.05)
    sample5.xgrid.grid_line_color = None
    sample5.y_range.start = 0
    # show(sample5)

    show(row(sample1, sample2, sample3, sample4, sample5))

#This function changes values in the matrix into tuples, so the gene labels aren't lost
#on the numbers when we do quantile normalization
def generateQuantileTupleMatrix(log_matrix):
    quantile_tuple_matrix = []
    for lst in range(len(log_matrix)):
        inner_lst = []
        for index in range(len(log_matrix[0])):
            # inner_lst.append(qtObj(gl_read_matrix[lst+1][0], log_matrix[lst][index]))
            inner_lst.append((gl_read_matrix[lst+1][0], log_matrix[lst][index]))
        quantile_tuple_matrix.append(inner_lst)
    return quantile_tuple_matrix

#generates the quantile normalization of the dataset
def quantileNormalization(qtMatrix):
    qtNomralizedList = []
    mean = 0
    rowTotal = 0
    for i in range(len(qtMatrix[0])):
        #get the column we'll be sorting
        col = column(qtMatrix, i)
        # qtMatrix[i].sort(reverse = True)
        # sortedCol = insertionSort(col)
        sortedCol = Sort_Tuple(col)
        qtNomralizedList.append(sortedCol)
    qtObjMat = tupleToObjMat(qtNomralizedList)
    qtNomralizedList = np.transpose(qtObjMat)
    
    #qtMatrix is now presumably sorted in ascending order
    #now to calculate the mean, set that as the new values for all data points in each row
    for i in range(len(qtNomralizedList)):
        #go through each row
        rowTotal = 0
        for j in range(len(qtNomralizedList[0])):
            rowTotal+= qtNomralizedList[i][j].getLogVal()
        mean = rowTotal/len(qtNomralizedList[0])
        for j in range(len(qtNomralizedList[0])):
            qtNomralizedList[i][j].setLogVal(mean)
    
    return qtNomralizedList

def firstFiveQuantile(matrix):
    # tupMat = objToTupleMat(matrix)
    logValues = []
    for i in range(len(matrix)):
        sublist = []
        for j in range(len(matrix[0])):
            sublist.append(matrix[i][j].getLogVal())
        logValues.append(sublist)
    firstFive(logValues)

#takes too long
def insertionSort(col):
    for i in range(1, len(col)):
        j = i
        while( j>0 and col[j-1].getLogVal() < col[j].getLogVal()):
            temp = col[j-1]
            col[j-1] = col[j]
            col[j] = temp
            j = j-1
    return col

#The following function is from GeeksforGeeks
# Function to sort the list by second item of tuple
def Sort_Tuple(tup):
 
    # reverse = None (Sorts in Ascending order)
    # key is set to sort using second element of
    # sublist lambda has been used
    tup.sort(key = lambda x: x[1], reverse = True)
    return tup

def tupleToObjMat(matrix):  
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            matrix[i][j] = qtObj(matrix[i][j][0], matrix[i][j][1])
    return matrix

def objToTupleMat(matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            matrix[i][j] = (matrix[i][j].getName(), matrix[i][j].getLogVal())
    return matrix


#####Part 3 funcitons#######
def createWilcoxonList(logMatrix):
    groups = sortPatients(clinic)

    extractGeneList(logMatrix, groups)

def sortPatients(clnc):
    short_group = []
    long_group = []
    for row in clnc:
        patient = row.strip().split(",")
        id = patient[2]
        id = id[1:-1] #cut out quotes
        if id == "short":
            short_group.append(patient[1])
        elif id == "long":
            long_group.append(patient[1])
    return (short_group, long_group)

def extractGeneList(logMatrix, groups):
    wilcoxonData = []
    
    for gene in gl_genes:
        longVals = []
        shortVals = []
        for cidx in range(len(logMatrix[0])):
            for ridx in range(len(logMatrix)):
                if logMatrix[ridx][cidx].getName() == gene:
                    if gl_patients[cidx] in groups[0]:
                        shortVals.append(logMatrix[ridx][cidx].getLogVal())
                    elif gl_patients[cidx] in groups[1]:
                        longVals.append(logMatrix[ridx][cidx].getLogVal())
                    else:
                        print(f'Something went wrong!')
                    
        result = stats.ranksums(shortVals, longVals)
        print(result)
        wilcoxonData.append((gene, result.pvalue))
        fp_wilcoxon.write(f'{gene},{result.pvalue}\n')



def bonferroni():
    count = 0

    # Create a list of the adjusted p-values
    p_adjusted = multipletests(gl_pvals, alpha=0.05, method='bonferroni')
    booly = p_adjusted[0]

    #report number of significant genes
    for sig in booly:
        if sig == True:
            count += 1
    return count

def benjaminiHochberg():
    count = 0 
    p_adjusted = multipletests(gl_pvals, alpha=0.05, method='fdr_bh')
    booly = p_adjusted[0]

    #report number of significant genes
    for sig in booly:
        if sig == True:
            count += 1
    return count

def graphDESeq():
    threshold = []
    # p_adjusted = multipletests(gl_pvals, alpha=0.05, method='fdr_bh', returnsorted=True)
    # sortedBH = p_adjusted[1]
    sorted_pval = sorted(gl_pvals)

    first500 = sorted_pval[0:500]
    
    # threshold = i/n*alpha
    for i in range(len(first500)):
        threshold.append((i*0.05)/num_genes)
    
    plt.plot(first500, label="DESeq2-p-value")
    plt.plot(threshold, label="BH-threshold")
    plt.xlabel("index")
    plt.ylabel("pvalues")
    plt.title("BH threshold")
    plt.legend()

    plt.show()



#from stack overflow
def column(matrix, i):
    return [row[i] for row in matrix]

def main():
    # # shrtLng = find_short_long(clinic=clinic)
    # totalMatrix = getTotalsCSV(count)
    # # print(totalMatrix)
    # # graphTotals(totalMatrix)
    # normMatrix = normalizeMatrix(read_matrix=gl_read_matrix, totals=totalMatrix)
    # normTotals = getTotalsLST(normMatrix)
    # # graphNormalization(normTotals)

    # logMatrix = logNormalize(normMatrix)
    
    # (binList, binMultiples) = arrangeDatasetHistogram(logMatrix)
    # # createHistogram(binList, binMultiples)
    # # print(binList)
    # quantile_tuple_matrix = generateQuantileTupleMatrix(logMatrix)
    # quantile_normalized = quantileNormalization(quantile_tuple_matrix)


    # # firstFiveQuantile(quantile_normalized)

    # createWilcoxonList(quantile_normalized)

    x = bonferroni()
    y = benjaminiHochberg()
    print(f'bonferroni count: {x}\nBH count: {y}')
    graphDESeq()

    # test = [[qtObj("A", 1), qtObj("B", 2), qtObj("C", 3)], [qtObj("D", 4), qtObj("E", 5), qtObj("F", 6)], [qtObj("G", 7), qtObj("H", 8), qtObj("I", 9)]]
    # print(f'original list: {test}')
    # result = quantileNormalization(test)
    # print(f'result: {result}')

    
main()



