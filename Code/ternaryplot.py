import ternary

### Scatter Plot
scale = 1
figure, tax = ternary.figure(scale=scale)
# Set Axis labels and Title
fontsize = 20



tax.boundary(linewidth=1.0)
tax.gridlines(multiple=0.1, color="blue")
tax.gridlines(color="blue", multiple=0.05, linewidth=0.5)

# Code for extracting data
# type = 0 for simulation data and 1 for null data
# trial_number is file number, ranging from 1 to 500
def getSingleFileData(type, trial_number, n, K, maf):
    x_0 = []
    x_1 = []
    x_2 = []

    filepath = '/home/jiahuang/Dropbox/Princeton/Sophomore Summer/Summer Project/Summer2016/New Simulation Result/file_n=' + str(n) + '-maf=' + str(int(maf * 100)) + '-K=' + str(K) + '/'
    maf_formated = maf * 100
    if type == 0:
        file_type = 'simulation'
    else:
        file_type = 'nullSimulation'
    filename = filepath + file_type + str(trial_number) + '_n=' + str(n) + '-maf=' + str(int(maf_formated)) + '-K=' + str(K) + '.txt'

    with open(filename) as f:
        for line in f:
            raw_data = map(float, line.split())
            if raw_data[0] == 0:
                x_0.append(tuple(raw_data[1:K+1]))
            if raw_data[0] == 1:
                x_1.append(tuple(raw_data[1:K+1]))
            if raw_data[0] == 2:
                x_2.append(tuple(raw_data[1:K+1])) 
    return [x_0, x_1, x_2]

def getAllFileData(type, n, K, maf):
    x_0 = []
    x_1 = []
    x_2 = []

    for trial_number in range(1, 501):
        filepath = '/home/jiahuang/Dropbox/Princeton/Sophomore Summer/Summer Project/Summer2016/New Simulation Result/file_n=' + str(n) + '-maf=' + str(int(maf * 100)) + '-K=' + str(K) + '/'
        maf_formated = maf * 100
        if type == 0:
            file_type = 'simulation'
        else:
            file_type = 'nullSimulation'
        filename = filepath + file_type + str(trial_number) + '_n=' + str(n) + '-maf=' + str(int(maf_formated)) + '-K=' + str(K) + '.txt'

        with open(filename) as f:
            for line in f:
                raw_data = map(float, line.split())
                if raw_data[0] == 0:
                    x_0.append(tuple(raw_data[1:K+1]))
                if raw_data[0] == 1:
                    x_1.append(tuple(raw_data[1:K+1]))
                if raw_data[0] == 2:
                    x_2.append(tuple(raw_data[1:K+1])) 
    return [x_0, x_1, x_2]

# type, trial number, n, K, maf
type = 1
trial_no = 2
n = 100
K = 3
maf = 0.1
final_data = getSingleFileData(type, trial_no, n, K, maf)

if type == 0:
    file_type = 'simulation'
else:
    file_type = 'nullSimulation'
#final_data = getAllFileData(1, 10, 3, 0.1)

# Plot a few different styles with a legend
points =  final_data[0]
tax.scatter(points, marker='s', color='red', label="x=0")
points =  final_data[1]
tax.scatter(points, marker='D', color='green', label="x=1")
points =  final_data[2]
tax.scatter(points, marker='D', color='blue', label="x=2")

tax.set_title(file_type + ' n=' + str(n) + ' maf=' + str(float(maf)) + ' K=' + str(K), fontsize=fontsize)
tax.left_axis_label("Isoform 1", fontsize=fontsize)
tax.right_axis_label("Isoform 2", fontsize=fontsize)
tax.bottom_axis_label("Isoform 3", fontsize=fontsize)

tax.legend()
tax.ticks(axis='lbr', linewidth=1, multiple=0.1)
tax.clear_matplotlib_ticks()

tax.show()