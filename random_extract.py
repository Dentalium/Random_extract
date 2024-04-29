import random
import sys
import getopt
from Bio import SeqIO

inputfile = ""
threshold = 0
outputfile = "output.fa"

try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:t:o:",)
except getopt.GetoptError:
    print(sys.argv[0], "-i <inputfile> -t <threshold>")
    sys.exit(2)
for opt, arg in opts:
    if opt == "-h":
        print("Usage:", sys.argv[0], "-i <inputfile> -t <threshold> [-o outputfile]")
        sys.exit()
    elif opt == "-i":
        inputfile = arg
    elif opt == "-t":
        threshold = int(arg)
    elif opt == "-o":
        outputfile = arg

if threshold == 0:
    print("-t threshold must be determined!")
    sys.exit(2)


def random_sequence_extraction(fasta_file, threshold, outputfile):
    
    print("Loading input...")
    total_length = 0
    mysample = []
    sequences = list(SeqIO.parse(fasta_file, "fasta"))

    print("Calculating...")
    while total_length < threshold:
        # 随机抽样
        index = random.randrange(len(sequences))
        selected_seq = sequences[index]
        mysample.append(selected_seq)
        # 计算累积长度
        total_length += len(selected_seq)
        # 无放回抽样
        del sequences[index]

    # 选择更接近阈值的样本，若同等接近则随机选择
    flag = total_length - threshold - (threshold - total_length + len(selected_seq))
    if flag > 0 or (flag == 0 and random.random() < 0.5):
        mysample.pop()
        total_length -= len(selected_seq)
    
    print("Calculation finished!")
    # 输出至文件
    with open(outputfile, 'w') as file:
        # 重定向标准输出到文件
        sys.stdout = file
        for seq in mysample:
            print(">", seq.id, sep="")
            print(str(seq.seq))
    # 关闭文件后，恢复标准输出
    sys.stdout = sys.__stdout__
    file.close()
        
    print("Actural total length is", total_length)

# input fasta path and threshold
random_sequence_extraction(inputfile, threshold, outputfile)
