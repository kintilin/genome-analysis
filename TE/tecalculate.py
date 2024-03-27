import os
import sys

def getfile(fold):
    for names in os.walk(fold):
        file_list = []
        for filename in names:
            if type(filename) == str:
                path = filename
            if type(filename) == list:
                if filename != []:
                    for i in filename:
                        if i.endswith("meg"):
                            new_filename = path  + i
                            file_list.append(new_filename)
        return file_list
def gettime(file):
    with open(file, "r") as fi:
        ltr_pair = ""
        ltr_time_dic = {}
        ltr_dic = {}
        f = fi.readlines()
        for line in f:
            lin = line.strip().split(" ")
            if line.startswith("[1]"):
                if len(lin) == 2:
                    ltr1 = lin[1].strip("#")
            if line.startswith("[2]"):
                if len(lin) == 2:
                    ltr2 = lin[1].strip("#")
                    ltr_pair = ltr1 + ";" + ltr2
                    ltr_dic[ltr1] = ltr2
                if len(lin) == 3:
                    distance = float(lin[2])
                    r2 = 1.99*10**(-9) #####这里是碱基突变率#####
                    time = (distance/(2*r2))/1000000
                    ltr_time_dic[ltr_pair] = time
                    out_line = ltr_pair + "\t" + str(time)
                    #print(ltr_pair, time)
                    #print(time, ltr_time_dic)
    return out_line  ###ltr_time_dic

def writelis(lis, fil):
    with open(fil, "w") as out_f:
        for it in lis:
            line1 = it + "\n"
            out_f.write(line1)
    out_f.close()
f_list = getfile(r"output2/")
new_list = []
for item in f_list:
    time_line = gettime(item)
    new_list.append(time_line)
writelis(new_list, "ltr_insertime_out.txt")

