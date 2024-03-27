#coding UTF-8
# 读取fasta文件
def read_fasta(file):
    seq_list = {}
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                seq_id = line.strip()[1:].split('#')[0]
                seq_list[seq_id] = ''
            else:
                seq_list[seq_id] += line.strip()
    return seq_list


def read_gff(file):
    gff_list = []
    with open(file, 'r') as f:
        for line in f:
            gff_list.append([line.split('\t')[8].split(';')[2][5:].split(':')[0],line.split('\t')[3],line.split('\t')[4]])
    return gff_list


#计算序列截取
def seq_intercept(seq_list,gff_list,LTR):
    #print(seq_list.keys())
    with open(f'{LTR}.fa', 'w') as f:      
        for seq_n in range(len(gff_list)):
            #print(gff_list[seq_n][0])
            #print(int(gff_list[seq_n][1])-int(gff_list[seq_n][3]))
            #print(len(seq_list[gff_list[seq_n][0]][int(gff_list[seq_n][1])-int(gff_list[seq_n][3]):int(gff_list[seq_n][2])-int(gff_list[seq_n][3])+1]))
            print(f'chr:{gff_list[0]}\tSTART:{int(gff_list[seq_n][1])}\tEND:{int(gff_list[seq_n][2])}\tlen1:{int(gff_list[seq_n][2])-int(gff_list[seq_n][1])}\tlen2:{len(seq_list[gff_list[seq_n][0]][int(gff_list[seq_n][1])-1:int(gff_list[seq_n][2])])}')
            f.write(f'>{LTR}:{seq_n}:\n{seq_list[gff_list[seq_n][0]][int(gff_list[seq_n][1])-1:int(gff_list[seq_n][2])]}\n')





# 主程序
if __name__ == '__main__':

    cds_file = 'genome.fa'
    lLTR = 'lLTR.gff'
    rLTR = 'rLTR.gff'
    seq_dict = read_fasta(cds_file)
    lgff_list = read_gff(lLTR)
    seq_intercept(seq_dict,lgff_list,"lLTR")
    rgff_list = read_gff(rLTR)
    seq_intercept(seq_dict,rgff_list,"rLTR")
    print("END")

    
