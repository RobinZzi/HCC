import sys
from Bio import SeqIO

if len(sys.argv)==1:
    print ("//\ntype Python getPromoter.py -h for details\n//")
    sys.exit()

inputGenome = ""
inputgff = ""
outputfile = ""
name_detail = False
n = 1000

for i in range(1,len(sys.argv)):
    if sys.argv[i] == "-fa":
        inputGenome = sys.argv[i+1]
    elif sys.argv[i] == "-g":
        inputgff = sys.argv[i+1]
    elif sys.argv[i] == "-out":
        outputfile = sys.argv[i+1]
    elif sys.argv[i] == "-n":
        n = int(sys.argv[i+1])
    elif sys.argv[i] == "-detail":
        name_detail = True
    elif sys.argv[i] == "-h":
        print ("\n// designed by zby\n")
        print ("python getPromoter.py [-fa genome.fa] [-gff cdna.gff3] [-out promoter.fa] [-n 2000]")
        print ("-fa <str> Genome File. FASTA format")
        print ("-g <str> Gene Annotation File")
        print ("-out <str> Export promoter sequence")
        print ("-n <int> Promoter sequence length")
        print ("-detail <Boolean> whether fasta file's header contains detailed information. DEFAULT:False")
        print ("-h Print this page")
        print ("//")

        sys.exit()

if inputGenome=="" or inputgff=="" or outputfile=="":
    print ("//\nthe -fa -gff -out are required parameter.\n//")
    sys.exit()

gfftype = inputgff.split(".")[-1]

print ("Reading gff file ...")
gff_file = {}
for line in open(inputgff):
    if line.startswith("#") or line.strip()=="":
        pass
    else:
        l = line.strip().split("\t")
        if l[2] == "gene":
            ID = ""
            if gfftype=="gff" or gfftype=="gff3":
                ID = l[-1].split(";")[0].split("=")[1]
            elif gfftype=="gtf":
                ID = l[-1].split(";")[0].strip().split()[1][1:-1]
            else:
                print ("Please keep the input gff file ending with gff/gff3/gtf")
                sys.exit()
            if ID not in gff_file.keys():
                start = int(l[3])-1
                end = int(l[4])
                gff_file[ID] = [l[0],start,end,l[6]]

print （"Reading genome file ..."）
seq_file = {}
for query in SeqIO.parse(inputGenome,'fasta'):
    chr = query.id
    seq = str(query.seq)
    seq_file[chr] = seq

def rev(seq):
    base_trans = {"A":"T","C":"G","T":"A","G":"C","a":"t","t":"a","c":"g","g":"c","n":"n","N":"N","R":"R"}
    rev_seq_list = [base_trans[k] for k in seq]
    rev_seq = "".join(rev_seq_list)
    return rev_seq

print ("Intercept the sequence and output ...")
output_fasta = open(outputfile,"w")
for key,value in gff_file.items():
    s = value[1]-n
    e = s+n
    if s<0:
        s = 0
    else:
        s = s
    if name_detail:
        head = ">"+key+" ["+value[0]+" "+str(s)+" "+str(e)+" "+value[3]+" "+str(n)+"bp"+"]"
    else:
        head = ">"+key

    seq = seq_file[value[0]][s:e]
    if value[3] == "+":
        seq = seq
    elif value[3] == "-":
        seq = rev(seq)
    output_fasta.write(head+"\n"+seq+"\n")
output_fasta.close()
print ("finished!")
