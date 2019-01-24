import sys
import csv


def main():
    samfile = sys.argv[1]
    outprefix = sys.argv[2]
    
    outfile1 = outprefix + "_lt_100.sam"
    outfile2 = outprefix + "_100_to_200.sam"
    outfile3 = outprefix + "_gt_200.sam"
    
    out1 = open(outfile1, 'w')
    out2 = open(outfile2, 'w')
    out3 = open(outfile3, 'w')
    
    fh = open(samfile,'r')
    samreader = csv.reader(fh,delimiter='\t')
    for line in samreader:
        template_len = abs(int(line[8]))
        if template_len < 100:
            out1.write("\t".join(line)+"\n")
        elif template_len < 200:
            out2.write("\t".join(line)+"\n")
        else:
            out3.write("\t".join(line)+"\n")
    
    out1.close()
    out2.close()
    out3.close()
    fh.close()

if __name__ == "__main__":
    main()
