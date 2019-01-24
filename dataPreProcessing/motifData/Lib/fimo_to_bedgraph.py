import sys
import csv
import os

def main():
    infile = sys.argv[1]
    outdir = sys.argv[2]
    
    try:
        os.mkdir(outdir)
    except OSError:
        pass
    
    fh = open(infile,'r')
    reader = csv.reader(fh,delimiter="\t")
    reader.next() # skip header
    
    cur_tf = ""
    ofh = None
    merged = []
    cur_chr = ""
    
    for line in reader:
        tfid = line[0]
        chrom = line[1]
        start = int(line[2])
        end = int(line[3])
        score = float(line[5]) #log-odds score
        if score < 0:
            continue
        if tfid != cur_tf:
            if ofh is not None:
                for m in merged:
                    ofh.write("%s\t%d\t%d\t%f\n" % m)
                ofh.close()
            outfile = outdir + "/" + tfid + ".bedgraph"
            ofh = open(outfile,'w')
            cur_tf = tfid
            cur_chr = chrom
            merged = []
            print >> sys.stderr, "Merging %s" % tfid
        if chrom != cur_chr:
            for m in merged:
                ofh.write("%s\t%d\t%d\t%f\n" % m)
            merged = []
            cur_chr = chrom
        if not merged:
            merged.append((chrom,start,end,score))
        else:
            assert merged[-1][0] == chrom # make sure the chrs are the same
            last_loc = merged[-1]
            last_start = last_loc[1]
            last_end = last_loc[2]
            last_score = last_loc[3]
            while start < last_start and score > last_score:
                del merged[-1]
                last_loc = merged[-1]
                last_start = last_loc[1]
                last_end = last_loc[2]
                last_score = last_loc[3]
            if last_end > start: # there is some overlap
                #print >> sys.stderr, "last_end %d is after start %d" % (last_end,start)
                if last_score == score: # straight merge
                    #print >> sys.stderr, "scores %f are equal, straight merge" % (score)
                    #print >> sys.stderr, "modifying last one to %s, %d, %d, %f" % (chrom,last_start,end,score)
                    merged[-1] = (chrom,last_start,end,score)
                elif last_score > score: # trim off the beginning of the cur
                    #print >> sys.stderr, "last score %f gt score %f" % (last_score,score)
                    if last_end < end:
                        #print >> sys.stderr, "last end %d lt end %d, so move start to after last_end" % (last_end,end)
                        start = last_end + 1
                        #print >> sys.stderr, "adding %s, %d, %d, %f" % (chrom,start,end,score)
                        merged.append((chrom,start,end,score))
                elif last_score < score: # trim off the end of the last
                    #print >> sys.stderr, "last score %f lt score %f" % (last_score,score)
                    #print >> sys.stderr, "modifying last one to %s, %d, %d, %f" % (chrom,last_start,start-1,score)
                    merged[-1] = (chrom,last_start,start-1,last_score)
                    #print >> sys.stderr, "adding %s, %d, %d, %f" % (chrom,start,end,score)
                    merged.append((chrom,start,end,score))
            else:
                merged.append((chrom,start,end,score))
        #print >> sys.stderr, merged
    if ofh is not None:
        for m in merged:
            ofh.write("%s\t%d\t%d\t%f\n" % m)
        ofh.close()
            

if __name__ == "__main__":
    main()
