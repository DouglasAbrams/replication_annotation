import pyBigWig as wigutils
import pandas as pd
from vcf.parser import _Info as VcfInfo, field_counts as vcf_field_counts



def query_position(chrom, pos, bw):
    return bw.values(chrom, pos, pos+1)[0]

def annotate_vcf(in_vcf, out_vcf, bw, ann_name):

    #add annotation field to header
    in_vcf.infos[ann_name] = VcfInfo(
        ann_name, vcf_field_counts['A'], 'Float',
        'Replication Value', source="douglas", version="none")

    #annotate, read by read
    for i, record in enumerate(in_vcf):

        chrom = "chr" + record.CHROM
        pos = record.POS

        #get replication val at base from bw
        val = query_position(chrom, pos, bw)

        record.add_info(ann_name, val)
        out_vcf.write_record(record)


def annotate_breakpoints_csv(breakpoints, outfile, bigwig):

    breakpoints["position_1_replication_val"] = breakpoints.apply(lambda row:
                                                  query_position("chr" + row.chromosome_1,
                                                                 row.position_1, bigwig), axis=1)
    breakpoints["position_2_replication_val"] = breakpoints.apply(lambda row:
                                                  query_position("chr" + row.chromosome_2,
                                                                 row.position_2, bigwig), axis=1)

    breakpoints.to_csv(outfile, sep="\t", index=False)

def test():
    bw = wigutils.open("wgEncodeUwRepliSeqBg02esG1bPctSignalRep1.bigWig")
    # reader = vcf.Reader(filename="Sample_T1_filtered_consensus_calls.csv")
    # reader.infos['REPLICATION'] = VcfInfo(
    #     'REPLICATION', vcf_field_counts['A'], 'Float',
    #     'Replication Value', source="douglas", version="none")
    # writer = vcf.Writer(open('annotated.vcf', 'w'), reader,
    #                     lineterminator='\n')
    breakpoints = pd.read_csv("Sample_T1_filtered_consensus_calls.csv")
    annotate_breakpoints_csv(breakpoints, "Sample_T1_filtered_consensus_calls_rep_annotated.tsv", bw)
    breakpoints = pd.read_csv("Sample_T2-A_filtered_consensus_calls.csv")
    annotate_breakpoints_csv(breakpoints, "Sample_T2-A_filtered_consensus_calls_rep_annotated.tsv", bw)
    breakpoints = pd.read_csv("Sample_T2-E_filtered_consensus_calls.csv")
    annotate_breakpoints_csv(breakpoints, "Sample_T2-E_filtered_consensus_calls_rep_annotated.tsv", bw)
    # print(breakpoints.columns)
    # annotate_vcf(breakpoints, bw, "REPLICATION")

def main():
    test()




main()




