import pyBigWig as wigutils
import vcf
from vcf.parser import _Info as VcfInfo, field_counts as vcf_field_counts


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
        val = bw.values(chrom, pos, pos+1)[0]

        record.add_info(ann_name, val)
        out_vcf.write_record(record)


def test():
    bw = wigutils.open("wgEncodeUwRepliSeqBg02esG1bPctSignalRep1.bigWig")
    reader = vcf.Reader(filename="Sample_T2-E_samtools_germline_small.vcf")
    reader.infos['REPLICATION'] = VcfInfo(
        'REPLICATION', vcf_field_counts['A'], 'Float',
        'Replication Value', source="douglas", version="none")
    writer = vcf.Writer(open('annotated.vcf', 'w'), reader,
                        lineterminator='\n')

    annotate_vcf(reader, writer, bw, "REPLICATION")

def main():
    test()




main()




