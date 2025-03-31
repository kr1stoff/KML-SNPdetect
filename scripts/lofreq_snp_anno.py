import sys


rs_id_file = sys.argv[1]
depth_file = sys.argv[2]
vcf_file = sys.argv[3]
out_file = sys.argv[4]

# "A175loci-hg19.freq-anno.tsv"
# "align/8435A3/8435A3.aligned.recalibrated.bam.depth"
# "result/vcf/8435A3/8435A3.lofreq.vcf"
# "result/vcf/8435A3/8435A3.snp_anno.tsv"


# chr-pos - rs_id 字典
rs_dict = {}
with open(rs_id_file) as f:
    for line in f:
        lns = line.strip().split("\t")
        """
        {('chr1', '14155402'): ['rs7520386', 'G', 'A'],
        """
        rs_dict[(lns[0], lns[1])] = (lns[2:])
# rs_dict


# 深度字典
dp_dict = {}
with open(depth_file) as f:
    for line in f:
        lns = line.strip().split("\t")
        # depth 和 vcf 位置一致
        dp_dict[(lns[0], lns[1])] = lns[2]


vcf_dict = {}
with open(vcf_file) as f:
    for line in f:
        if line.startswith("##"):
            continue
        elif line.startswith("#CHROM"):
            headers = line.strip().split("\t")
            hdr = {hie[1]: hie[0] for hie in enumerate(headers)}
        else:
            lns = line.strip().split("\t")
            chrom = lns[hdr["#CHROM"]]
            pos = lns[hdr["POS"]]
            ref = lns[hdr["REF"]]
            alt = lns[hdr["ALT"]]
            info = lns[hdr["INFO"]]
            info_dict = {i.split("=")[0]: i.split("=")[1] for i in info.split(";")}
            info_dp = info_dict["DP"]
            info_af = info_dict["AF"]
            vcf_dict[(chrom, pos)] = [chrom, pos, ref, alt, info_dp, info_af]


with open(out_file, "w") as f:
    f.write("#RS_ID\tCHROM\tPOS\tREF\tALT\tDP\tAF\n")
    for k in rs_dict:
        if k not in vcf_dict:
            if k not in dp_dict:
                dp_dict[k] = "0"
            f.write('\t'.join([rs_dict[k][0]] + list(k) +
                    rs_dict[k][1:] + [dp_dict[k], "0"]) + "\n")
        else:
            rs_id = rs_dict[k]
            f.write('\t'.join([rs_dict[k][0]] + vcf_dict[k]) + "\n")
