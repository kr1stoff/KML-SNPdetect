import sys


rs_id_file = sys.argv[1]
vcf_file = sys.argv[2]
out_file = sys.argv[3]

# "A175loci-hg19.freq-anno.tsv"
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


def get_alt_ao_af(alt, ao, chrom, pos, dp):
    rs_alt = rs_dict[(chrom, pos)][-1]
    # 处理多等位基因情况
    if ',' in alt:
        alt_list = alt.split(",")
        ao_list = ao.split(",")
        ao = ao_list[alt_list.index(rs_alt)] if rs_alt in alt_list else "0"
        alt = rs_alt
    # 处理其他情况
    elif alt == "." or alt != rs_alt:
        alt = rs_alt
        ao = "0"
    af = format(int(ao) / int(dp), ".6f")
    return alt, ao, af


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
            dp = info_dict["DP"]
            ao = info_dict["AO"] if "AO" in info_dict else "0"
            alt, ao, af = get_alt_ao_af(alt, ao, chrom, pos, dp)
            vcf_dict[(chrom, pos)] = [chrom, pos, ref, alt, dp, af]


with open(out_file, "w") as f:
    f.write("#RS_ID\tCHROM\tPOS\tREF\tALT\tDP\tAF\n")
    for k in rs_dict:
        rs_id = rs_dict[k]
        f.write('\t'.join([rs_dict[k][0]] + vcf_dict[k]) + "\n")
