#!/usr/bin/env python3
import sys
import gzip

def open_maybe_gzip(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def main(vcf_path):
    csq_fields = None
    sample_names = []
    out = sys.stdout

    with open_maybe_gzip(vcf_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("##INFO=<ID=CSQ"):
                desc = line.split("Description=")[1]
                if "Format:" in desc:
                    fmt_part = desc.split("Format:")[1]
                    fmt_part = fmt_part.strip('"> ')
                    csq_fields = [x.strip() for x in fmt_part.split("|")]
            elif line.startswith("#CHROM"):
                header_cols = line.lstrip("#").split("\t")
                if len(header_cols) > 9:
                    sample_names = header_cols[9:]
                else:
                    sample_names = []
                out_header = [
                    "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                    "SAMPLE", "GT", "DP", "AD", "GQ",
                    "SYMBOL", "Transcript", "Consequence", "Impact",
                    "HGVSc", "HGVSp", "Existing_variation",
                    "gnomADg_AF", "CLIN_SIG"
                ]
                out.write("\t".join(out_header) + "\n")
                break

        if csq_fields is None:
            sys.stderr.write("ERROR: Could not find CSQ definition in VCF header.\n")
            sys.exit(1)

        def idx(name):
            return csq_fields.index(name) if name in csq_fields else None

        i_allele   = idx("Allele")
        i_cons     = idx("Consequence")
        i_impact   = idx("IMPACT")
        i_symbol   = idx("SYMBOL")
        i_feat     = idx("Feature")       
        i_hgvsc    = idx("HGVSc")
        i_hgvsp    = idx("HGVSp")
        i_exist    = idx("Existing_variation")
        i_gnomadg  = idx("gnomADg_AF")
        i_clinsig  = idx("CLIN_SIG")

        for line in f:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt, qual, flt, info, fmt = fields[:9]
            samples = fields[9:]

            info_dict = {}
            for entry in info.split(";"):
                if "=" in entry:
                    k, v = entry.split("=", 1)
                    info_dict[k] = v
                else:
                    info_dict[entry] = True

            csq_raw = info_dict.get("CSQ")
            if not csq_raw:
                continue

            csq_entries = csq_raw.split(",")
            chosen_csq = None
            for e in csq_entries:
                cols = e.split("|")
                if i_symbol is not None and len(cols) > i_symbol and cols[i_symbol]:
                    chosen_csq = cols
                    break
            if chosen_csq is None:
                chosen_csq = csq_entries[0].split("|")

            def get_csq(idx_):
                if idx_ is None or idx_ >= len(chosen_csq):
                    return ""
                return chosen_csq[idx_] if chosen_csq[idx_] is not None else ""

            csq_symbol  = get_csq(i_symbol)
            csq_trans   = get_csq(i_feat)
            csq_cons    = get_csq(i_cons)
            csq_imp     = get_csq(i_impact)
            csq_hgvsc   = get_csq(i_hgvsc)
            csq_hgvsp   = get_csq(i_hgvsp)
            csq_exist   = get_csq(i_exist)
            csq_gnomadg = get_csq(i_gnomadg)
            csq_clin    = get_csq(i_clinsig)

            gt = dp = ad = gq = ""
            if samples:
                fmt_keys = fmt.split(":")
                sample_vals = samples[0].split(":")
                fmt_map = {k: v for k, v in zip(fmt_keys, sample_vals)}

                gt = fmt_map.get("GT", "")
                dp = fmt_map.get("DP", "")
                ad = fmt_map.get("AD", "")
                gq = fmt_map.get("GQ", "")

            out_row = [
                chrom, pos, vid, ref, alt, qual, flt,
                sample_names[0] if sample_names else "",
                gt, dp, ad, gq,
                csq_symbol, csq_trans, csq_cons, csq_imp,
                csq_hgvsc, csq_hgvsp, csq_exist,
                csq_gnomadg, csq_clin
            ]
            out.write("\t".join(out_row) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write(f"Usage: {sys.argv[0]} <vep_annotated.vcf[.gz]>\n")
        sys.exit(1)
    main(sys.argv[1])
