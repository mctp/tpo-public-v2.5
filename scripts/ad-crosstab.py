#!/usr/bin/env python3
from moke import *
import pysam
import numpy as np

AD_MISSING = (None,)

def process_record(inp_rec, cases, noise_lo, noise_hi, out_head):
    nalts = len(inp_rec.alts)
    ADXS = []
    for alt_idx, alt in enumerate(inp_rec.alts):
        ADX[:] = 0
        for tsample, nsample in cases:
            trec = inp_rec.samples[tsample]
            nrec = inp_rec.samples[nsample]
            tads = np.array([0 if a is None else a for a in trec["AD"]])
            nads = np.array([0 if a is None else a for a in nrec["AD"]])
            if (trec["AD"] == AD_MISSING) or (None in trec["GT"]):
                tad = 0
            else:
                tad = tads[alt_idx+1]
            if nrec["AD"] == AD_MISSING or (None in nrec["GT"]):
                nad = 0
            else:
                nad = nads[alt_idx+1]
            if tad == 0:
                tadi = 0
            elif tad <= noise_lo:
                tadi = 1
            elif tad <= noise_hi:
                tadi = 2
            else:
                tadi = 3
            if nad == 0:
                nadi = 0
            elif nad <= noise_lo:
                nadi = 1
            elif nad <= noise_hi:
                nadi = 2
            else:
                nadi = 3
            ADX[tadi,nadi] += 1
        ADXS.append(":".join(ADX.flatten().astype('str')))
    
    # prepare output record
    out_rec = out_head.new_record(contig=inp_rec.contig, start=inp_rec.start, stop=inp_rec.stop)
    out_rec.ref = inp_rec.ref
    out_rec.alts = inp_rec.alts
    out_rec.info["ADX"] = ADXS
    return(out_rec)

@task
def main(inp_fn, out_fn, noise_lo=4, noise_hi=8, chrom=None):
    """
    Tumor-vs-Normal Allele Depth Cross-Tabulation.

    - inp_fn(``str``) input VCF/BCF file.
    - out_fn(``str``) output VCF/BCF file.
    - noise_lo(``int``) lower (inclusive) noise threshold [4]
    - noise_hi(``int``) higher (inclusive) noise threeshold [8]
    - chrom(`str``) chromosome to process [optional]

    """
    ## open input file
    inp_vcf = pysam.VariantFile(inp_fn)
    ## initialize global ADX array
    global ADX
    ## [0], [1,4], [5, 8], [9, Inf]
    ADX = np.zeros((4)*(4), dtype="int").reshape((4, 4))
    ## initialize output file
    out_head = pysam.VariantHeader()
    for k,v in inp_vcf.header.contigs.items():
        out_head.contigs.add(k, v.length)
    out_head.info.add("ADX", number="A", type="String", description="Allele Depth Cross-Tabulated")
    out_vcf = pysam.VariantFile(out_fn, 'w', header=out_head)
    ## process records
    with np.errstate(invalid='ignore'):
        samples = tuple(inp_vcf.header.samples)
        cases = tuple(zip(samples[0::2], samples[1::2])) # TODO
        # iterate over each record in the file
        for inp_rec in inp_vcf.fetch(chrom):
            out_rec = process_record(inp_rec, cases, noise_lo, noise_hi, out_head)
            out_vcf.write(out_rec)
    out_vcf.close()
    inp_vcf.close()

if __name__ == "__main__":
   task()
