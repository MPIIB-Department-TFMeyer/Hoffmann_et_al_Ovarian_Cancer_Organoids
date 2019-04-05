#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May  5 10:23:32 2017

@author: hilmar
"""

import pysam
import re

# order taken from http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf
annotype_sort_order = """chromosome_number_variation
exon_loss_variant
frameshift_variant
stop_gained
stop_lost
start_lost
splice_acceptor_variant
splice_donor_variant
rare_amino_acid_variant
missense_variant
inframe_insertion
disruptive_inframe_insertion
inframe_deletion
disruptive_inframe_deletion
5_prime_UTR_truncation+exon_loss_variant
3_prime_UTR_truncation+exon_loss
splice_branch_variant
splice_region_variant
splice_branch_variant
stop_retained_variant
initiator_codon_variant
synonymous_variant
initiator_codon_variant+non_canonical_start_codon
stop_retained_variant
coding_sequence_variant
5_prime_UTR_variant
3_prime_UTR_variant
5_prime_UTR_premature_start_codon_gain_variant
upstream_gene_variant
downstream_gene_variant
TF_binding_site_variant
regulatory_region_variant
miRNA
custom
sequence_feature
conserved_intron_variant
intron_variant
intragenic_variant
conserved_intergenic_variant
intergenic_region
coding_sequence_variant
non_coding_exon_variant
nc_transcript_variant
gene_variant
chromosome"""

annotype_order_hash = {}
tmp = annotype_sort_order.strip().split("\n")
for ii in tmp:
    annotype_order_hash[ii]=tmp.index(ii)
annotype_order_hash_max = max(annotype_order_hash.values())+1

impact_order = """HIGH
MODERATE
LOW
MODIFIER
"""

impact_order_hash = {}
tmp = impact_order.strip().split("\n")
for ii in tmp:
    impact_order_hash[ii]=tmp.index(ii)


def import_vcf(ifile):
    """
    read VCF and extract relevant fields. 
    Extract highest impact changes for each variant site
    """

    vcf_reader = pysam.VariantFile(filename=ifile)
    
      
#    header_info_tags = [e for e in vcf_reader.header.info]
    
    header_hash = {}
    for e in vcf_reader.header.records:
        items = e.items()
        if len(items)<1:
            continue
        
        tmp = {}
        for f in e.iteritems():
            tmp[f[0]]=f[1]
                
        eid = tmp["ID"] if "ID" in tmp else None
        
        if eid is None:
            print tmp
            sys.stop("No ID")
        
        header_hash[(e.type, eid)] = tmp
    
    # Get all annotation fields from ANN provided by SnpEff
    ann_entry = header_hash[('INFO','ANN')]["Description"].replace('"','')
    af_list =  ann_entry.replace("Functional annotations:","").strip(" '").split(" | ")
    af_ind_hash = {}
    for a in af_list:
        af_ind_hash[a] = af_list.index(a)

    rec_cnt = 0
    
    # global overview on relevant variants per gene
    per_gene_variants = {'LOF': {}, 'NMD': {}, 'ANN': {}}

    unknown_annotation_types = {}
    
    all_sites = {}
    
    for record in vcf_reader.fetch():
        
#        if rec_cnt < 5:
#            print record

        chrom = record.contig
        pos = record.pos
        #rid = record.ID
        ref = record.ref
        alt = record.alts
        qual = record.qual

        # Replace missing FILTER value with empty string instead of None
        filter_stat = record.filter
        if not filter_stat is None:
            filter_str = "PASS" if len(filter_stat)==0 else ";".join(filter_stat)
        else:
            filter_str = "."
        
        polymorphism_site="."
        if not record.id is None and record.id[0:2]=="rs":
            polymorphism_site = record.id
        
        
        try:
            dp = record.info["DP"]
        except:
            dp = -1
        
        all_sites[(chrom, pos, ref)] = {'basic':{'QUAL':qual, 'FILTER':filter_str, 'SNP':polymorphism_site, 'DP':dp}}
        
        
        # collect all called alleles
        all_alleles = [ref]
        all_alleles.extend([str(e) for e in alt])
        
        # map alleles to samples
        sample_alleles = {}
        
        for sn in record.samples:
            #s_alleles = record.samples[sn].alleles
             
            
            # exclude duplicate reference samples - Strelka specific feature
            if sn[1]==":":
                continue
            
            gt_raw = record.samples[sn]["GT"]
            if gt_raw == (None, None):
                gt_ind=(0,0)
            else:
                gt_ind = gt_raw 
                # Freebayes seems to produce variant quals with lower QUAL that have 0/0 as GT. Change this in 0/1
                # This should not be mixed up with "no calls" inserted by bctools which are coded as empty entries (.:.:. etc)
                if gt_raw == (0,0):
                    gt_ind = (0,1)

            # extract the allele index
            if gt_ind is None: # this can only happen if the GT field is read as None
                var_allele_ind = None
            else:
                var_allele_ind = [e for e in gt_ind if e > 0]
                
            
            if not var_allele_ind is None and len(var_allele_ind)>0:
                
                for vi in var_allele_ind:
                    curr_allele = all_alleles[vi]
                    if not curr_allele in sample_alleles:
                        sample_alleles[curr_allele] = {}
                    try:
                        depth = record.samples[sn]['DP'] 
                    except:
                        depth = record.samples[sn]['DP'] 
                        
                    ad = record.samples[sn]['AD']
                    if len(ad) != len(all_alleles):
                        sys.stop("Number of AD entries and alleles does not match.")
                    sample_alleles[curr_allele][sn] = (depth, ad[vi], record.samples[sn]['GQ'], -1) 
                    #sample_alleles[sn][all_alleles[vi]] = (depth, -1, -1, -1) 
        
        all_sites[(chrom, pos, ref)]['alleles_per_sample'] = sample_alleles
        
        # treat the SnpEff info parameters
        # Loss of function info
        """ From MacArthur et al Science 2012 p. 823-828        try:
            record.info["LOF"]
        except:
            has_lof = False
            We adopted a definition for LoF variants expected to correlate with complete loss of function of the affected transcripts: 
            stop codon–introducing (nonsense) or 
            splice site–disrupting single-nucleotide variants (SNVs), 
            insertion/deletion (indel) variants predicted to disrupt a transcript’s reading frame,
            or larger deletions removing either the first exon or more than 50% of the protein-coding sequence of the affected transcript. 
            """
    
        has_lof = True
        try:
            record.info["LOF"]
        except:
            has_lof = False
        if has_lof:
            lof_record = record.info['LOF']
            for rec_str in lof_record:
                gg = rec_str.strip("()").split("|", -1)
                gene_symbol = gg[0]
                lof_result = (gg[1], gg[2], gg[3]) # gene id, number of transcripts affected, % of tx affected
                if not gene_symbol in per_gene_variants['LOF']:
                    per_gene_variants['LOF'][gene_symbol]= [lof_result]
                else:
                    per_gene_variants['LOF'][gene_symbol].append(lof_result)
        
        # nonsense-mediated decay
        # for definition see Nat Rev Mol Cell Biol. 2004 Feb;5(2):89-99.
        has_nmd = True
        try:
            record.info["NMD"]
        except:
            has_nmd = False
        if has_nmd:
            nmd_record = record.info['NMD']
            for rec_str in nmd_record:
                gg = rec_str.strip("()").split("|", -1)
                gene_symbol = gg[0]
                nmd_result = (gg[1], gg[2], gg[3]) # gene id, number of transcripts affected, % of tx affected
                if not gene_symbol in per_gene_variants['NMD']:
                    per_gene_variants['NMD'][gene_symbol]= [nmd_result]
                else:
                    per_gene_variants['NMD'][gene_symbol].append(nmd_result)

        curr_site_anno = {}
        curr_site_anno2 = {}
        
        
        
        # complete annotation
        has_ann = True
        try:
            record.info["ANN"]
        except:
            has_ann = False
        if has_ann:
            ann_record = record.info['ANN']
            #print ann_record
            # Loop through all annotation records for this site
            # SnpEff can annotate different variant alleles for the same site
            # Here we split all annotations for different alleles and for all listed transcripts effects per allele
            # Records are comma separated and will automatically split by pyVCF
            for rec_str in ann_record:
                # Split fields within record
                ar = rec_str.split("|", -1)
                allele = ar[af_ind_hash['Allele']]
                atype = ar[af_ind_hash['Annotation']]
                impact = ar[af_ind_hash['Annotation_Impact']]
                gsymbol = ar[af_ind_hash['Gene_Name']]
                
                if atype == "structural_interaction_variant":
                    continue
                
                # reorder effect type by importance
                if not atype in annotype_order_hash:
                    if not atype in unknown_annotation_types:
                         unknown_annotation_types[atype] = 1
                    else:
                         unknown_annotation_types[atype] += 1
                
                # Store records for the current site by allele, gene symbol, impact and effect type
                if not allele in curr_site_anno:
                    curr_site_anno[allele] = {}
                
                if not gsymbol in curr_site_anno[allele]:
                    curr_site_anno[allele][gsymbol]={}
                
                if not impact in curr_site_anno[allele][gsymbol]:
                    curr_site_anno[allele][gsymbol][impact] = {}
                
                if not atype in curr_site_anno[allele][gsymbol][impact]:
                    curr_site_anno[allele][gsymbol][impact][atype] = []
                
                # Store/append the full record
                curr_site_anno[allele][gsymbol][impact][atype].append(ar)
                
                #############################################
                # Store/append the full record under a shared impact/type level
                # This makes sorting effects by impact easier
                impact_id = (impact, atype)
                if not allele in curr_site_anno2:
                    curr_site_anno2[allele] = {}
                
                if not gsymbol in curr_site_anno2[allele]:
                    curr_site_anno2[allele][gsymbol]={}
                
                if not impact_id in curr_site_anno2[allele][gsymbol]:
                    curr_site_anno2[allele][gsymbol][impact_id] = []

                curr_site_anno2[allele][gsymbol][impact_id].append(ar)


        per_allele_highest_impact = {}
        # Identify the highest impact effect for each allele                
        for c_allele in curr_site_anno2.keys():
            for c_gene in curr_site_anno2[c_allele].keys():
                sorted_impacts = sorted(curr_site_anno2[c_allele][c_gene].keys(), key=lambda x: (impact_order_hash[x[0]], annotype_order_hash[x[1]] if x[1] in annotype_order_hash else annotype_order_hash_max))
                #print c_allele, c_gene, sorted_impacts

                highest_impact = sorted_impacts[0]
                gene_symbol = c_gene
        
                if not gene_symbol in per_gene_variants['ANN']:
                    per_gene_variants['ANN'][gene_symbol]= {}
                
                if not highest_impact in per_gene_variants['ANN'][gene_symbol]:
                    per_gene_variants['ANN'][gene_symbol][highest_impact] = []
                
                per_gene_variants['ANN'][gene_symbol][highest_impact].append((chrom, pos, ref, c_allele))

                if not c_allele in per_allele_highest_impact:
                    per_allele_highest_impact[c_allele] = {}
                
                per_allele_highest_impact[c_allele][c_gene] = highest_impact
                

        all_sites[(chrom, pos, ref)]['annotation'] = (curr_site_anno, per_allele_highest_impact)
        
        #print(unknown_annotation_types)

        rec_cnt += 1

        all_samples = list(vcf_reader.header.samples)

    return((all_sites, per_gene_variants, af_ind_hash, all_samples))

def export_annotations_by_gene(variant_data, ofile):
    """
    Export all variants affecting each gene together with highest impact effect 
    and allelelic frequencies or VAF for each sample.
    """

    per_gene_variants = variant_data[1]
    all_samples = variant_data[3]
    anno_ind_hash = variant_data[2]

    anno_impact_ind = anno_ind_hash["Annotation_Impact"]
    anno_etype_ind = anno_ind_hash["Annotation"]
    anno_hgsvc_ind = anno_ind_hash["HGVS.c"]
    anno_hgsvp_ind = anno_ind_hash["HGVS.p"]
    anno_tx_ind = anno_ind_hash["Feature_ID"]
    anno_bt_ind = anno_ind_hash["Transcript_BioType"]

    ofile.write("\t".join(["Gene","Impact","Effect","Chromosome","Position","RefAllele","VariantAllele","Coverage","TranscriptID","TranscriptType", "HGVS.c", "HGVS.p"]) )
    ofile.write("\t" + "\t".join(all_samples))
    ofile.write("\n")

    # loop through all genes
    for gene in sorted(per_gene_variants["ANN"].keys()):
        
        # for each gene consider all highest impact categories identified
        for impact in per_gene_variants["ANN"][gene].keys():
            impact_level, etype = impact
            
            # Identify variants for each highest impact category
            for v in per_gene_variants["ANN"][gene][impact]:
                chrom, pos, ref, alt_allele = v
                
                var_details = variant_data[0][(chrom, pos, ref)]
                
                vd_alleles_per_sample = var_details["alleles_per_sample"]
                vd_annotation_all = var_details["annotation"][0][alt_allele][gene][impact_level][etype]
                vd_annotation_highest_impact = var_details["annotation"][1][alt_allele][gene]
                
                highest_impact_effects = [e for e in vd_annotation_all if e[anno_impact_ind]==vd_annotation_highest_impact[0] and e[anno_etype_ind]==vd_annotation_highest_impact[1] ]
                h_i_hgvsc = [e[anno_hgsvc_ind] for e in highest_impact_effects]
                h_i_hgvsp = [e[anno_hgsvp_ind] for e in highest_impact_effects]
                h_i_tx = [e[anno_tx_ind] for e in highest_impact_effects]
                h_i_bt = [e[anno_bt_ind] for e in highest_impact_effects]
                
                vd_basic = var_details["basic"]
        
                ofile.write("%s\t%s\t%s\t%s\t%d\t%s\t%s\t" % (gene, impact_level, etype, chrom, pos, ref, alt_allele) )
                                
                ofile.write("%d\t%s\t%s\t%s\t%s" % (vd_basic["DP"], ", ".join(h_i_tx), ", ".join(h_i_bt), ", ".join(h_i_hgvsc), ", ".join(h_i_hgvsp)) )
                
                for aa in vd_alleles_per_sample.keys():
                    sample_affected = vd_alleles_per_sample[aa]
                    
                    for s in all_samples:
                        if s in sample_affected:
                            ii = vd_alleles_per_sample[aa][s]
                            vaf = float(ii[1])/ii[0]
                            ofile.write("\t%0.3g" % vaf)
                        else:
                            ofile.write("\t-1")
                ofile.write("\n")

if __name__ == '__main__':

    import argparse, bz2file, gzip, fileinput, sys
    import time
    parser = argparse.ArgumentParser(description='Convert SnpEff-annotated VCF to tall-skinny table format')
    parser.add_argument('input_vcf', type=str, help='Input vcf')

    parser.add_argument('-o','--output-file', action='store', type=str, default=None,
                       help='File where output should go to. Defaults to STDOUT')

    args = parser.parse_args()

    ifile = args.input_vcf

    if args.output_file is None:    
        ofile = sys.stdout
    else:
        ofile = open(args.output_file,"w")
       

    variant_data = import_vcf(ifile)

    export_annotations_by_gene(variant_data, ofile)
    ofile.close()


