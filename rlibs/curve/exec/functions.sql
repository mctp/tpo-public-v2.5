CREATE OR REPLACE FUNCTION get_gxp(rid int)
  RETURNS table(count integer,cpm real,rpkm real,gene_id text, gene_name text, c int) AS $$
select count,cpm,rpkm,gene_id,gxp_genes.gene_name,cnv.c from
(select unnest(count) as count,unnest(cpm) as cpm,unnest(rpkm) as rpkm, generate_subscripts(count,1) as idx from gxp where runid=rid) a
left join
gxp_genes 
on a.idx=gxp_genes.idx
left join (
  select c,unnest(genes) as gene from cnv where runid=rid
  ) as cnv
on gxp_genes.gene_id=cnv.gene;
$$ 
LANGUAGE SQL;

CREATE OR REPLACE FUNCTION get_sideload_data(nm text, rid int)
  RETURNS table(id text, gene text, value numeric) AS $$
select id,gene,value from
(select unnest(value) as value, generate_subscripts(value,1) as idx from sideload_data where id=rid and name=nm) d
left join
(select unnest(id) as id,unnest(gene) as gene, generate_subscripts(id,1) as idx from sideload_meta where name=nm) m 
on d.idx=m.idx
$$ 
LANGUAGE SQL;

CREATE OR REPLACE FUNCTION get_gene(g text)
  RETURNS table(runid bigint,tp text, alignid text, cpm real,rpkm real, patient text, cohort text) AS $$
select gxp.runid, gxp.tp, gxp.alignid,cpm,rpkm,patient,cohort from 
(select runid,tp,alignid,cpm[(select idx from gxp_genes where gene_id in (g))],rpkm[(select idx from gxp_genes where gene_id in (g))] from gxp) as gxp
left join groupings
on gxp.alignid=groupings.alignid_tr;
$$
LANGUAGE SQL;

CREATE OR REPLACE FUNCTION get_fus(rid int)
  RETURNS table(
   runid bigint, 
   locus_id_5_1 text,
   cyt_5_1 text,
   locus_id_3_1 text,
   cyt_3_1 text,
   locus_id_5_2 text,
   cyt_5_2 text,
   locus_id_3_2 text,
   cyt_3_2 text,
   gene_ids_5_1 text,
   gene_names_5_1 text,
   gene_ids_3_1 text,
   gene_names_3_1 text,
   gene_ids_5_2 text,
   gene_names_5_2 text,
   gene_ids_3_2 text,
   gene_names_3_2 text,
   type integer,
   chr_5 text,
   pos_5 integer,
   str_5 text,
   chr_3 text,
   pos_3 integer,
   str_3 text,
   qname text,
   rep_5 integer,
   rep_3 integer,
   seg_5_pos integer,
   seg_5_cig text,
   seg_3_pos integer,
   seg_3_cig text,
   ftype text,
   clp_5_l integer,
   clp_5_r integer,
   clp_3_l integer,
   clp_3_r integer,
   seg_5_as integer,
   seg_5_nm integer,
   seg_5_len integer,
   seg_3_as integer,
   seg_3_nm integer,
   seg_3_len integer,
   ovr_5 real,
   ovr_3 real,
   err_5 real,
   err_3 real,
   seq_5 text,
   seq_3 text,
   low_5 real,
   low_3 real,
   unq_jnc real,
   hq_jnc boolean,
   sum_jnc integer,
   hq_sum_jnc integer,
   unq_sum_jnc real,
   avg_err_5 real,
   avg_err_3 real,
   avg_low_5 real,
   avg_low_3 real,
   max_ovr_5 real,
   max_ovr_3 real,
   unq_ovr_5 integer,
   unq_ovr_3 integer,
   cls_5 text,
   cls_3 text,
   ss_5 text,
   ss_3 text,
   frm_5 real,
   frm_3 real,
   orf boolean,
   d2a boolean,
   mot_5 text,
   mot_3 text,
   tx_dst integer,
   tx_pos_5 integer,
   tx_pos_3 integer,
   gx_dst integer,
   bp_dst integer,
   art_5 text,
   art_3 text,
   art boolean,
   sum_bpt integer,
   dst text,
   topo text,
   l1 text,
   l2 text,
   l3 text,
   hq_bpt boolean,
   hi_bpt boolean,
   hr_bpt boolean,
   rec_5 integer,
   unq_rec_5 integer,
   rec_3 integer,
   unq_rec_3 integer,
   hc_bpt boolean,
   sl_link integer,
   bs_link integer,
   ts_link integer,
   sv_link integer,
   sl_chain integer,
   bs_chain integer,
   ts_chain integer,
   sv_chain integer,
   mm2_valid boolean,
   gmap_valid boolean,
   tot_enc integer,
   tot_jnc integer,
   ts_warn boolean,
   ctg_seq text,
   tot_sp_jnc_5 integer,
   tot_sp_jnc_3 integer,
   tot_jnc_5 integer,
   tot_jnc_3 integer,
   csqf5 text,
   csqf3 text,
   csqr5 text,
   csqr3 text,
   var_id text,
   triage text,
   pair text,
   recur integer
  ) AS $$
select fus.*,COALESCE(recur,1) from 
(select *,gene_ids_5_1||'-'||gene_ids_3_1 as pair from fusions where runid in (rid)) as fus
left join
(select pair,recur::integer from fusion_recur) as recur
on fus.pair=recur.pair;
$$ 
LANGUAGE SQL;


--Returns integrated Variants, copynumber and expression for dna,rna grouping
CREATE OR REPLACE FUNCTION get_int_vars(rid int)
  RETURNS table(
        iddb bigint,
        var_id text,
        chr text,
        pos int,
        ref text,
        alt text,
        context text,
        sbs96 text,
        gene_name text,
        gene_id text,
        transcript text,
        exon text,
        hgvsc text,
        hgvsp text,
        consequence text,
        impact text,
        sift text,
        warnings text,
        hgvs_equivalent text,
        aft real,
        afn real,
        adt integer,
        dpt integer,
        adn integer,
        dpn integer,
        adt_fwd integer,
        adt_rev integer,
        str boolean,
        str_ru text,
        str_len int,
        str_diff int,
        ecnt int,
        tlod real,
        nlod real,
        flod real,
        xfet numeric,
        sor real, 
        mqrs real,
        mqs real,
        oxog real,
        dbsnp int,
        clinid int,
        clinvar text,
        cosmic_cnt integer,
        gnomad_af real,
        kg_af real,
        pon text, --problematic
        rmsk_hit text, -- problematic
        sblr real,
        pid int,
        ccf real,
        m int,
        is_indel boolean,
        all_pon int,
        runid int,
        triage text,
        cn integer,
        k integer,
        cpm real,
        rpkm real,
        recur bigint
    ) AS $$
SELECT var.*,cn,k,cpm,rpkm,COALESCE(recur,0) from
(select * from somatic where runid = rid) as var
LEFT JOIN
(select var_id,recur from somatic_recur) as sr
on var.var_id=sr.var_id
LEFT JOIN
(select runid,c as cn,k,seg,unnest(genes) as gene from cnv where runid=rid) as cnv
on var.gene_id=cnv.gene
LEFT JOIN
(select gene_name,gene_id,cpm,rpkm from get_gxp(rid)) as gxp
on var.gene_id=gxp.gene_id
$$
LANGUAGE SQL;


CREATE OR REPLACE FUNCTION get_int_ger(rid int)
  RETURNS table(
        iddb bigint,
        var_id text,
        chr text,
        pos int,
        ref text,
        alt text,
        context text,
        sbs96 text,
        gene_name text,
        gene_id text,
        transcript text,
        exon text,
        hgvsc text,
        hgvsp text,
        consequence text,
        impact text,
        sift text,
        warnings text,
        hgvs_equivalent text,
        aft real,
        afn real,
        adt integer,
        dpt integer,
        adn integer,
        dpn integer,
        adt_fwd integer,
        adt_rev integer,
        str boolean,
        str_ru text,
        str_len int,
        str_diff int,
        ecnt int,
        tlod real,
        nlod real,
        flod real,
        xfet numeric,
        sor real, 
        mqrs real,
        mqs real,
        oxog real,
        dbsnp int,
        clinid int,
        clinvar text,
        cosmic_cnt integer,
        gnomad_af real,
        kg_af real,
        pon text, --problematic
        rmsk_hit text, -- problematic
        sblr real,
        pid int,
        ccf real,
        m int,
        is_indel boolean,
        all_pon int,
        runid int,
        triage text,
        cn integer,
        k integer,
        cpm real,
        rpkm real,
        recur bigint
    ) AS $$
SELECT var.*,cn,k,cpm,rpkm,COALESCE(recur,0) from
(select * from germline where runid = rid) as var
LEFT JOIN
(select var_id,recur from germline_recur) as sr
on var.var_id=sr.var_id
LEFT JOIN
(select runid,c as cn,k,seg,unnest(genes) as gene from cnv where runid=rid) as cnv
on var.gene_id=cnv.gene
LEFT JOIN
(select gene_name,gene_id,cpm,rpkm from get_gxp(rid)) as gxp
on var.gene_id=gxp.gene_id;
$$
LANGUAGE SQL;
