DROP TABLE IF EXISTS groupings CASCADE;
CREATE TABLE groupings (
        id serial PRIMARY KEY,
        patient text,
        cohort text,
        alignid_t text,
        alignid_n text,
        alignid_tr text,
        alignid_nr text,
        varid text,
        dna_cap text,
        time timestamp DEFAULT NOW(),
        UNIQUE (alignid_t, alignid_n, alignid_tr, alignid_nr)
        );

DROP TABLE IF EXISTS meta CASCADE;
CREATE TABLE meta (
      id int references groupings (id),
      assembly text,
      dna_cap text,
      sex text,
      purity real,
      ploidy real,
      user_meta jsonb,
      time timestamp DEFAULT NOW()
);

DROP TABLE IF EXISTS dna_qc CASCADE;
CREATE TABLE dna_qc (
        id serial PRIMARY KEY,
        alignid text, --reference?
        genotype text,
        total_seq integer,
        aln_pct real,
        junc integer
);

DROP TABLE IF EXISTS dna_qc CASCADE;
CREATE TABLE dna_qc (
        id serial PRIMARY KEY,
        alignid text, -- reference?
        mean_target_coverage real,
        pct_selected_bases real,
        pct_target_bases_100x real,
        pf_reads int,
        percent_duplication real,
        mean_read_length real,
        pct_pf_reads_aligned real,
        strand_balance real,
        at_dropout real,
        gc_dropout real,
        genotype text,
        virus text
);

-- Data Tables
DROP TABLE IF EXISTS somatic CASCADE;
CREATE TABLE somatic (
        iddb bigserial PRIMARY KEY,
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
        runid int REFERENCES groupings (id),
        triage text
        );
CREATE INDEX runid_somatic_idx ON somatic using hash (runid);
CREATE INDEX id_somatic_idx ON somatic (var_id);
CREATE INDEX symbol_somatic_idx ON somatic (gene_name);

DROP MATERIALIZED VIEW IF EXISTS somatic_recur;
CREATE MATERIALIZED VIEW somatic_recur AS
select var_id,
  count(var_id) as recur,
  (select count(distinct runid) from somatic) as total
FROM somatic 
GROUP BY var_id
HAVING count(var_id)>1;
CREATE UNIQUE INDEX id_recur ON somatic_recur (var_id);

DROP TABLE IF EXISTS germline;
CREATE TABLE germline (LIKE somatic INCLUDING ALL);
CREATE INDEX runid_germline_idx ON germline using hash (runid);
CREATE INDEX id_germline_idx ON germline (var_id);
CREATE INDEX symbol_germline_idx ON germline (gene_name);

DROP MATERIALIZED VIEW IF EXISTS germline_recur;
CREATE MATERIALIZED VIEW germline_recur AS
select var_id,
  count(var_id) as recur,
  (select count(distinct runid) from germline) as total
FROM germline
GROUP BY var_id
HAVING count(var_id)>1;
CREATE UNIQUE INDEX ger_id_recur ON germline_recur (var_id);

DROP TABLE IF EXISTS gxp CASCADE;
CREATE TABLE gxp (
        runid bigint REFERENCES groupings (id), 
        tp text,
        alignid text unique,
        count integer[],
        cpm real[],
        rpkm real[]
        );
CREATE INDEX gxp_crispid_idx ON gxp (alignid);

DROP TABLE IF EXISTS tmb;
CREATE TABLE tmb (
  runid integer REFERENCES groupings (id),
  mutations integer,
  tmb real,
  msi_score real,
  msi_call boolean,
  mask boolean
);

DROP TABLE IF EXISTS gxp_genes;
CREATE TABLE gxp_genes (
  idx serial PRIMARY KEY,
  gene_id text,
  gene_name text
  );
CREATE INDEX gxp_genes_gene_idx ON gxp_genes (gene_name);
CREATE INDEX gxp_genes_id_idx ON gxp_genes (gene_id);

--A view with the gxp percentile 
DROP MATERIALIZED VIEW IF EXISTS gxp_pctile;
CREATE MATERIALIZED VIEW gxp_pctile as
SELECT gene_id,gene_name,pctile.* FROM   
(select
  idx as fooidx,
  percentile_disc(0.2) within group (order by foo.rpkm) as pct_20,
  percentile_disc(0.3) within group (order by foo.rpkm) as pct_30,
  percentile_disc(0.4) within group (order by foo.rpkm) as pct_40,
  percentile_disc(0.5) within group (order by foo.rpkm) as pct_50,
  percentile_disc(0.75) within group (order by foo.rpkm) as pct_75,
  percentile_disc(0.90) within group (order by foo.rpkm) as pct_90,
  percentile_disc(0.95) within group (order by foo.rpkm) as pct_95,
  percentile_disc(0.99) within group (order by foo.rpkm) as pct_99,
  percentile_disc(0.999) within group (order by foo.rpkm) as pct_999
  from 
    (select generate_subscripts(rpkm,1) as idx,unnest(rpkm) as rpkm from gxp where tp='tumor') as foo --long, unnested
  group by idx) as pctile
left join gxp_genes on pctile.fooidx=gxp_genes.idx;

DROP TABLE IF EXISTS cnv;
CREATE TABLE cnv (
  runid bigint REFERENCES groupings (id), 
  tp text,
  seg integer,
  chr text,
  seg_start integer,
  seg_end integer,
  strand text,
  c integer,
  k integer,
  lr real,
  tl real,
  al real,
  d real,
  anom real,
  mse real,
  nlr integer,
  naf integer,
  sc real,
  var_id text,
  genes text[]
);

DROP TABLE IF EXISTS cnv_ger;
CREATE TABLE cnv_ger (like cnv including all);
ALTER TABLE cnv_ger ADD CONSTRAINT cnv_ger_runid_fkey FOREIGN KEY(runid) REFERENCES groupings (id);
--CREATE INDEX varid_germline_idx ON germline (varid);

DROP TABLE IF EXISTS fusions CASCADE;
CREATE TABLE fusions (
   runid bigint REFERENCES groupings (id), 
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
   triage text
);
CREATE INDEX fus_gene_5_1_idx ON fusions (gene_ids_5_1);
CREATE INDEX fus_gene_5_2_idx ON fusions (gene_ids_5_2);
CREATE INDEX fus_gene_3_1_idx ON fusions (gene_ids_3_1);
CREATE INDEX fus_gene_3_2_idx ON fusions (gene_ids_3_2);

DROP MATERIALIZED VIEW IF EXISTS fusion_recur;
CREATE MATERIALIZED VIEW fusion_recur AS
select gene_ids_5_1||'-'||gene_ids_3_1 as pair,
  count(gene_ids_5_1||'-'||gene_ids_3_1) as recur,
  (select count(*) from (select distinct runid from fusions) as tmp) as total
FROM fusions
GROUP BY gene_ids_5_1||'-'||gene_ids_3_1
HAVING count(gene_ids_5_1||'-'||gene_ids_3_1)>1;
CREATE UNIQUE INDEX fus_id_recur ON fusion_recur (pair);


DROP TABLE IF EXISTS STRUCTURAL;
CREATE TABLE STRUCTURAL (
  runid bigint REFERENCES groupings (id), 
  bnd_id text,
  chr1 text,
  chr2 text,
  pos1 integer,
  pos2 integer,
  topo text,
  insert text,
  qual real,
  aft real,
  adt real,
  dpt real,
  manta_filter text,
  manta_score integer,
  manta_contig text,
  gene_name_1 text,
  gene_id_1 text,
  transcript1 text,
  exon1 text,
  intron1 text,
  strand1 text,
  end1 text,
  gene_name_2 text,
  gene_id_2 text,
  transcript2 text,
  exon2 text,
  intron2 text,
  strand2 text,
  end2 text,
  csq1 text,
  csq2 text,
  b1 text,
  b2 text,
  b12 text,
  exons1 text,
  exons2 text,
  var_id text,
  triage text
);

DROP TABLE IF EXISTS sideload_meta;
CREATE TABLE sideload_meta (
  name text unique,
  id text[],
  gene text[]
);

DROP TABLE IF EXISTS sideload_data;
CREATE TABLE sideload_data (
  name text REFERENCES sideload_meta (name),
  id bigint REFERENCES groupings (id),
  value numeric[],
  UNIQUE (name,id)
);

DROP TABLE IF EXISTS last_refresh;
CREATE TABLE last_refresh (
  somatic_recur timestamp,
  germline_recur timestamp,
  fusion_recur timestamp,
  gxp_pctile timestamp
);

--this needs a value for automation to start
insert into last_refresh values 
  (
  TO_TIMESTAMP('2021/01/01/12:00:00', 'YYYY/MM/DD/HH24:MI:ss'),
  TO_TIMESTAMP('2021/01/01/12:00:00', 'YYYY/MM/DD/HH24:MI:ss'),
  TO_TIMESTAMP('2021/01/01/12:00:00', 'YYYY/MM/DD/HH24:MI:ss'),
  TO_TIMESTAMP('2021/01/01/12:00:00', 'YYYY/MM/DD/HH24:MI:ss')
  );
