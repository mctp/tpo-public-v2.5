test_that("dtVault Creation",{ 
  vault <- readRDS(v)
  genes <- '/data/tporoot/tpo/refs/grch38/ensembl/grch38.108.clean.gtf'

  #make the dtVault
  dtanno <- dtVaultAnno(genes)
  dtv <- dtVaultUpdate(vault, dtanno)
  #filter the dtVault
  dtvt <- vaultTriage(dtv)
  dtvf <- vaultFilter(dtvt)
  #expected variants? CTNNA2, BRAF, MET
  expect_true(all(
    c(
    'chr2:80581735_C/A',
    'chr7:140800456_C/T',
    'chr7:116774880_G/T'
    ) %in% 
    dtvf$tables$somatic$var_id
    ))
  })
