library(testthat)

#set v <- 'path/to/vault.rds'

test_that("Vault Structure",{
  vault <- readRDS(v)
  expect_type(vault, "list")
  expect_named(vault, c('maps','tables','meta','qc'))
})


test_that("dtVault Creation",{ 
  #load the vault
  vault <- readRDS(v)
  genes <- '/data/tporoot/tpo/refs/grch38/ensembl/grch38.108.clean.gtf'

  #make the dtVault
  dtanno <- dtVaultAnno(genes)
  dtv <- dtVaultUpdate(vault, dtanno)
  expect_type(dtv,"list")
  expect_named(dtv, c('maps','tables','meta','anno','storage','format'))

  #filter the dtVault
  dtvt <- vaultTriage(dtv)
  dtvf <- vaultFilter(dtvt)
  expect_type(dtvf,"list")
  expect_named(dtvf, c('maps','tables','meta','anno','storage','format'))

  #Filters should remove rows
  expect_lt(nrow(dtvf$tables$somatic),nrow(dtv$tables$somatic))
  expect_lt(nrow(dtvf$tables$germline),nrow(dtv$tables$germline))
  expect_lt(nrow(dtvf$tables$fusion),nrow(dtv$tables$fusion))

})

