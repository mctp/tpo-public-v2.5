#' Format variants for display
#'
#' @param var Variants from vault or database
#' @param dict The curve data dictionary
#' @param detail Include columns from the show_detail column rather than the show_table
#'
#' @return out$table Variants formatted for display
#' @return out$hide index of columns to hide
#' @return out$renderpct index of columns to show as percentages
#'
#' @export
format_variants <- function(vars,dict,tp,detail=FALSE){
        if(detail){
          .dict <- dict %>% filter(table==tp & show_detail)
        }else{
          .dict <- dict %>% filter(table==tp & show_table)
        }
        rnm <- as.character(.dict$display_name)
        names(rnm) <- .dict$db_name
        fmt <- vars %>% select(.dict$db_name) %>% S4Vectors::rename(rnm)

        hide=c()
        renderpct=which(.dict$db_name %in% c('aft','afn','gnomad_af','kg_af','xfet'))

        return(list(table=fmt,hide=hide,renderpct=renderpct))
}

#' Format gene expression for display
#'
#' @param gxp Gene expression from vault or DB
#' @param dict The curve data dictionary
#'
#' @return out$table Gene expression table formatted for display
#' @return out$hide index of columns to hide
#'
#' @export
format_gxp <- function(gxp,dict,detail=FALSE){
        if(detail){
          .dict <- dict %>% filter(table=='gxp' & show_detail)
        }else{
          .dict <- dict %>% filter(table=='gxp' & show_table)
        }
        rnm <- as.character(.dict$display_name)
        names(rnm) <- .dict$db_name
        fmt <- gxp %>% select(.dict$db_name) %>% S4Vectors::rename(rnm)

        hide=c()
        renderpct=c()

        return(list(table=fmt,hide=hide,renderpct=renderpct))
}

#' Format CNV Tables for display
#'
#' @param d cnv segments from the database
#' @param dict The curve data dictionary
#'
#' @return out$table cnv table for display
#' @return out$hide index of columns to hide
#'
#' @export
format_cnv <- function(d,dict,table='segments'){
        d$k[is.na(d$k)] <- -1
        d$c[is.na(d$c)] <- -1
        d$event <- 'Neutral'
        d[d$k==-1 | d$c==-1,'event'] <- 'Unknown'
        d[d$c>2,'event'] <- 'Gain'
        d[d$c<2,'event'] <- 'Loss'
        d[d$c==0,'event'] <- 'Hom. Deletion'
        d[d$c>=2 & d$k==0,'event'] <- 'LOH'
        if(table=='segments'){
          .dict <- dict %>% filter(table=='cnv' & show_table)
          rnm <- as.character(.dict$display_name)
          names(rnm) <- .dict$db_name
          fmt <- d %>% select(.dict$db_name) %>% S4Vectors::rename(rnm)
        }else if(table=='genes'){
          .dict <- dict %>% filter(table=='cnv' & show_table)
          rnm <- as.character(.dict$display_name)
          names(rnm) <- .dict$db_name
          fmt <- d %>% unnest(genes) %>% 
                    select(gene_id=genes,.dict$db_name) %>% S4Vectors::rename(rnm) %>%
                    merge(gxp_to_id, by='gene_id') %>%
                    select(Gene=gene_name,ID=gene_id,everything())
          
        }

        hide=c()
        renderpct=c()

        return(list(table=fmt,hide=hide,renderpct=renderpct))
}

#' Format SV tables for display
#'
#' @param sv structural variants from the database
#' @param dict The curve data dictionary
#'
#' @return out$table sv table for display
#' @return out$hide index of columns to hide
#'
#' @export
format_sv <- function(sv,dict,detail=FALSE){
        w <- which(is.na(sv$gene_name_1))
        sv$gene_name_1[w] <- paste0(sv$chr1[w],':',sv$pos1[w])
        w <- which(is.na(sv$gene_name_2))
        sv$gene_name_2[w] <- paste0(sv$chr2[w],':',sv$pos2[w])
        if(detail){
          .dict <- dict %>% filter(table=='struc' & show_detail)
        }else{
          .dict <- dict %>% filter(table=='struc' & show_table)
        }
          rnm <- as.character(.dict$display_name)
          names(rnm) <- .dict$db_name
          fmt <- sv %>% select(.dict$db_name) %>% S4Vectors::rename(rnm)

        hide=c()
        renderpct=which(.dict$db_name %in% c('aft'))

        return(list(table=fmt,hide=hide,renderpct=renderpct))
}

#' Format Fusion tables for display
#'
#' @param fus fusions from the database
#' @param dict The curve data dictionary
#'
#' @return out$table fus table for display
#' @return out$hide index of columns to hide
#'
#' @export
format_fus <- function(fus,dict,detail=FALSE){
        # prefer names to cytoband
        w <- which(fus$gene_names_3_1=="")
        fus$gene_names_3_1[w] <- fus$gene_names_3_2[w]
        w <- which(fus$gene_names_3_1=="")
        fus$gene_names_3_1[w] <- fus$cyt_3_1[w]
        w <- which(fus$gene_names_3_1=="")
        fus$gene_names_3_1[w] <- fus$cyt_3_2[w]
        #
        w <- which(fus$gene_names_5_1=="")
        fus$gene_names_5_1[w] <- fus$gene_names_5_2[w]
        w <- which(fus$gene_names_5_1=="")
        fus$gene_names_5_1[w] <- fus$cyt_5_1[w]
        w <- which(fus$gene_names_5_1=="")
        fus$gene_names_5_1[w] <- fus$cyt_5_2[w]

        fus$pos_3 <- paste0(fus$chr_3,':',fus$pos_3)
        fus$pos_5 <- paste0(fus$chr_5,':',fus$pos_5)

        if(detail){
          .dict <- dict %>% filter(table=='fusion' & show_detail)
        }else{
          .dict <- dict %>% filter(table=='fusion' & show_table)
        }
          rnm <- as.character(.dict$display_name)
          names(rnm) <- .dict$db_name
          fmt <- fus %>% select(.dict$db_name) %>% S4Vectors::rename(rnm)

        hide=c()
        renderpct=c()

        return(list(table=fmt,hide=hide,renderpct=renderpct))
}
