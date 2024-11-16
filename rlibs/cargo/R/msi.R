mi_msi <- function(somatic, homopolymers, targets) {

    # limit homopolymer ranges to only those overlapping with the capture region
    hits <- findOverlaps(homopolymers, targets)
    homopolymer_ranges <- homopolymers[queryHits(hits)]

    # Initialize metadata column to count indels
    homopolymer_ranges$indels <- 0

    # Only keep poly-N indels
    rep_indels <- somatic[pass=="pass" & str==TRUE & str_ru %in% c("A", "T", "G", "C"),
        .(ID, chr, pos, ref, alt, str, str_ru, str_len, str_diff)
    ]

    # Convert to GRange
    if (nrow(rep_indels) > 0) {
        rep_indel_ranges <- GRanges(seqnames = rep_indels$chr,
                                    strand = "*",
                                    ranges = IRanges(start = rep_indels$pos, width = 2))
    } else {
        rep_indel_ranges <- GRanges(seqnames = NULL,
                                    strand = NULL,
                                    ranges = NULL)
    }

    # Find overlap between indels and homopolymer ranges
    indel_homopolymer_overlap <- findOverlaps(rep_indel_ranges, homopolymer_ranges)

    # Add a count to each homopolymer with an observed indel
    if (length(indel_homopolymer_overlap) > 0) {
        homopolymer_ranges[subjectHits(indel_homopolymer_overlap)]$indels <- 1
    }

    # Tally indels based on length of homopolymer
    if (length(indel_homopolymer_overlap) > 0) {
        homopolymer_indel_counts <- width(ranges(homopolymer_ranges[homopolymer_ranges$indels > 0])) %>% table()
    } else {
        homopolymer_indel_counts <- tibble(length = character(), indels = numeric())
    }

    # Tally total number of sites of each length
    homopolymer_location_counts <- width(ranges(homopolymer_ranges)) %>% table()

    # Join to get frequency of indels based on length
    homopolymer_freqs <- dplyr::tibble(
                        length = names(homopolymer_location_counts),
                        sites = as.numeric(homopolymer_location_counts)) %>%
            dplyr::left_join(dplyr::tibble(
                        length = names(homopolymer_indel_counts),
                        indels = as.numeric(homopolymer_indel_counts)),
                        by="length"
                        ) %>%
        dplyr::mutate(indels = tidyr::replace_na(indels, 0)) %>%
        dplyr::mutate(freq = indels/sites) %>%
        dplyr::mutate(length = as.numeric(length))

    # Score to report
    # % of homopolymers with length 7 to 20 that have an observed indel. Ignoring homopolymers
    # less than length 7 because they have lower mutation rates in MSI, but are more numerous
    # so dominate the calculation if included. Homopolymers greater than length 20 are extremely
    # rare in coding sequences so they contribute little and mess up the plots.
    homopolymer_freqs_cut <- homopolymer_freqs[homopolymer_freqs$length %in% 7:20, ]

    msi_score <- (sum(homopolymer_freqs_cut$indels, na.rm = TRUE))/(sum(homopolymer_freqs_cut$sites))*100
    msi <- list(mimsi_score=msi_score, mimsi_data=homopolymer_freqs_cut)
    return(msi)
}

msiVault <- function(vault, homopolymers) {
    somatic <- vault$tables$somatic
    targets <- import(vault$meta$targets)
    mi_msi(somatic, homopolymers, targets)
}

msiMetaVault <- function(mv) {
    homopolymers <- mv$anno$homopolymers
    foreach(id.sel=mv$meta$id) %do% {
        mv$meta[id==id]
        somatic <- mv$tables$somatic[id==id.sel]
        targets <- import(mv$meta[id==id.sel]$targets)
        mi_msi(somatic, homopolymers, targets)
    }
}