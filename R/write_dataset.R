format.ASV.tab <- function(ps) {
  if (as.logical(class(phyloseq::otu_table(ps))[1] == "otu_table") &&
      as.logical(taxa_are_rows(phyloseq::otu_table(ps)) == TRUE)) {
    asv.tab <- as.matrix(phyloseq::otu_table(ps))
  } else {
    asv.tab <- as.matrix(t(phyloseq::otu_table(ps)))
  }
  asv.tab
}

write.dataset.biom <-
  function(ps,
           filePATH,
           filePREFIX,
           writeFASTA = TRUE,
           rename = FALSE,
           useREFSEQ = TRUE) {
    #pull seqs from refseq slot or extract from ASV ID for fasta format
    if (isTRUE(useREFSEQ)) {
      #from phyloseq refseq slot
      f.onames <- phyloseq::refseq(ps)
    } else {
      f.onames <- phyloseq::taxa_names(ps)
    }
    
    if (isTRUE(rename)) {
      phyloseq::taxa_names(ps) <-
        paste("ASV", 1:length(phyloseq::taxa_names(ps)), sep = "")
      names(f.onames) <- paste0(">", phyloseq::taxa_names(ps))
    } else {
      names(f.onames) <- paste0(">", phyloseq::taxa_names(ps))
    }
    
    #generate biom file
    suppressWarnings(
      ps.b <- biomformat::make_biom(
        data = format.ASV.tab(ps),
        sample_metadata = as.data.frame(phyloseq::sample_data(ps)),
        observation_metadata = as.data.frame(phyloseq::tax_table(ps)),
        matrix_element_type = "int"
      )
    )
    
    #create output string
    if (isTRUE(writeFASTA)) {
      fa <- print(paste0(filePATH, filePREFIX, "_ASVs.fasta"))
    }
    bo <- print(paste0(filePATH, filePREFIX, "_ASV_table.biom"))
    
    #write output
    if (isTRUE(writeFASTA)) {
      write.table(
        x = f.onames,
        file = fa,
        quote = FALSE,
        sep = "\n",
        col.names = FALSE
      )
    }
    #biom export
    biomformat::write_biom(x = ps.b, biom_file = bo)
    
    #return phyloseq object with taxa renamed to ASV1, etc., if desired
    if (isTRUE(rename)) {
      return(ps)
    }
  }

write.dataset <-
  function(ps,
           filePATH,
           filePREFIX,
           writeFASTA = TRUE,
           rename = FALSE,
           useREFSEQ = TRUE) {
    #pull seqs from refseq slot or extract from ASV ID for fasta format
    if (isTRUE(useREFSEQ)) {
      #from phyloseq refseq slot
      f.onames <- phyloseq::refseq(ps)
    } else {
      f.onames <- phyloseq::taxa_names(ps)
    }
    
    if (isTRUE(rename)) {
      phyloseq::taxa_names(ps) <-
        paste("ASV", 1:length(phyloseq::taxa_names(ps)), sep = "")
      names(f.onames) <- paste0(">", phyloseq::taxa_names(ps))
    } else {
      names(f.onames) <- paste0(">", phyloseq::taxa_names(ps))
    }
    
    
    #generate asv table formatted for biom generation
    asv.tab <- format.ASV.tab(ps)
    suppressWarnings(asv.tab <- as.matrix(asv.tab))
    cb <- as.matrix(cbind(rownames(asv.tab), asv.tab))
    rcb <- as.matrix(rbind(colnames(cb), cb))
    rcb[1, 1] <- "#ASVID"
    rownames(rcb) <- NULL
    colnames(rcb) <- NULL
    
    #generate tax table formatted for biom generation
    tax.tab <- as.data.frame(phyloseq::tax_table(ps))
    tax.tab$taxonomy <-
      tidyr::unite(tax.tab, "out", c(colnames(tax.tab)), sep = ";")
    cbt <- as.matrix(cbind(rownames(tax.tab), tax.tab$taxonomy))
    rcbt <- as.matrix(rbind(c("#ASVID", "taxonomy"), cbt))
    rownames(cbt) <- NULL
    colnames(cbt) <- NULL
    
    #generate sampledf table formatted for biom generation
    samdf <- suppressWarnings(as.matrix(phyloseq::sample_data(ps)))
    cbs <- as.matrix(cbind(rownames(samdf), samdf))
    rcbs <- as.matrix(rbind(colnames(cbs), cbs))
    rcbs[1, 1] <- "#SampleID"
    rownames(rcbs) <- NULL
    colnames(rcbs) <- NULL
    
    #create output string
    if (isTRUE(writeFASTA)) {
      fa <- print(paste0(filePATH, filePREFIX, "_ASVs.fasta"))
    }
    otb <- print(paste0(filePATH, filePREFIX, "_ASV_table.txt"))
    ttb <- print(paste0(filePATH, filePREFIX, "_ASV_taxonomy.txt"))
    stb <- print(paste0(filePATH, filePREFIX, "_sample_data.txt"))
    
    
    #write output
    #ASV fasta
    if (isTRUE(writeFASTA)) {
      write.table(
        x = f.onames,
        file = fa,
        quote = FALSE,
        sep = "\n",
        col.names = FALSE
      )
    }
    #asv.tab
    write.table(
      x = rcb,
      file = otb,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
    #tax.tab
    write.table(
      x = rcbt,
      file = ttb,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
    #sampledf
    write.table(
      x = rcbs,
      file = stb,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
    
    #return phyloseq object with taxa renamed to ASV1, etc., if desired
    if (isTRUE(rename)) {
      return(ps)
    }
  }