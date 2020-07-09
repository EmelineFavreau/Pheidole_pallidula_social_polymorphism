# Load libraries; install from scratch if needed
libraries <- c("taxize")
for (lib in libraries) {
    if (require(package = lib, character.only = TRUE)) {
        print("Successful")
    } else {
        print("Installing")
        install.packages(lib)
        library(lib, character.only = TRUE )
    }
}

# define function to lookup phylum from ncbi using taxid (move to another file)
ranks_from_NCBI_taxid <- function(identifier) {
    # ensure that we have a usable identifier
    if (is.null(identifier) | is.na(identifier)) {
        return(c(NA, NA, NA, NA)) 
    }
    if (identifier == 0 | identifier == "") {
        return(c(NA, NA, NA, NA)) 
    }
    identifier <- as.character(identifier)
    if (identifier != as.numeric(identifier)) { 
        stop("Converting to and from numeric leads to ambiguity. ", 
             "There must be a problem with id formatting for: ", identifier) 
    }
    if (length(identifier) != 1) {
        stop("Just one id at a time please.")
    }
    
    # get full taxonomic hierarchy from NCBI
    ncbi_info <- classification(x = identifier, db = "ncbi")
    ncbi_info_table <- data.frame(ncbi_info[1]) # rest is confusing
    if (nrow(ncbi_info_table) <= 1) {
        warning("Got empty table for: ", identifier, "\n")
        return(c(NA, NA, NA, NA))
    }
    colnames(ncbi_info_table) <- c("name", "rank", "id")
    
    # extract phylum class and order (and sanity check we only have 1 of each)
    phylum <- ncbi_info_table$name[ ncbi_info_table$rank == "phylum" ]
    class <- ncbi_info_table$name[ ncbi_info_table$rank == "class" ]
    order <- ncbi_info_table$name[ ncbi_info_table$rank == "order" ]
    
    
    if (length(phylum) != 1) {
        warning("Didn't get 1 phylum as I expected for: ", identifier, "\n")
        phylum <- NA
    } 
    if (length(class) != 1) {
        warning("Didn't get 1 class as I expected for: ", identifier, "\n")
        class <- NA
    }
    if (length(order) != 1) {
        warning("Didn't get 1 order as I expected for: ", identifier, "\n")
        order <- NA
    }
    
    # other will be the last row in the table if there is no phylum, class or order in the table
    if (is.na(phylum) & is.na(class) & is.na(order)) {
        warning("No phylum, class or order for: ", identifier, ", returning lowest common taxon\n")
        other <- ncbi_info_table$name[nrow(ncbi_info_table)]
    } else {
        other <-  NA
    }
    
    return(c(phylum, class, order, other))
}
