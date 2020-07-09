# aim: understand if any id groupings are due to regional clusters
# import the names of the samples, based on their id
pallidula_samples <- read.csv(file = "result/pallidula_samples.txt", header = FALSE, sep = "\t")
colnames(pallidula_samples) <- c("samples")
pallidula_samples$id <- rep("pallidula", times = nrow(pallidula_samples)) 
sp_samples <- read.csv(file = "result/sp_samples.txt", header = FALSE, sep = "\t")
colnames(sp_samples) <- c("samples")
sp_samples$id <- rep("sp", times = nrow(sp_samples))
missing_id_samples <- read.csv(file = "result/samples_wo_id.txt", header = FALSE, sep = "\t")
colnames(missing_id_samples) <- c("samples")
missing_id_samples$id <- rep("missing_id", times = nrow(missing_id_samples))
ambiguous_pallidula <- read.csv(file = "result/ambiguous_samples/ambiguous_pallidula.txt", header = FALSE, sep = "\t")
colnames(ambiguous_pallidula) <- c("samples")
ambiguous_pallidula$id <- rep("ambiguous_pallidula", times = nrow(ambiguous_pallidula))
ambiguous_sp <- read.csv(file = "result/ambiguous_samples/ambiguous_sp.txt", header = FALSE, sep = "\t")
colnames(ambiguous_sp) <- c("samples")
ambiguous_sp$id <- rep("ambiguous_sp", times = nrow(ambiguous_sp))

# put all df together
df <- rbind(pallidula_samples, sp_samples, missing_id_samples, ambiguous_pallidula, ambiguous_sp)
write.csv(df, file = "identification_results")
df$samples <- as.character(df$samples)
df$colony <- gsub(pattern = "-[A-Z]", x = df$samples, replacement = "")
View(df)

# import the location data
location_df <- read.csv(file = "S1_pheidole_paper.csv", header = TRUE, sep = ",")
location_df <- location_df[, 1:8]
location_df$Sample.ID <- as.character(location_df$Sample.ID)
colnames(location_df) <- c("colony", "Social.type", "Country.Origin", "Region.Origin", "Location", "Latitude", "Longitude", "Elevation..m.")

# make the wanted table
# samples |  Region.Origin | id
new_df <- merge(x = df, y = location_df, by = intersect(names(df), names(location_df)))
colnames(new_df)

# where are all the sp? Not in France: "Tuscany" "Corsica" 
unique(as.character(new_df[new_df$id == "sp", 6]))

# where are all the ambiguous_sp? "Corsica" "Tuscany"
unique(as.character(new_df[new_df$id == "ambiguous_sp", 6]))

# where are all the pallidula? "Occitanie" "Andalucia"
unique(as.character(new_df[new_df$id == "pallidula", 6]))

# where are all the ambiguous_pallidula? "Occitanie" "Tuscany"
unique(as.character(new_df[new_df$id == "ambiguous_pallidula", 6]))

# where are all the missing id? "Catalonia" "Andalucia"
unique(as.character(new_df[new_df$id == "missing_id", 6]))

# In summary, the barcode "sp" flags only some Italian and Corsican samples
# the barcode "pallidula" flags samples from Bruniquel and Spain
# the ambiguous barcodes flag both sides of the Alps.
# There is some indication that the barcode "pallidula" comes from France (as expected)
# and "sp" from Italy (as expected)
