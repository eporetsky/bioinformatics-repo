# Using TIC (Total ion chromatogram, also known as TIS, total ion signal) 
# Divide each feature intensity by the sum of all features in that sample
# Then scaling by mean-centering and dividing by the STD of each variable

################################ LOAD DATA ################################
filename <- "data.tsv"
df <- read.table(filename, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
col_names <- colnames(df)           # Original column names in order to group rows
rep_names <- unname(unlist(df[1,])) # Original reps to save the transformed and scaled df

# Sometimes when transposing data or if data is transposed columns might be of character type
# To convert columns to numeric ([,2:x] to select numeric type): 
#     apply(df[,2:5], 1, function(x){as.numeric(x)})

######################## START DATA TRANSFORMATION ########################
# Log2 transform (data + 1)
df <- t(log2(df+1))
# Divide each column by sum of the column
df <- sweep(df, 2, colSums(df), FUN = '/')
# Transpose because scaling works on columns
df <- scale(t(df), center = TRUE, scale = TRUE)

########################### CALCULATE GROUP MEAN ##########################
dft <- as.data.frame(t(df))  # Easier to re-assign columns with data.frame type
dft["Rep"] <- col_names      # Easier to group_by with unique line names

# Groupby and genotype names and calculate the mean of every group for each column
#https://www.datanovia.com/en/blog/dplyr-how-to-compute-summary-statistics-across-multiple-columns/
dfg <- dft %>%
  group_by(Rep) %>%
  summarise(across(
    .cols = is.numeric, 
    .fns = list(Mean = mean), na.rm = TRUE, 
    .names = "{col}"
  ))

######################### WRITE DATA TO TSV FILE ##########################
# Use this write.table method to save data as TASSEL phenotype format for GWAS
# Useful if you want to rename the first cell of the table to anything
# https://stackoverflow.com/questions/2478352/write-table-writes-unwanted-leading-empty-column-to-header-when-has-rownames
write.table(data.frame("<trait>"=rownames(dfg),dfg, check.names = FALSE),"output.tsv", 
            row.names=FALSE, sep="\t", quote=FALSE, )

# With dplyr group_by function change the name of the "Rep" column to "<trait>"
colnames(dfg)[1] <- "<trait>"
write.table(dfg, "output.tsv", row.names=FALSE, sep="\t", quote=FALSE, )
