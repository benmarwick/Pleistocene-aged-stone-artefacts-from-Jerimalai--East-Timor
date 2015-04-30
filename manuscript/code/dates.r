# code for getting calibrations of dates from OxCal

# from https://gist.github.com/benmarwick/7916697


# read data in, three columns Name = lab code, Date = radiocarbon age, Uncertainty = error


dates <- read.table(header = TRUE, text = 
"Name Date Uncertainty
Wk-19228  124	32
Wk-19229	4962	50
Wk-19230	4580	42
Wk-18157	4867	42
Wk-18158	5595	43
Wk-18159	5694	45
Wk-17832	5939	45
Wk-19316	6118	41
Wk-19231	8879	78
Wk-30500	6223	26
Wk-30501	5575	27
Wk-30502	13901	45
Wk-30503	9457	32
Wk-18160	14007	146
Wk-30504	13778	43
Wk-19232	35387	534
Wk-17833	37267	453")


# construct OxCal format
oxcal_format <- paste0('R_Date(\"',  gsub("^\\s+|\\s+$", "", dates$Name), '\",', dates$Date, ',', dates$Uncertainty, ');')
# inspect
cat(oxcal_format)

# write formatted dates to text file
write.table(oxcal_format, file = 'oxcal_format.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)

# find location of text file
getwd()

# Now 95% ready for pasting into OxCal batch conversion
# at https://c14.arch.ox.ac.uk/oxcal/OxCal.html

# In OxCal, File -> New then click the table view button (fourth along the top, no tool tips on hover sadly!) and paste in the dates between this (without #)

# Plot()
# {
# ...insert dates here...
# };

# Then edit format and settings to return results in BP (before running) and get median and sigma (after running)
# using IntCal 13 curve
# then click File -> Run, then File -> Save As to get calibrated dates in a CSV
