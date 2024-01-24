## Working on the DArT data
## the SNP data file
## Reading DArT Files into a genlight Object using read.dart()

### first, install the dartR package
install.packages("devtools")
library(devtools)
install_github("green-striped-gecko/dartR")
library(dartR)
gl.install.vanilla.dartR() 

##1. Filtering (95% reproducibility, >80% call rate,  monomorphic loci, and removing individuals with < 50% call rate )
gl <- gl.read.dart(filename = "Report_DAfg23-2937_SNP_2.csv", ind.metafile = "pop_1.csv")
gl
gl <- gl.compliance.check(gl) # renders the genlight object compliant with dartR
nLoc(gl) # 10219 initial SNPS
gl2 <- gl.filter.secondaries(gl) # filter out SNPs that share a sequence tag, except one retained at random

nLoc(gl2) #produced 9819 SNPs
gl3 <- gl.filter.rdepth(gl2) # filter out loci with exceptionally low or high read depth (coverage)
nLoc(gl3) #produced 8549 SNPs

gl4 <- gl.filter.reproducibility(gl3, threshold = 0.95) # filter out loci for which the reproducibility (strictly repeatability) is less than threshold = 0.95
gl4 <- gl.recalc.metrics(gl4) # recalculating the locus metrics after filtering because the initial call rate parameter will no longer be accurate after the filtering out individuals or populations
nLoc(gl4) #produced 7530 SNPs

gl5 <- gl.filter.callrate(gl4, method = "ind", threshold = 0.50) #Filter individuals on call rate (threshold =50% )deleted individuals. Individuals deleted (CallRate <=  0.5 ):AG025[Tesso_Borite_p], AG041[Wonsho_Abbo_p], AG057[Wonsho_Gudumale_p], AG002[Tesso_Sodicho_p], AG090[Wondogenet_college_p], AG026[Tesso_Borite_p], AG042[Wonsho_Abbo], AG050[Wonsho_Abbo], AG058[Wonsho_Gudumale_p], AG083[Wondogenet_college_p], AG027[Tesso_Borite], AG043[Wonsho_Abbo], AG060[Wonsho_Gudumale_p], AG029[Tesso_Borite_p], AG037[Wonsho_Abbo_p], AG053[Wonsho_Abbo_p], AG022[Tesso_Borite_p], AG030[Tesso_Borite_p], AG038[Wonsho_Abbo_p], AG079[Wondogenet_college_p], AG048[Wonsho_Abbo], AG064[Wonsho_Gudumale_p], AG103[Auger_p], AG127[ME_Gedam_p, AG135[ME_Gedam], AG143[Deban_Amba_den_p], AG151[Deban_Amba_den], AG128[ME_Gedam], AG136[ME_Gedam, AG169[Bera_Tedicho_p], AG137[ME_Gedam_p], AG106[Auger_p], AG138[ME_Gedam_p], AG154[Deban_Amba_den], AG131[ME_Gedam], AG155[Deban_Amba_den_p], AG172[Bera_Tedicho_p], AG140[ME_Gedam_p], AG148[Deban_Amba_den], AG173[Bera_Tedicho_p], AG181[Bera-Tedicho], AG109[Auger_p], AG141[Deban_Amba_den], AG149[Deban_Amba_den_p], AG174[Bera_Tedicho_p], AG142[Deban_Amba_den], AG150[Deban_Amba_den_p]

gl5 <- gl.recalc.metrics(gl5)
nLoc(gl5) # produced 7530 SNPs
nInd(gl5) # 138 individuals, reduced from the initial 185

gl6 <- gl.filter.callrate(gl5, method="loc", threshold = 0.8)  #filter out loci with callrate lower than 80%
gl6 <- gl.recalc.metrics(gl6)

nLoc(gl6) # produced 2198 SNPs

gl7 <- gl.filter.monomorphs(gl6) #filter out monomorphic loci
gl7 <- gl.recalc.metrics(gl7)
gl7
nLoc(gl7) # produced 1900 SNPs

m <- as.matrix(gl7) #save genlight object as matrix
write.csv(m,file="filtered_loci_metadata.csv") # saves as .csv
gl.save(gl7,file="afrocarpus.Rdata") # saves the the data as an R object filtered with 95% reproducibility, 80 % call rate, deleted individuals with less than 50% call rate, and monomorphic loci filtered out. it resulted in 1900SNPS 

