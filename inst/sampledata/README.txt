
Sample data for the R/qtl package

These files contain sample data in several formats, so that the user
may better understand how data may be formatted for import via the
read.cross function.  These are the same as the "listeria" data
set included with the R/qtl package.  

Note: Replace the "..." in the directory string to the appropriate
location of the sampledata directory (for example,
"/usr/local/lib/R/library/qtl" or "c:/R/rw1081/library/qtl").


1. "csv" format

    File:

        listeria.csv

    Data import:

        listeria.a <- read.cross("csv", ".../sampledata", "listeria.csv")

2. "csvr" format (rotated "csv")

   File:

        listeria_rot.csv

   Data import:

        listeria.a2 <- read.cross("csvr", ".../sampledata", "listeria_rot.csv")


3. "csvs" format (like "csv", but with separate files for phenotype
   and genotype data)

   Files:

        listeria_gen.csv         Genotype data
        listeria_phe.csv         Phenotype data

   Data import:

        listeria.a3 <- read.cross("csvs", ".../sampledata",
                                  "listeria_gen.csv", "listeria_phe.csv")


4. "csvsr" format (like "csvr", but both files are rotated)

   Files:

        listeria_gen_rot.csv     Genotype data
        listeria_phe_rot.csv     Phenotype data

   Data import:

        listeria.a4 <- read.cross("csvsr", ".../sampledata",
                                  "listeria_gen_rot.csv", 
                                  "listeria_phe_rot.csv")


5. "mm" (mapmaker) format

   Files:

       listeria_raw.txt   "raw" file (with genotype and phenotype data)
       listeria_map.txt    Genetic map information (markers must be in
                           order; map positions are not required)
       listeria_maps.txt   Genetic map information, as produced by 
                           Mapmaker/exp

   Data import:

       listeria.b <- read.cross("mm", ".../sampledata",
                                "listeria_raw.txt",,"listeria_map.txt")

       listeria.bb <- read.cross("mm", ".../sampledata",
                                 "listeria_raw.txt",,"listeria_maps.txt")


6. "qtx" (Mapmanager QTX) format

   File:

       listeria.qtx

   Data import:

       listeria.c <- read.cross("qtx", ".../sampledata", "listeria.qtx")


7. "qtlcart" (QTL Cartographer) format

   Files:

        listeria_qc_cro.txt    Genotype/phenotype data
        listeria_qc_map.txt    Genetic map information

   Data import:

        listeria.d <- read.cross("qtlcart", ".../sampledata", 
                                 "listeria_qc_cro.txt", "listeria_qc_map.txt")


8. "karl" format

    Files:

        gen.txt    Genotype data
        phe.txt    Phenotype data
        map.txt    Genetic map information

    Data import:

        listeria.e <- read.cross("karl", ".../sampledata",
                                 genfile="gen.txt", phefile="phe.txt",
                                 mapfile="map.txt")

    or just the following:  

        listeria.e <- read.cross("karl", ".../sampledata")
