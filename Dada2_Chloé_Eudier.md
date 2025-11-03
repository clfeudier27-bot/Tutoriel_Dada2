Tutorial DADA2
================
Chloé Eudier
2025-10-23

# Chargement des données

### Chargement de la librarie dada2 pour utiliser les commandes nécessaires à l’analyse de données de séquençage.

``` r
library(dada2)
```

    ## Loading required package: Rcpp

``` r
packageVersion("dada2")
```

    ## [1] '1.28.0'

### Chargement des données et création de la variable path pour indiquer le chemin d’accès aux lectures utilisées.

``` r
path<-"~/MiSeq_SOP"
list.files(path)
```

    ##  [1] "F3D0_S188_L001_R1_001.fastq"   "F3D0_S188_L001_R2_001.fastq"  
    ##  [3] "F3D1_S189_L001_R1_001.fastq"   "F3D1_S189_L001_R2_001.fastq"  
    ##  [5] "F3D141_S207_L001_R1_001.fastq" "F3D141_S207_L001_R2_001.fastq"
    ##  [7] "F3D142_S208_L001_R1_001.fastq" "F3D142_S208_L001_R2_001.fastq"
    ##  [9] "F3D143_S209_L001_R1_001.fastq" "F3D143_S209_L001_R2_001.fastq"
    ## [11] "F3D144_S210_L001_R1_001.fastq" "F3D144_S210_L001_R2_001.fastq"
    ## [13] "F3D145_S211_L001_R1_001.fastq" "F3D145_S211_L001_R2_001.fastq"
    ## [15] "F3D146_S212_L001_R1_001.fastq" "F3D146_S212_L001_R2_001.fastq"
    ## [17] "F3D147_S213_L001_R1_001.fastq" "F3D147_S213_L001_R2_001.fastq"
    ## [19] "F3D148_S214_L001_R1_001.fastq" "F3D148_S214_L001_R2_001.fastq"
    ## [21] "F3D149_S215_L001_R1_001.fastq" "F3D149_S215_L001_R2_001.fastq"
    ## [23] "F3D150_S216_L001_R1_001.fastq" "F3D150_S216_L001_R2_001.fastq"
    ## [25] "F3D2_S190_L001_R1_001.fastq"   "F3D2_S190_L001_R2_001.fastq"  
    ## [27] "F3D3_S191_L001_R1_001.fastq"   "F3D3_S191_L001_R2_001.fastq"  
    ## [29] "F3D5_S193_L001_R1_001.fastq"   "F3D5_S193_L001_R2_001.fastq"  
    ## [31] "F3D6_S194_L001_R1_001.fastq"   "F3D6_S194_L001_R2_001.fastq"  
    ## [33] "F3D7_S195_L001_R1_001.fastq"   "F3D7_S195_L001_R2_001.fastq"  
    ## [35] "F3D8_S196_L001_R1_001.fastq"   "F3D8_S196_L001_R2_001.fastq"  
    ## [37] "F3D9_S197_L001_R1_001.fastq"   "F3D9_S197_L001_R2_001.fastq"  
    ## [39] "filtered"                      "HMP_MOCK.v35.fasta"           
    ## [41] "Mock_S280_L001_R1_001.fastq"   "Mock_S280_L001_R2_001.fastq"  
    ## [43] "mouse.dpw.metadata"            "mouse.time.design"            
    ## [45] "stability.batch"               "stability.files"

### Lecture du nom des fichiers et création d’une liste des lectures forward (fnFS) et reverse (fnRs) par la manipulation de leur chaîne de caractères variables.

``` r
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
```

### Création d’une variable afin d’extraire le nom des échantillons auxquels appartiennent les couples de lectures.

``` r
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

# Etude des profils de qualité

#### Visualisation de la qualité de chaque séquence forward (fwd) grâce au score qualité de chaque nucléotide la composant.

``` r
plotQualityProfile(fnFs[1:2])
```

![](CC1_ADM_Chloé_Eudier_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

##### On obtient une heat map nous indiquant la fréquence de chaque score de qualité pour chaque position. Sur les graphiques on observe le score de qualité moyen par la ligne verte, les quartiles de la distributions des scores par les lignes oranges et la proportion des lectures allant jusqu’à au moins cette position.

##### On peut observer que les séquences fwd sont de bonnes qualité en regardant la ligne verte et les lignes oranges qui se maintiennent à un score entre 30 et 40. La qualité diminue cependant sur les dernières bases des séquences donc il faudra couper les 10 derniers nucléotides des séquences (à partir de la position 240) pour augmenter le score.

#### Visualisation de la qualité de chaque séquence reverse (rs) grâce au score qualité de chaque nucléotide la composant.

``` r
plotQualityProfile(fnRs[1:2])
```

![](CC1_ADM_Chloé_Eudier_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

##### En revanche, on observe que les séquences rs sont de moins bonne qualité. Il faudra les couper à partir de la position 160 qui est celle à partir de laquelle la qualité des séquences chute.

# Tronquage et Filtrage

### Création d’objets (filtFs et filtRs) servant à stocker les séquences filtrées et dont les noms sont rattachés à ceux des échantillons.

``` r
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

### Filtrage des données en jouant sur différents paramètres de sorte à obtenir un meilleur score de qualité de séquence.

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
```

    ##                               reads.in reads.out
    ## F3D0_S188_L001_R1_001.fastq       7793      7113
    ## F3D1_S189_L001_R1_001.fastq       5869      5299
    ## F3D141_S207_L001_R1_001.fastq     5958      5463
    ## F3D142_S208_L001_R1_001.fastq     3183      2914
    ## F3D143_S209_L001_R1_001.fastq     3178      2941
    ## F3D144_S210_L001_R1_001.fastq     4827      4312

##### Grâce à cet outil, une petite partie des séquences a été retirée pour augmenter le score de qualité de celles-ci. Donc comme expliqué précédemment, les séquences fwd ont été coupées sur les 10 derniers nucléotides et les rs sur les 90 derniers car le score de qualité des bases a chuté bien avant. On obtient ainsi la liste des lectures avant et après filtrage.

# Aprentissage des taux d’erreurs

### Estimation du nombre d’erreurs de séquençage dans les échantillons pour différencier les séquences mutées de celles erronées. Le modèle d’erreurs est calculé en alternant l’estimation des taux d’erreur et l’inférence de la composition de l’échantillon jusqu’à ce qu’ils convergent vers une solution cohérente.

#### Nombre de bases utilisées par le modèle pour les séquences fwd

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 33514080 total bases in 139642 reads from 20 samples will be used for learning the error rates.

#### Nombre de bases utilisées par le modèle pour les séquences rs

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 22342720 total bases in 139642 reads from 20 samples will be used for learning the error rates.

#### Visualisation des erreurs estimées (pour chaque transition de base possible).

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    ## Transformation introduced infinite values in continuous y-axis

![](CC1_ADM_Chloé_Eudier_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

##### Les taux d’erreurs observés pour chaque score de qualité sont représentés par un point noir. La ligne noire indique les taux d’erreur estimés après convergence de l’algorithme d’apprentissage automatique. La ligne rouge indique les taux d’erreur attendus selon la définition nominale du score de qualité.

##### Comme il est possible de le voir sur les graphiques, les taux d’erreurs estimés correspondent bien aux taux observés car les points noirs suivent la ligne rouge pour chaque transition et les taux d’erreurs diminuent avec l’augmentation de la qualité.

# Inférence des échantillons

### Détection des inférences dans les séquences filtrées afin de différencier les vraies séquences biologiques (variants biologiques) des erreurs de séquençage en utilisant un modèle d’erreurs paramétrique.

#### Inférences dans les séquences fwd filtrées

``` r
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1979 unique sequences.
    ## Sample 2 - 5299 reads in 1639 unique sequences.
    ## Sample 3 - 5463 reads in 1477 unique sequences.
    ## Sample 4 - 2914 reads in 904 unique sequences.
    ## Sample 5 - 2941 reads in 939 unique sequences.
    ## Sample 6 - 4312 reads in 1267 unique sequences.
    ## Sample 7 - 6741 reads in 1756 unique sequences.
    ## Sample 8 - 4560 reads in 1438 unique sequences.
    ## Sample 9 - 15637 reads in 3590 unique sequences.
    ## Sample 10 - 11413 reads in 2762 unique sequences.
    ## Sample 11 - 12017 reads in 3021 unique sequences.
    ## Sample 12 - 5032 reads in 1566 unique sequences.
    ## Sample 13 - 18075 reads in 3707 unique sequences.
    ## Sample 14 - 6250 reads in 1479 unique sequences.
    ## Sample 15 - 4052 reads in 1195 unique sequences.
    ## Sample 16 - 7369 reads in 1832 unique sequences.
    ## Sample 17 - 4765 reads in 1183 unique sequences.
    ## Sample 18 - 4871 reads in 1382 unique sequences.
    ## Sample 19 - 6504 reads in 1709 unique sequences.
    ## Sample 20 - 4314 reads in 897 unique sequences.

#### Inférences dans les séquences rs filtrées

``` r
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1660 unique sequences.
    ## Sample 2 - 5299 reads in 1349 unique sequences.
    ## Sample 3 - 5463 reads in 1335 unique sequences.
    ## Sample 4 - 2914 reads in 853 unique sequences.
    ## Sample 5 - 2941 reads in 880 unique sequences.
    ## Sample 6 - 4312 reads in 1286 unique sequences.
    ## Sample 7 - 6741 reads in 1803 unique sequences.
    ## Sample 8 - 4560 reads in 1265 unique sequences.
    ## Sample 9 - 15637 reads in 3414 unique sequences.
    ## Sample 10 - 11413 reads in 2522 unique sequences.
    ## Sample 11 - 12017 reads in 2771 unique sequences.
    ## Sample 12 - 5032 reads in 1415 unique sequences.
    ## Sample 13 - 18075 reads in 3290 unique sequences.
    ## Sample 14 - 6250 reads in 1390 unique sequences.
    ## Sample 15 - 4052 reads in 1134 unique sequences.
    ## Sample 16 - 7369 reads in 1635 unique sequences.
    ## Sample 17 - 4765 reads in 1084 unique sequences.
    ## Sample 18 - 4871 reads in 1161 unique sequences.
    ## Sample 19 - 6504 reads in 1502 unique sequences.
    ## Sample 20 - 4314 reads in 732 unique sequences.

#### Evaluation du nombre d’ASVs dans l’échantillon n°1

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 128 sequence variants were inferred from 1979 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

##### L’algorithme a trouvé 128 variants de séquence réels à partir des 1979 séquences du premier échantillon.

# Fusion des séquences paires

### Fusion des séquences R1 et R2 pour obtenir les amplicons en entier et observer leur fiabilité. Pour que la fusion se fasse, il faut que les séquences se chevauchent sur au moins 12 nucléotides identiques entre les 2 séquences.

``` r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

    ## 6540 paired-reads (in 107 unique pairings) successfully merged out of 6891 (in 197 pairings) input.

    ## 5028 paired-reads (in 101 unique pairings) successfully merged out of 5190 (in 157 pairings) input.

    ## 4986 paired-reads (in 81 unique pairings) successfully merged out of 5267 (in 166 pairings) input.

    ## 2595 paired-reads (in 52 unique pairings) successfully merged out of 2754 (in 108 pairings) input.

    ## 2553 paired-reads (in 60 unique pairings) successfully merged out of 2785 (in 119 pairings) input.

    ## 3646 paired-reads (in 55 unique pairings) successfully merged out of 4109 (in 157 pairings) input.

    ## 6079 paired-reads (in 81 unique pairings) successfully merged out of 6514 (in 198 pairings) input.

    ## 3968 paired-reads (in 91 unique pairings) successfully merged out of 4388 (in 187 pairings) input.

    ## 14233 paired-reads (in 143 unique pairings) successfully merged out of 15355 (in 352 pairings) input.

    ## 10528 paired-reads (in 120 unique pairings) successfully merged out of 11165 (in 278 pairings) input.

    ## 11154 paired-reads (in 137 unique pairings) successfully merged out of 11797 (in 298 pairings) input.

    ## 4349 paired-reads (in 85 unique pairings) successfully merged out of 4802 (in 179 pairings) input.

    ## 17431 paired-reads (in 153 unique pairings) successfully merged out of 17812 (in 272 pairings) input.

    ## 5850 paired-reads (in 81 unique pairings) successfully merged out of 6095 (in 159 pairings) input.

    ## 3716 paired-reads (in 86 unique pairings) successfully merged out of 3894 (in 147 pairings) input.

    ## 6865 paired-reads (in 99 unique pairings) successfully merged out of 7191 (in 187 pairings) input.

    ## 4426 paired-reads (in 67 unique pairings) successfully merged out of 4603 (in 127 pairings) input.

    ## 4576 paired-reads (in 101 unique pairings) successfully merged out of 4739 (in 174 pairings) input.

    ## 6092 paired-reads (in 109 unique pairings) successfully merged out of 6315 (in 173 pairings) input.

    ## 4269 paired-reads (in 20 unique pairings) successfully merged out of 4281 (in 28 pairings) input.

``` r
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                       sequence
    ## 1 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAGTGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTATCGAACAGG
    ## 2 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACAGG
    ## 3 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTGTTAAGTCAGCGGTCAAATGTCGGGGCTCAACCCCGGCCTGCCGTTGAAACTGGCGGCCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ## 4 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTTTTAAGTCAGCGGTAAAAATTCGGGGCTCAACCCCGTCCGGCCGTTGAAACTGGGGGCCTTGAGTGGGCGAGAAGAAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCCTTCCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCGAACAGG
    ## 5 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGACTCTCAAGTCAGCGGTCAAATCGCGGGGCTCAACCCCGTTCCGCCGTTGAAACTGGGAGCCTTGAGTGCGCGAGAAGTAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCCTACCGGCGCGCAACTGACGCTCATGCACGAAAGCGTGGGTATCGAACAGG
    ## 6 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGATGCCAAGTCAGCGGTAAAAAAGCGGTGCTCAACGCCGTCGAGCCGTTGAAACTGGCGTTCTTGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1       579       1       1    148         0      0      1   TRUE
    ## 2       470       2       2    148         0      0      2   TRUE
    ## 3       449       3       4    148         0      0      1   TRUE
    ## 4       430       4       3    148         0      0      2   TRUE
    ## 5       345       5       6    148         0      0      1   TRUE
    ## 6       282       6       5    148         0      0      2   TRUE

##### On obtient ainsi la liste de données pour chaque échantillon qui comprend : la séquence entière, son abondance et les indices de variants de séquence (ASV) fwd et rs qui ont été fusionnées. De plus, les lectures appariées qui ne se chevauchaient pas exactement ont été supprimées pour réduire le risque de résultats erronés.

# Création d’une table d’ASV

### Construction d’une table d’ASV pour observer l’abondance des séquences dans les échantillons

#### Evaluation du nombre d’ASV parmi l’ensemble des échantillons

``` r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

    ## [1]  20 293

#### Etude de la distribution de la longueur des séquences

``` r
table(nchar(getSequences(seqtab)))
```

    ## 
    ## 251 252 253 254 255 
    ##   1  88 196   6   2

##### Comme on peut le voir dans le tableau, on obtient 293 ASVs de différentes longueurs parmi les 20 échantillons. En effet, on trouve 1 séquence de 251 nucléotides (NT), 88 de 252, 196 de 253 NT, 6 de 254 NT et 2 de 255 NT.

##### Les longueurs des séquences fusionnées se situent toutes dans la fourchette attendue pour cet amplicon V4.

# Eliminitaion des ASVs chimériques

### Suppression des ASVs chimériques (séquences non biologiques). Les séquences chimériques sont identifiées si elles peuvent être reconstruites exactement en combinant un segment gauche et un segment droit provenant de deux séquences « parentales » plus abondantes.

#### Identification du nombre d’ASVs chimériques

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 61 bimeras out of 293 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]  20 232

##### Ici,il y’a 61 ASVs qui sont en réalité des chimères sur les 293 trouvées précédemment.

#### Evaluation de la proportion des chimères parmi l’ensemble ASVs

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.9640374

##### Les chimères représentent donc environ 21 % des variants de séquences fusionnées, mais si l’abondance de ces variants est prise en compte, les chimères ne représentent qu’environ 4 % des lectures de séquence fusionnées puisque 96% des séquences ont été conservées.

# Tableau de suivi

### Observation du nombre de séquences éliminées à chaque étape pour voir si un maximum de lectures présentent au départ ont été conservées

``` r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

    ##        input filtered denoisedF denoisedR merged nonchim
    ## F3D0    7793     7113      6976      6979   6540    6528
    ## F3D1    5869     5299      5227      5239   5028    5017
    ## F3D141  5958     5463      5331      5357   4986    4863
    ## F3D142  3183     2914      2799      2830   2595    2521
    ## F3D143  3178     2941      2822      2868   2553    2519
    ## F3D144  4827     4312      4151      4228   3646    3507

##### Une grande majorité des séquences a été conservée ce qui indique qu’il n’y a pas de problèmes.

# Assigniation taxonomique des séquences

### Assigniation de taxonomie pour chaque séquences afin d’obtenir des informations sur les différents rangs (Domaine jusqu’à l’espèce si c’est possible) par la comparaison avec des séquences de référence connues

``` r
taxa <- assignTaxonomy(seqtab.nochim,"~/Data/silva_nr99_v138.2_toSpecies_trainset.fa.gz", multithread=TRUE)
```

#### Affichage des résultats de l’assignation

``` r
taxa.print <- taxa 
rownames(taxa.print) <- NULL
head(taxa.print)
```

    ##      Kingdom    Phylum         Class         Order           Family          
    ## [1,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [2,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [3,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [4,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [5,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Bacteroidaceae"
    ## [6,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ##      Genus         Species     
    ## [1,] NA            NA          
    ## [2,] NA            NA          
    ## [3,] NA            NA          
    ## [4,] NA            NA          
    ## [5,] "Bacteroides" "caecimuris"
    ## [6,] NA            NA

##### La majorité des séquences des différents échantillons ont été assignées au phylum des Bacteroidetes. Dans le cas où l’outil ne peut pas assigner précisement des séquences à une espèce spécifique, il va seulement l’assigner aux rangs supérieurs tels que le genre ou la famille.

# Evaluation de la précision Dada2

### Comparaison de la composition de l’échantillon en ASV inferées à la composition attendue de la communauté de référence (Mock community)

``` r
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE)
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
```

    ## DADA2 inferred 20 sample sequences present in the Mock community.

``` r
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
```

    ## Of those, 20 were exact matches to the expected reference sequences.

##### La Mock community contient 20 souches bactériennes et l’algorithme a identifié 20 ASV qui correspondent exactement aux génomes de référence des membres attendus de la communauté. Le taux d’erreur résiduel après le pipeline DADA2 pour cet échantillon est donc de 0 %.
