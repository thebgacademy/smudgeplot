# Understanding k-mers and ploidy using Smudgeplot

#BGA24/sessions #GitPod #Tools #Kmers

This session is part of [**Biodiversity Genomics Academy 2024**](https://thebgacademy.org)

[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/thebgacademy/smudgeplot) 

## Session Leader(s)

Kamil S. Jaron; Amjad Khalaf

Tree of Life, Wellcome Sanger Institute

[Google slides](https://docs.google.com/presentation/d/10g_TTIYl-v_F9jzZrKu5w7aTJs7tNWjHoL-6pr9NSjA/edit?usp=sharing)

## Description

By the end of this session you will be able to:

1. Understand how Smudgeplot estimates ploidy
2. Appreciate strengths and weakness of Smudgeplot
3. Run Smudgeplot and understand the input parameters
4. Critically evaluate a Smudgeplot

## Prerequisites

1. Understanding of linux command line basics
2. Knowledge of basic genome biology
3. (optional) read the smudgeplot sections of <https://www.nature.com/articles/s41467-020-14998-3>

!!! warning "Please make sure you MEET THE PREREQUISITES and READ THE DESCRIPTION above"

    You will get the most out of this session if you meet the prerequisites above.

    Please also read the description carefully to see if this session is relevant to you.
    
    If you don't meet the prerequisites or change your mind based on the description or are no longer available at the session time, please email tol-training at sanger.ac.uk to cancel your slot so that someone else on the waitlist might attend.

## Tutorial


Have you ever sequenced something not-well studied? Something that might show strange genomic signatures? Smudgeplot is a visualisation technique for whole-genome sequencing reads from a single individual. The visualisation techique is based on the idea of het-mers. Het-mers are k-mer pairs that are exactly one nucleotide pair away from each other, while forming a unique pair in the sequencing dataset. These k-mers are assumed to be mostly representing two alleles of a heterozygous, but potentially can also show pairing of imperfect paralogs, or sequencing errors paired up with a homozygous genomic k-mer. Nevertheless, the predicted ploidy by smudgeplot is simply the ploidy with the highest number of k-mer pairs (if a reasonable estimate must be evaluated for each individual case!).

### Constructing a FastK database and finding hetmers

These k-mer analyses operate on raw, or trimmed sequencing reads. From those we generate a k-mer database using [FastK](https://github.com/thegenemyers/FASTK). FastK is currently the fastest k-mer counter out there and the only supported by the lastest version of smudgeplot*. This database contains an index of all the k-mers and their coverages in the sequencing readset.

*Note: The previous versions of smudgeplot (up to 2.5.0) were operating on k-mer "dumps" flat files you can generate with any counter you like. You can imagine that text files are very inefficient to operate on. The new version is operating directly on the optimised k-mer database instead.

To learn how to build a FastK database and learn how to find hetmers, we will use relatively small yeast data.

```
SAMPLE=SRR3265401
FastK -v -t1 -k31 -M120 -T4 smudgeplot/data/Scer/"$SAMPLE"* -Nsmudgeplot/data/Scer/FastK_Table_"$SAMPLE"
```

Now, that you have a database, you can search for k-mer pairs, but I would advice to take a moment and look at a k-mer spectra first. You can get k-mer spectra from the database using `Histex`, a different tool from the same suite. The generated k-mer histogram can be used to fit GenomeScope model.

```
Histex -G smudgeplot/data/Scer/FastK_Table_"$SAMPLE" > smudgeplot/data/Scer/"$SAMPLE"_k31.hist
```

(optional) run GenomeScope on that k-mer histogram. You can use preinstalled genomescope, or upload the histogram to genomescope2 webserver: http://genomescope.org/genomescope2.0/. Looking at a k-mer histogram; you should be able to see what is the coverage of the possible genomic k-mers. 

To find all the k-mer pairs, we must chose a theshold for excluding low frequencing k-mers that will be considered errors with the aim of processing mostly real genomic k-mers. That choice is not too difficult to make by looking at the k-mer spectra. 

![image](https://github.com/user-attachments/assets/29b4100b-70aa-4a45-a93b-6fdd108e2ccc)

In this example, a meaningful error threshold would be 10. As a rule of thumb, no dataset should have this threshold much below <10 unless it's a very clean sequencing run (you can separate error and genomic peaks). Furthermore, it is not the end of the world if we lose a bit of the real genomic k-mers (as long as there is enough signal). However these are just some gudances, what is sensible really depends on each individual datasets!

Now we can finally find k-mer pairs by running `smudgeplot.py hetmers`. This command will interanlly call a C-kernel optimised for the searched designed by Gene Myers. We specify `-L 10` (the coverage threshold for genomic k-mers) and use 4 threads (`-t`). The parameter `-o` just specifies the pattern for the output files

```
smudgeplot.py hetmers smudgeplot/data/Scer/FastK_Table_SRR3265401 -L 10 -t 4 -o smudgeplot/data/Scer/smudgeplot
_pairs
```

Now a `smudgeplot/data/Scer/smudgeplot_pairs_text.smu` was created. You can take a look how it looks like

```
(smudgeplot) gitpod /workspace $ head smudgeplot/data/Scer/smudgeplot_pairs_text.smu 
10      10      300
10      11      566
10      12      530
11      11      254
10      13      576
11      12      552
10      14      552
11      13      546
12      12      328
10      15      572
```

### Infering coverage and plotting smudgeplot

Plotting smudge requires knowing and infering 1n coverage. This coverage is the same 1n coverage estimated by GenomeScope and the tools must give a consistent estimate if we belive that both converged well. By default, smudgeplot will infer coverage as well as estimate sizes of all smudges. User can specify the limits for coverage inference as well as many other things (see `smudgeplot.py all -h` for all the options). For the yeast we will just run the default fit 

```
smudgeplot.py all smudgeplot/data/Scer/smudgeplot_pairs_text.smu -o smudgeplot/data/Scer/smudgeplot
```

and generate this smudgeplot

<img width="668" alt="Screenshot 2024-10-18 at 01 49 52" src="https://github.com/user-attachments/assets/1c17caf2-7848-4691-9595-ce5ba41a424b">

For more resolution on smaller smudges, look at the log version of the plot. What ploidy you think this yeast is? 

<details>
<summary><b> Answer </b></summary>
Tetraploid, specifically of `AAAB` type. Notably, this constitution does not necesarily indicate one of the haplotypes is more divergent to others because we the B k-mers can be on different haplotype for each individual k-mer pair possibly making the 4 haplotypes equidistant. We can refute a hypothesis of two and two haplotypes that are closer to each genomes as subgenomes, as those would generate prominent AABB smudge, hence the genome is quite possibly autotetraploid.  
</details>

### Run smudgeplot on lots of data

We provide you with a fine selection of 34 species saved in indexed subdirectories in `smudgeplot/data/`. For all those species we provide you with a k-mer histogram `*.hist` and with a file with all the k-mer pairs `*.smu.txt`. Both can be generated using FastK batabase build from raw sequencing reads, how to do that, look again at the previous section if you are not sure.

We oranised the datasets in a [table of species](https://docs.google.com/document/d/1Ad_FSMoSUEtG0aBEJCa5jmZvAkjxkek9qiTVxSA4fE0/edit?usp=sharing), we will use the same document to upload our results too; Pick the one that is not taken yet, add your name to the document and analyse the dataset. 

It you see three columns, it's a good sign. You can proceed to finally plot the smudgeplot. I would encourage to run `smudgeplot plot -h` to see all the options and understand what they mean, but a minimilistic command like this should do:

```
smudgeplot.py plot -t SRR926341 -o SRR926341_k31_smudgeplot SRR926341_k31_pairs_text.smu
```

How does the smudgeplot look? You shold see something like this:

![smudgeplot](https://user-images.githubusercontent.com/8181573/267332563-1b9d8bc1-6241-4ebb-a92a-32d02c7c38d1.png)

A plot with a bunch of smudges, and annotations that are overlapping the smudges. In the top right panel you see proportions of kmer pairs in the individual smudges sorted by frequency. In the bottom right corner you see the 1n coverage estimate for the dataset. This is the same 1n coverage as was infered by GenomeScope, these two numbers need to be the same for the model and smudgeplot be telling the same story. If they are substantially different, one should investigate why. In different genomes smudgeplot or genomescope are better in figuring out the coverage, and usually the diffences are in factor of 2. If your think it's smudgeplot coverage estimate that is off, rerun smudgeplot with paramter '-n' and provide a number corresponding to the 1n peak in your genomescope plot such as '-n 50'.

Once you are happy with your smudgeplot, upload it to shared docs with results.

We will discuss the results and then hear from Amjad, what is the actual biological story.

### Where to go next

- [Smudgeplot v2.5.0](https://github.com/KamilSJaron/smudgeplot/wiki) documentation: most of it applies the same for this development version (2.9.9) and there are plenty of useful things to learn in there
- original [Genomescope & Smudgeplot paper](https://www.nature.com/articles/s41467-020-14998-3): this describes a lot older version of the software (0.1.3), but the general idea applies.
- [OH-KNOW k-mer workshop learning materials](https://github.com/KamilSJaron/oh-know/wiki/Characterization-of-polyploid-genomes-using-k-mer-spectra-analysis).
- If you would be interested finding out more about research in Jaron group, visit our [website](https://www.sanger.ac.uk/group/jaron-group/).


