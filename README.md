# Easy Import

Code to make it easy to import heterogeneous data into an [EnsEMBL](http://ensembl.org)
  database.

## 1. Server setup

### Step 1.1: Install dependencies

Recycling ensembl-easy for now - may be too many dependencies?

```bash
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install git
cd ~
git clone --recursive https://github.com/lepbase/easy-import ei # with gff-parser and ensembl-easy as submodules
cd ~/ei/ee
sudo ./install-dependencies.sh ../conf/setup.ini
```

### Step 1.2: Setup database connections

Again, recycle ensembl-easy

```bash
cd ~/ei/ee
./setup-databases.sh ../conf/setup-db.ini
```

### Step 1.3: git clone any essential ensembl repositories

Recycle ensembl-easy (could be submodules?)

```bash
cd ~/ei/ee
./update-ensembl-code.sh ../conf/setup.ini
```


## 2. Core import

### Step 2.1: (optional) Fetch/summarise assembly/annotation files

```bash
mkdir ~/import
cd ~/import
perl ../ei/core/00_summarise_files.pl ../ei/conf/test.ini
```

```bash
summary/GFF.Obru_genes.gff.sorted.attribute_counts.txt (END)
CDS     ID      16912
CDS     Parent  16912
exon    ID      79932
exon    Parent  79932
five_prime_UTR  ID      759
five_prime_UTR  Parent  759
gene    Alias   16912
gene    ID      16912
gene    Name    16912
mRNA    Alias   16912
mRNA    ID      16912
mRNA    Name    16912
mRNA    Parent  16912
three_prime_UTR ID      35
three_prime_UTR Parent  35
summary/GFF.Obru_genes.gff.sorted.attribute_counts.txt (END)
```

### Step 2.2: create database and load sequence data

```bash
cd ~/import
perl ../ei/core/10_import_sequences.pl ../ei/conf/test.ini
perl ../ei/core/12_import_sequence_synonyms.pl ../ei/conf/test.ini
```

### Step 2.3: Prepare the gff file for import

```bash
cd ~/ei/core
perl ../ei/core/20_prepare_gff.pl ../ei/conf/test.ini
```

TODO - Handle any exceptions

### Step 2.4: import gff from modified file

```bash
cd ~/ei/core
perl ../ei/core/30_import_gene_models.pl ../ei/conf/test.ini
```

### Optional: import additional annotations

TODO - test these:

```bash
cd ~/ei/core
perl ../ei/core/60_import_blastp.pl ../ei/conf/test*.ini
perl ../ei/core/62_import_repeatmasker.pl ../ei/conf/test*.ini
perl ../ei/core/64_import_interproscan.pl ../ei/conf/test*.ini
```

### Optional: export files

```bash
cd ~/ei/core
```

### Optional: verify import

```bash
cd ~/ei/core
```

### Optional: generate search index

```bash
cd ~/ei/core
```

### Optional: generate files for web

```bash
cd ~/ei/core
```


## 3. Compara import


## 4. Web site configuration

### Step 4.1: Update Ensembl webcode

Return to step 3 of ensembl-easy

```bash
cd ~/ei/ee
./update-ensembl-code.sh example.ini
```

### Step 4.2: Reload Ensembl website

Step 4 of ensembl-easy

```bash
cd ~/ei/ee
./reload-ensembl-site.sh my.ini
```
