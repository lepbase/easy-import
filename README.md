# Easy Import

Code to make it easy to import heterogeneous data into an [EnsEMBL](http://ensembl.org)
  database.

## 1. Server setup

### Step 1.1: Install dependencies

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

```bash
cd ~/ei/ee
./setup-databases.sh ../conf/setup-db.ini
```

### Step 1.3: git clone any essential ensembl repositories

```bash
cd ~/ei/ee
./update-ensembl-code.sh ../conf/setup.ini
```


## 2. Core import

``example.ini`` will install a new core database for *Operophtera brumata*

### Step 2.1: (optional) Fetch/summarise assembly/annotation files

```bash
mkdir ~/import
cd ~/import
perl ../ei/core/summarise_files.pl ../ei/conf/example.ini
```

### Step 2.2: create database and load sequence data

```bash
cd ~/import
perl ../ei/core/import_sequences.pl ../ei/conf/example.ini
perl ../ei/core/import_sequence_synonyms.pl ../ei/conf/example.ini
```

### Step 2.3: Prepare the gff file for import

```bash
cd ~/import
perl ../ei/core/prepare_gff.pl ../ei/conf/example.ini
```

### Step 2.4: import gff from modified file

```bash
cd ~/import
perl ../ei/core/import_gene_models.pl ../ei/conf/example.ini
```

### Optional: import additional annotations

```bash
cd ~/import
perl ../ei/core/import_blastp.pl ../ei/conf/test.ini ../ei/conf/example-extra.ini
perl ../ei/core/import_repeatmasker.pl ../ei/conf/test.ini ../ei/conf/example-extra.ini
perl ../ei/core/import_interproscan.pl ../ei/conf/test.ini ../ei/conf/example-extra.ini
```

### Optional: export files

```bash

```

### Optional: verify import

```bash

```

### Optional: generate search index

```bash
cd ~/import
perl ../ei/core/generate_file_stats.pl ../ei/conf/example.ini
```

### Optional: generate files for web

```bash
cd ~/import
perl ../ei/core/generate_file_stats.pl ../ei/conf/example.ini
```


## 3. Compara import


## 4. Web site configuration

### Step 4.1: Update Ensembl webcode

edit ``setup.ini`` to add ``operophtera_brumata_v1_core_31_84_1`` to ``[DATA_SOURCE] SPECIES_DBS``

```bash
cd ~/ei/ee
./update-ensembl-code.sh ../conf/setup.ini
```

### Step 4.2: Reload Ensembl website

```bash
cd ~/ei/ee
./reload-ensembl-site.sh ../conf/setup.ini
```
