# Easy Import

Code to make it easy to import heterogeneous data into an [EnsEMBL](http://ensembl.org) database.

The instructions below will help you get an Ensembl database and website up and running in an afternoon - with four Lepidopteran genomes mirrored from Ensembl Metazoa plus a fresh import of the genome of the winter moth *Operophtera brumata* direct from publicly hosted ``.gff`` and ``.fasta`` files.

This is a sister project to [easy-mirror](https://github.com/lepbase/easy-mirror) (included as a submodule), which makes it possible to set up a mirror of any Ensembl or Ensembl Genomes (including Bacteria, Metazoa, Fungi, Plants and Protists) species in four simple steps that can be run in less than an hour on a fresh Ubuntu installation.

The latest and most complete documentation for both projects is available at [easy-import.readme.io](http://easy-import.readme.io)

## Stage 1. Server setup

### Step 1.1: Install dependencies

```bash
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install git
cd ~
git clone --recursive https://github.com/lepbase/easy-import ei
cd ~/ei/em
sudo ./install-dependencies.sh ../conf/setup.ini
```

### Step 1.2: Setup database connections

```bash
cd ~/ei/em
./setup-databases.sh ../conf/setup-db.ini
```

### Step 1.3: git clone any essential ensembl repositories

```bash
cd ~/ei/em
./update-ensembl-code.sh ../conf/setup.ini
```


## Stage 2. Core import

Using ``core-import.ini`` will install a new core database for the winter moth *Operophtera brumata*

### Step 2.1: (optional) Fetch/summarise assembly/annotation files

```bash
mkdir ~/import
cd ~/import
perl ../ei/core/summarise_files.pl ../ei/conf/core-import.ini
```

### Step 2.2: create database and load sequence data

```bash
cd ~/import
perl ../ei/core/import_sequences.pl ../ei/conf/core-import.ini
perl ../ei/core/import_sequence_synonyms.pl ../ei/conf/core-import.ini
```

### Step 2.3: Prepare the gff file for import

```bash
cd ~/import
perl ../ei/core/prepare_gff.pl ../ei/conf/core-import.ini
```

### Step 2.4: import gff from the prepared file

```bash
cd ~/import
perl ../ei/core/import_gene_models.pl ../ei/conf/core-import.ini
```

### Optional: Step 2.5: import additional annotations

```bash
cd ~/import
perl ../ei/core/import_blastp.pl ../ei/conf/example.ini ../ei/conf/core-import-extra.ini
perl ../ei/core/import_repeatmasker.pl ../ei/conf/example.ini ../ei/conf/core-import-extra.ini
perl ../ei/core/import_interproscan.pl ../ei/conf/example.ini ../ei/conf/core-import-extra.ini
perl ../ei/core/import_cegma_busco.pl ../ei/conf/example.ini ../ei/conf/core-import-extra.ini
```

### Optional: Step 2.6: export files

```bash
cd ~/import
perl ../ei/core/export_sequences.pl ../ei/conf/core-import.ini
perl ../ei/core/export_json.pl ../ei/conf/core-import.ini
```

### Optional: Step 2.7: verify import

```bash
cd ~/import
perl ../ei/core/verify_translations.pl ../ei/conf/core-import.ini
```

### Optional: Step 2.8: generate search index

```bash
cd ~/import
perl ../ei/core/index_database.pl ../ei/conf/core-import.ini
```


## Stage 3. Web site configuration

### Step 3.1: Update Ensembl webcode

edit ``setup.ini`` to add ``operophtera_brumata_v1_core_31_84_1`` to ``[DATA_SOURCE] SPECIES_DBS``

```bash
cd ~/ei/em
./update-ensembl-code.sh ../conf/setup.ini
```

### Step 3.2: Reload Ensembl website

```bash
cd ~/ei/em
./reload-ensembl-site.sh ../conf/setup.ini
```
