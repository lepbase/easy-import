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
git clone https://github.com/lepbase/easy-import ei # with gff-parser and ensembl-easy as submodules
cd ~/ei/ee
sudo ./install-dependencies.sh import.ini
```

### Step 1.2: Setup database connections

Again, recycle ensembl-easy

```bash
cd ~/ei/ee
./setup-databases.sh import-databases.ini
```

## 2. Core import

### Step 2.1: Fetch assembly/annotation files

```bash
cd ~/ei/core
```

### Step 2.2: create database and load sequence data

```bash
cd ~/ei/core
```

### Step 2.3: Prepare the gff file for import

Handle any exceptions

```bash
cd ~/ei/core
```

### Step 2.4: import gff from modified file

```bash
cd ~/ei/core
```

### Step 2.5: import additional annotations

- interpro
- blastp
- repeats

```bash
cd ~/ei/core
```

### Step 2.6: export files, verify and summarise

```bash
cd ~/ei/core
```

## 3. Compara import


## 4. Web site configuration

### Step 4.1: Update Ensembl code

Return to step 3 of ensembl-easy

```bash
cd ~/ei/ee
./update-ensembl-code.sh example.ini
```

### Step 4.2: Reload Ensembl site

Step 4 of ensembl-easy

```bash
cd ~/ei/ee
./reload-ensembl-site.sh my.ini
```
