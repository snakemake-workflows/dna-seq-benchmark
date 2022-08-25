# Changelog

## [1.4.0](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.3.6...v1.4.0) (2022-08-25)


### Features

* collect FP and FN calls and compare them across callsets ([#18](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/18)) ([cca4571](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/cca45715f6a0f373c109486fd44c013ff66efbe8))


### Performance Improvements

* update snakemake wrappers to v1.7.2 ([#19](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/19)) ([65ffae3](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/65ffae30204e3337a807b254a31a91451a27550b))

### [1.3.6](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.3.5...v1.3.6) (2022-07-20)


### Bug Fixes

* fix genotype-ignoring recall computation ([a83f00d](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/a83f00dc24132bee3dcf1452ceea20794118c0f8))

### [1.3.5](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.3.4...v1.3.5) (2022-06-08)


### Bug Fixes

* fix samtools sort argument when retrieving reads ([05abcae](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/05abcae07941ed0f00f5d51adac8b9a2428e2cd6))

### [1.3.4](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.3.3...v1.3.4) (2022-05-09)


### Bug Fixes

* fix bwa index path ([56c1a29](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/56c1a29efce6da048bdc97cc21a6833d8ec58d9f))

### [1.3.3](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.3.2...v1.3.3) (2022-05-05)


### Bug Fixes

* fix param ([7e82d84](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/7e82d843debdc960b91b7af1139ccebc93334c5e))

### [1.3.2](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.3.1...v1.3.2) (2022-05-05)


### Bug Fixes

* fix bwa wrapper inputs ([63e92f2](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/63e92f23cb06b9dac6a2001877300496d3a28347))

### [1.3.1](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.3.0...v1.3.1) (2022-05-05)


### Bug Fixes

* use latest bwa wrapper to avoid issues with mamba pulling an old r-base package ([9a53ad0](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/9a53ad0a877f47562eaf043f56dfad7b43c1abc2))

## [1.3.0](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.2.0...v1.3.0) (2022-05-05)


### Features

* report columns ([#9](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/9)) ([dd8f776](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/dd8f776a3037c2c15359d6ea7a4e891ea0880dd7))


### Bug Fixes

* define target bed as local input file if not a URL ([#8](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/8)) ([a2c9844](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/a2c9844659ec8e1f625fe2966a4a706bd04508c8))

## [1.2.0](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.1.0...v1.2.0) (2022-03-15)


### Features

* Generalize towards being able to use other genomes than NA12878 ([#6](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/6)) ([741d8df](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/741d8df8bca1753d830835a3aa07b966678fa297))

## [1.1.0](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.0.1...v1.1.0) (2022-02-10)


### Features

* generalization to allow multiple callsets and arbitrary benchmark datasets ([#4](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/4)) ([52cef3a](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/52cef3abf2bb13321f2b7deda5da9cf5b7c06dd6))

### [1.0.1](https://www.github.com/snakemake-workflows/benchmark-giab/compare/v1.0.0...v1.0.1) (2021-12-20)


### Bug Fixes

* allow missing variant-calls key ([1c9c853](https://www.github.com/snakemake-workflows/benchmark-giab/commit/1c9c85349f5e2530fe71fa315a3dd7e2b74ae441))

## 1.0.0 (2021-12-20)


### Features

* allow to limit workflow to given regions ([6e3792a](https://www.github.com/snakemake-workflows/benchmark-giab/commit/6e3792a31ba6abbef82ab2b782a3076e807d25b2))
* allow to specify custom reads instead of downloading public data in the workflow ([cee1a05](https://www.github.com/snakemake-workflows/benchmark-giab/commit/cee1a05c78402953d0f7de50f9efde0ff2558bcd))
* plot numbers in precision-recall plot ([6899ee1](https://www.github.com/snakemake-workflows/benchmark-giab/commit/6899ee1c489a631d3ae400135db12ea2f5464a72))
* support multiple callsets ([65268e4](https://www.github.com/snakemake-workflows/benchmark-giab/commit/65268e4c8c7de6451b4c4901c405458aca42cb09))


### Bug Fixes

* custom read handling ([04057c1](https://www.github.com/snakemake-workflows/benchmark-giab/commit/04057c181a9ecbf0884f995c9eb32058c188e9ae))
* ensure that stratified regions do not contain duplicate entries ([5dbafc9](https://www.github.com/snakemake-workflows/benchmark-giab/commit/5dbafc952cdf864b9fc05c5f12ef16ffa7aa35af))
* fixed error with uninitialized variable ([1209946](https://www.github.com/snakemake-workflows/benchmark-giab/commit/1209946947f0082e366fe928a3c83bd0a8f944a7))
* handling of custom reads ([1bad2fc](https://www.github.com/snakemake-workflows/benchmark-giab/commit/1bad2fca185532e75e3c5c0e56f65e422746ead6))
* merge overlapping regions ([5e4c928](https://www.github.com/snakemake-workflows/benchmark-giab/commit/5e4c92863d940704570885756a8a9467f26eada6))
* mosdepth output path ([d6af6a3](https://www.github.com/snakemake-workflows/benchmark-giab/commit/d6af6a31dbb043cd62d90ff0e0ad30b96b20e191))
* removed superfluous argument ([ec259d1](https://www.github.com/snakemake-workflows/benchmark-giab/commit/ec259d16edda6a3a110a359531612e1df3d040cb))
* sort regions before merging them ([9a3c8af](https://www.github.com/snakemake-workflows/benchmark-giab/commit/9a3c8afd8a4e9506388ab2014808dfd5af4920ee))
* typo ([283e05f](https://www.github.com/snakemake-workflows/benchmark-giab/commit/283e05ff28ae89521673e0d3a9a463921c351976))
* typo ([68ae86b](https://www.github.com/snakemake-workflows/benchmark-giab/commit/68ae86bf23eb6c36bc08cb361af8e3ddb837e64e))
* typo in rule all ([894775f](https://www.github.com/snakemake-workflows/benchmark-giab/commit/894775fb51b8321dbeb72dc267d711530dfcc15b))
