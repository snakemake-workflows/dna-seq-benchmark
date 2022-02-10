# Changelog

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
