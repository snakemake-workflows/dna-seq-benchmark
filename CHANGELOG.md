# Changelog

## [1.13.0](https://github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.12.2...v1.13.0) (2025-08-01)


### Features

* add number of expected variants to report and improve the views ([#130](https://github.com/snakemake-workflows/dna-seq-benchmark/issues/130)) ([4ada0c1](https://github.com/snakemake-workflows/dna-seq-benchmark/commit/4ada0c1f16d9e4ebdee104b112e066574aa31c5a))


### Bug Fixes

* link to snakemake workflow catalog ([#129](https://github.com/snakemake-workflows/dna-seq-benchmark/issues/129)) ([1c49f2f](https://github.com/snakemake-workflows/dna-seq-benchmark/commit/1c49f2f5b9777d8b46d691efcd97d063e6afc911))
* use buffered and script based get-reads implementation ([#126](https://github.com/snakemake-workflows/dna-seq-benchmark/issues/126)) ([b9d4360](https://github.com/snakemake-workflows/dna-seq-benchmark/commit/b9d43605e736a3e919a74d811dadf93fad8c5cef))


### Performance Improvements

* add bot to automatically update github actions monthly ([#134](https://github.com/snakemake-workflows/dna-seq-benchmark/issues/134)) ([e1ac2d5](https://github.com/snakemake-workflows/dna-seq-benchmark/commit/e1ac2d5b0c67444adeff8af667f97249ef9073eb))
* update tools in conda environments ([#132](https://github.com/snakemake-workflows/dna-seq-benchmark/issues/132)) ([6d2dbc1](https://github.com/snakemake-workflows/dna-seq-benchmark/commit/6d2dbc1fa4ce97cfc3754e9050a266a9b499e3e5))

### [1.12.2](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.12.1...v1.12.2) (2025-04-15)


### Bug Fixes

* add dirs with results to rule all ([#120](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/120)) ([4dfdf6d](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/4dfdf6dd33c1e84d82f9647b24b8872e5d688f41))

### [1.12.1](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.12.0...v1.12.1) (2025-04-11)


### Bug Fixes

* add .get() for optional values in config instead of dict access ([#118](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/118)) ([94d9294](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/94d92945732fd6b30aa08ff8e267f8929098ce22))

## [1.12.0](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.11.7...v1.12.0) (2025-04-03)


### Features

* Add heatmap coloring for coverage stratification in fp fn report, hide class label. ([#112](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/112)) ([e0d3ea9](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/e0d3ea9ddf3a01ab8e02cf9a699186943e3d2c17))
* automate picking of rename file based on reference genome used for calling + setup for more than 2 reference genomes ([#113](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/113)) ([2bbb9ae](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/2bbb9ae5844982fb733fe30dbbf340cf462249b0))
* Report FP / FN negative variants in single table ([#107](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/107)) ([0464700](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/04647001dd78648309aa71d9c44152ffba4c73b6))
* add rule merge_callsets for snv/indel callsets ([#115](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/115)) ([7e5fdde](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/7e5fdde887cf8d693d2b39e2615c6ee0fea1cbac))


### Bug Fixes
* catch when FORMAT field is present in truth ([#114](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/114)) ([0e8e4be](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/0e8e4be88ba41510588060df779e14334f0b42a8))
* change vaf_field to dict ([#110](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/110)) ([e78abbb](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/e78abbbe6f999d24dedf9f934271c0335578f883))
* fix missing index as input for restrict_to_reference_contigs ([#102](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/102)) ([c55c503](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/c55c50322ea9d836702cfed4c753b87358f32ae9))
* link to truth vcf for na12878-somatic ([#116](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/116)) ([724ce2d](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/724ce2d22de2ec95000447db1deb6642e2c81e91))


### Performance Improvements

* update datavzrd wrapper ([#104](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/104)) ([db68793](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/db68793142b2fc98d444a18112f3da8802f777c0))

### [1.11.7](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.11.6...v1.11.7) (2024-09-12)


### Bug Fixes

* Only append labels in first iteration of storage loop; calculate chi2 only on non empty feature matrix. ([#96](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/96)) ([20a469a](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/20a469a6f661ac5c55ec561d8da5462710ec68d9))
* restrict calls to reference contigs ([#101](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/101)) ([0413580](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/0413580a4ce030fcf38b35298d5cf1f20664c29c))

### [1.11.6](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.11.5...v1.11.6) (2024-08-06)


### Bug Fixes

* hide sort_index in results table ([#94](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/94)) ([bb72674](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/bb72674f0530f55c5e7418c29c89955885eb2363))

### [1.11.5](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.11.4...v1.11.5) (2024-07-12)


### Performance Improvements

* update to latest datvzrd ([2ae57cc](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/2ae57ccf39e68080695b8dea409f7f44459d9e9e))

### [1.11.4](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.11.3...v1.11.4) (2024-07-08)


### Bug Fixes

* update to latest datavzrd wrappers ([ba6b039](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/ba6b039d8c8967421e4f6c836ec78f45f1e996e0))

### [1.11.3](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.11.2...v1.11.3) (2024-06-28)


### Bug Fixes

* fix f* calculation ([#89](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/89)) ([1afb621](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/1afb621e21a86210dfeb77c08b4984d70947d40e))
* Split intersection with target region and normalize calls in separate rules. ([#90](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/90)) ([795a82b](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/795a82b1536e6269cd031f6cee0903a50b44c8b6))

### [1.11.2](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.11.1...v1.11.2) (2024-06-24)


### Performance Improvements

* update to latest datavzrd (includes sorting and column hiding) ([#87](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/87)) ([182edbb](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/182edbb2cbe51e38e301e30c22feb9b09a013cc0))

### [1.11.1](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.11.0...v1.11.1) (2024-06-23)


### Bug Fixes

* fix F* calculation ([#85](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/85)) ([e684a33](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/e684a3301161680106cffddecbf459a2bb32163b))

## [1.11.0](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.10.0...v1.11.0) (2024-06-20)


### Features

* add F*-measure ([#84](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/84)) ([92f4903](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/92f4903bd9af9e875ec4a15ddb485e2313c2b513))
* add seqc2 wgs and ffpe to presets ([#82](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/82)) ([c475b0f](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/c475b0f5d263a7195946a1ec467050394696af08))
* add version of benchmark as key in presets and include in report ([#80](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/80)) ([111b4ec](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/111b4ece8413342f18c608c304e691ca8144aa93))
* offer download of datavzrd reports as excel ([#81](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/81)) ([6f5055e](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/6f5055ecdd4190782567965f77eb6bb112923a4c))

## [1.10.0](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.9.0...v1.10.0) (2024-06-11)


### Features

* improve coverage stratification for datasets with high coverage ([#75](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/75)) ([225113c](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/225113c4c5154cf1884dc211f73928159c370456))


### Bug Fixes

* change bm names to adhere to naming scheme in pm4onco ([#78](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/78)) ([20ad65e](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/20ad65e59d9b91013a4d07845885eab26cea4399))

## [1.9.0](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.8.5...v1.9.0) (2024-05-17)


### Features

* expand benchmarking datasets with IMGAG somatic validation data ([#61](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/61)) ([307e034](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/307e034ecad19d3d35a1bdc08aeb8a0aa1de815e))
* Liftover callsets from hg19/GRCh37 to hg38/GRCh38 ([#70](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/70)) ([28f96a8](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/28f96a8f53c1f03b62b00f153c6efa6ee520ff73))
* only consider PASS and no-filter variants in given callsets ([#66](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/66)) ([7c48a8c](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/7c48a8c623824a7c948376961b92ae232b95d239))
* stratification of VAFs ([#68](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/68)) ([8d2c121](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/8d2c121aa9a5ea30b1a39fe2dda0d74263aaf6fb))


### Bug Fixes

* fix failing giab testcase ([#72](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/72)) ([5f24c56](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/5f24c56646d533e621c1b6b3bda4d561016fc808))
* imgag-somatic and integrate seqc2 ([#63](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/63)) ([274d800](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/274d80085e2e6de07234c0954908c9e9a1ebe146))
* merge truthsets ([#71](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/71)) ([aa8ca0a](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/aa8ca0a79ac98cf274857e193a92d65394988b31))
* rename-contigs to include target bed intersection and add mapping files to resources ([#73](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/73)) ([65930c5](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/65930c56d6ad1abbba2036e1647f297b27c9b897))
* use latest datavzrd wrapper, integrating templating of config, thereby avoiding leaks of local storage paths into datavzrd config ([#74](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/74)) ([40a64ed](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/40a64ed1f0a54e56a580b921c69f7ecf510d56a5))

### [1.8.5](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.8.4...v1.8.5) (2023-11-13)


### Performance Improvements

* update to latest datavzrd ([#59](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/59)) ([d25bfda](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/d25bfdae77dad0fbf2b2834a0d8c2df94937bead))

### [1.8.4](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.8.3...v1.8.4) (2023-09-28)


### Bug Fixes

* avoid issue that occurs with single callset ([#57](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/57)) ([c41b1c3](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/c41b1c3cc3ad1de03ec17a20ccc8030fdb46268a))

### [1.8.3](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.8.2...v1.8.3) (2023-07-12)


### Bug Fixes

* work around issues with --all-temp and the checkpoint output ([1c8a11f](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/1c8a11f0599c8f18f5f08b0977cd47f1628644e2))

### [1.8.2](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.8.1...v1.8.2) (2023-07-11)


### Features

* limit the number of entries to be shown in the FP/FN tables ([7b2cc19](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/7b2cc193f69d90fe2ac8ed55c7cc4f80bc4148b0))


### Miscellaneous Chores

* release 1.8.2 ([faa44a6](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/faa44a6a9b5fdb06256277107341dc17bd24b760))

### [1.8.1](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.8.0...v1.8.1) (2023-07-11)


### Bug Fixes

* omit index creation in FP/FN collection ([aa090b1](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/aa090b10faf37eb52119b576574c2ad4bae16b9e))

### [1.8.1](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.8.0...v1.8.1) (2023-07-11)


### Bug Fixes

* omit index creation in FP/FN collection ([aa090b1](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/aa090b10faf37eb52119b576574c2ad4bae16b9e))

## [1.8.0](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.7.4...v1.8.0) (2023-07-11)


### Features

* chunked processing of fp/fn collection in order to save memory ([#51](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/51)) ([a36a13b](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/a36a13b96082d5caf4a42149de36f4dc980c0e49))

### [1.7.4](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.7.3...v1.7.4) (2023-07-03)


### Performance Improvements

* update to latest datavzrd ([fb0127c](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/fb0127cdfacc6aacafeacec3c11414a82205532c))

### [1.7.3](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.7.2...v1.7.3) (2023-06-21)


### Bug Fixes

* use a precision of 3 for displaying floats in datavzrd, use an FDR threshold of 0.25 in order to see more borderline cases ([237b86c](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/237b86c8bb7c4c1f744740eae75ae8f4819587d7))

### [1.7.2](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.7.1...v1.7.2) (2023-06-14)


### Performance Improvements

* update to latest datavzrd ([#46](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/46)) ([67171cf](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/67171cf6682d96487dcdb2088c5608b0d3721c26))

### [1.7.1](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.7.0...v1.7.1) (2023-05-30)


### Bug Fixes

* use 10% FDR threshold for correlation test ([b596fa3](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/b596fa3c3b2a6c35718039c424f99c4e0455fb77))
* use p-value <=0.05 for filtering chi-square results instead of FDR. While the FDR is still displayed, this better allows to check individual results without being influenced by the amount of multiple testing. ([096ffd0](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/096ffd0cfc1fb471ee037f8415d979fac202019f))

## [1.7.0](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.6.7...v1.7.0) (2023-05-25)


### Features

* add plot view of precision and recall ([28e4249](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/28e424916f0cc1b201e245a818639e3193277546))

### [1.6.7](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.6.6...v1.6.7) (2023-05-04)


### Bug Fixes

* define retries for fallible rules ([#41](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/41)) ([2f7dddf](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/2f7dddfdb894e3582bd1fb129179f20b5ab078c4))
* error out if a callset does not match the truth at all ([#43](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/43)) ([480b929](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/480b9296c84e5c71ffa2e560a9c50085836577fa))


### Performance Improvements

* update to latest datavzrd ([8022a46](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/8022a46ef81fa5365934f250cf4a5b1312388e33))

### [1.6.6](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.6.5...v1.6.6) (2023-04-27)


### Bug Fixes

* fixes for custom benchmarks and cleanup of chi2 output ([#39](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/39)) ([c8cdfb8](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/c8cdfb826f60749b9b8ab8cd218d45b991158764))

### [1.6.5](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.6.4...v1.6.5) (2023-03-20)


### Bug Fixes

* use Benjamini/Yekutieli method for FDR correction ([9391a1a](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/9391a1acbbab6a1d0f7b76060720997c36ffdaf2))

### [1.6.4](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.6.3...v1.6.4) (2023-02-12)


### Bug Fixes

* use resources to set sort threads ([ed5bca0](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/ed5bca0a49672b61a50b324128aaa81d03ffea65))

### [1.6.3](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.6.2...v1.6.3) (2023-02-12)


### Bug Fixes

* polish reports ([947a176](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/947a176fda49b64c2fcb098d07c4967a78f6d306))

### [1.6.2](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.6.1...v1.6.2) (2023-01-13)


### Bug Fixes

* correct columns in case of empty precision recall ([#33](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/33)) ([03e60da](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/03e60da54e5b55ed8036e3fa9416bb0a27817a05))
* display of false negatives ([#35](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/35)) ([b3feca6](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/b3feca6698f4279a1027f202ccf97a96992410d8))

### [1.6.1](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.6.0...v1.6.1) (2023-01-10)


### Bug Fixes

* sort stratifications ([efa4b7b](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/efa4b7bc3f168c9b02d3e93bf49f3652eeb2445d))


### Performance Improvements

* latest datavzrd ([ff8687c](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/ff8687c0c60ac8fe88ec8d6bd1a0de6267846908))
* update to latest datavzrd ([645bfb7](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/645bfb7e138f33ade64be9375db34e1a6f8ae2a9))

## [1.6.0](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.5.1...v1.6.0) (2023-01-10)


### Features

* replace hap.py with direct use of vcfeval ([#30](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/30)) ([75113fa](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/75113fafc8aaf797cdd05f3652c327ef0b301e26))

### [1.5.1](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.5.0...v1.5.1) (2022-11-11)


### Bug Fixes

* add common to source files ([f4c9d0b](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/f4c9d0bb646bb689c2bc34a8a357b95807f2408d))
* drop empty rows ([#29](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/29)) ([87b85dc](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/87b85dc62ca281faff5118c88d035bcaf26949af))
* Use common module as input for the correct rules. This fixes the deployment of this workflow as a module. ([7a00a2f](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/7a00a2fc7b722a66e6670f220e035ee46133a7bf))

## [1.5.0](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.4.1...v1.5.0) (2022-10-27)


### Features

* improved precision recall reporting ([#25](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/25)) ([5f5d18f](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/5f5d18f3bffa946174074e121d892d374b6ce51e))
* improved, interactive reporting (FP/FN and precision/recall) ([#27](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/27)) ([df6e430](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/df6e4308f259bf7c79a6e4960df9c2c7bd8de1c0))


### Performance Improvements

* upgrade datavzrd ([d844bce](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/d844bce31fea29821bfd86ba694926ca6428be77))

### [1.4.1](https://www.github.com/snakemake-workflows/dna-seq-benchmark/compare/v1.4.0...v1.4.1) (2022-08-25)


### Bug Fixes

* improved genome and benchmark parsing and collection ([#22](https://www.github.com/snakemake-workflows/dna-seq-benchmark/issues/22)) ([07934ad](https://www.github.com/snakemake-workflows/dna-seq-benchmark/commit/07934ad82b0276574bab603fa046f16ddaf97591))

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
