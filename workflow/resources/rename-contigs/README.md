# rename-contigs

The reference is denoted in Ensembl notation. If the caller is using USCS notation, the contigs need to be renamed like

```
chr12
```

becomes

```
12
```

Two files are provided: One for GRCh38 as reference genome and one for GRCh37. The path to those files needs to be denoted in the user config to enable the `rename_contigs` rule.

The files were taken from [dpryan79/ChromosomeMappings](https://github.com/dpryan79/ChromosomeMappings).
