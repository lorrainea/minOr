minOr
===

<b>Installation</b>: To install and compile minOr follow the instructions within the INSTALL file.

<b>INPUT</b>: A file containing a single text and a file with a list of orderings based on k.

<b>OUTPUT</b>: The number of minimizers identified.


```
Usage: 

./minOr <text> <ordering> <w> <k> <mode>

<text> - Name of input text file
<ordering> - Name of the input ordering file
<w> - Number of windows
<k> - Length of k-mer
<mode> - Mode to be used (s for standard, c for miniception, m for mod-sampling and r for rotational)
```

<b>Example</b>
```
 $ ./minOr ./data/ecoli ./orderings/standard/ecoli_standard_256_24_min 256 24 s
 $ ./minOr ./data/ecoli ./orderings/standard/ecoli_standard_256_24_max 256 24 s
```

