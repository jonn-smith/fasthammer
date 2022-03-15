# Fast Hammer
```
                                      )\   /|
                                   .-/'-|_/ |
                __            __,-' (   / \/          
            .-'"  "'-..__,-'""          -o.`-._   
           /                                   '/
   *--._ ./                                 _.-- 
         |   [===="                      _.-' 
         :     ||                    .-/   
          \    ||                )_ /
           \   ||           _)   / \(
             `.   /-.___.---'(  /   \\
              (  /   \\       \(     L\
               \(     L\       \\
                \\              \\
                 L\              L\
  [nabis]
```

Fast hammer is a tool to calculate all pairwise
hamming distances between words on a given list.

It will produce a distribution of the non-redundant
hamming distances between all given words, as 
well as a file containing the minimum hamming 
distance between a given word and any other word
on the list.

It is both multi-threaded and vectorized, so it's 
reasonably fast (though I'm sure it could be 
better).

## Compiling / Running:
```
make
./fastHammer testFile.txt 1 
```

Example:

```
time ./fastHammer testFile.txt 1
Reading in barcodes from testFile.txt...
Num barcodes read in: 104
Computing hamming distances in 1 threads for barcodes of length: 16
Progress %:
        10        20        30        40        50        60        70        80        90         100
........................................................................................................
# Hamming distance distribution:

0       20
1       32
2       0
3       0
4       0
5       0
6       0
7       0
8       0
9       0
10      0
11      0
12      0
13      0
14      0
15      64
16      5240

Writing per-barcode min hamming distance to: minimum_hamming_distances.tsv
DONE


real    0m0.010s
user    0m0.006s
sys     0m0.001s
```

