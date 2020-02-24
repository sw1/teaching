
# Alignment {#align}

## Longest Common Subsequence

We're going to tackle alignment programatically in three steps: (1) finding the longest common subsequence, (2) performing a global alignment, and finally (3) performing a local alignment. These steps build upon one another, so the order should seem natural.

A *subsequence* is simply an ordered set; **it need not be consecutive**. For example, if we had the sequence ABCDEFG, then ABC, ACF, and DFG would all be subsequences, whereas AGD would not. We aim to find the largest subsequence that two nucleotide sequences share. This boils down to essentially an alignment problem only involving insertions and deletions.

The algorithm is as follows:


```r
DEFINE LCS:

Given sequences x and y

1 Create empty an matrix S for scores and an empty matrix B for the backtracking path
2 For i in 2 to length x
..3 For j in 2 to length y
....4 Obtain the score of the upper S[i,j-1], upper-left diagonal S[i-1,j-1], 
....  and left S[i-1,j] cells.
...... if the current nucleotides in x and y are the same, such that x[i-1] == y[i-1], 
...... set s[i,j] to the maximum score among S[i,j-1], S[i-1,j-1] + 1, and S[i-1,j], 
...... where "+1" is the bonus for a match; 
...... if x[i-1] != y[i-1], then set s[i,j] to the maximum score among S[i,j-1] and S[i-1,j]
....6 Record the position of the maximum score (upper, upper-left, or left)
..End loop
End loop

7 Perform backtrack
```

We're going to create a function that computes the LCS based on the pseudocode above. It'll take two sequences, x and y. We'll call the function find_lcs, and it'll wrap our code like so:


```r
find_lcs <- function(x,y){
  
  ## code will go here
  
}
```

We'll start with the following two sequences:


```r
x <- 'AGCAGACACGTGAT'
y <- 'ATCACCGGTAT'
```

Now, the actual algorithm can be written a bunch of ways, but we're first going to split the sequence strings into character vectors, for example:


```r
unlist(strsplit(x,''))
```

```
##  [1] "A" "G" "C" "A" "G" "A" "C" "A" "C" "G" "T" "G" "A" "T"
```

so,


```r
x <- unlist(strsplit(x,''))
y <- unlist(strsplit(y,''))
```

We also will create two variables, m and n, that will store the length of our sequences:


```r
m <- length(x)
n <- length(y)
```

Lastly, we need to preinitialize our score and backtrack matrices. Recall that each is padded by an additional row and column that is filled in with either zeros or some type of penalty:


```r
s <- matrix(0,length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
b <- matrix('',length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
```

This gives us matrices that look like 


```
##     A T C A C C G G T A T
##   0 0 0 0 0 0 0 0 0 0 0 0
## A 0 0 0 0 0 0 0 0 0 0 0 0
## G 0 0 0 0 0 0 0 0 0 0 0 0
## C 0 0 0 0 0 0 0 0 0 0 0 0
## A 0 0 0 0 0 0 0 0 0 0 0 0
## G 0 0 0 0 0 0 0 0 0 0 0 0
## A 0 0 0 0 0 0 0 0 0 0 0 0
## C 0 0 0 0 0 0 0 0 0 0 0 0
## A 0 0 0 0 0 0 0 0 0 0 0 0
## C 0 0 0 0 0 0 0 0 0 0 0 0
## G 0 0 0 0 0 0 0 0 0 0 0 0
## T 0 0 0 0 0 0 0 0 0 0 0 0
## G 0 0 0 0 0 0 0 0 0 0 0 0
## A 0 0 0 0 0 0 0 0 0 0 0 0
## T 0 0 0 0 0 0 0 0 0 0 0 0
```

```
##      A  T  C  A  C  C  G  G  T  A  T 
##   "" "" "" "" "" "" "" "" "" "" "" ""
## A "" "" "" "" "" "" "" "" "" "" "" ""
## G "" "" "" "" "" "" "" "" "" "" "" ""
## C "" "" "" "" "" "" "" "" "" "" "" ""
## A "" "" "" "" "" "" "" "" "" "" "" ""
## G "" "" "" "" "" "" "" "" "" "" "" ""
## A "" "" "" "" "" "" "" "" "" "" "" ""
## C "" "" "" "" "" "" "" "" "" "" "" ""
## A "" "" "" "" "" "" "" "" "" "" "" ""
## C "" "" "" "" "" "" "" "" "" "" "" ""
## G "" "" "" "" "" "" "" "" "" "" "" ""
## T "" "" "" "" "" "" "" "" "" "" "" ""
## G "" "" "" "" "" "" "" "" "" "" "" ""
## A "" "" "" "" "" "" "" "" "" "" "" ""
## T "" "" "" "" "" "" "" "" "" "" "" ""
```

Putting it all together, our function should so far look like:


```r
find_lcs <- function(x,y){
  
  x <- unlist(strsplit(x,''))
  y <- unlist(strsplit(y,''))
  
  m <- length(x)
  n <- length(y)
  
  s <- matrix(0,length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
  b <- matrix('',length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
  
  ## more code will go here

}
```

Now, we have to loop through these matrices just like we would if we were doing this problem on paper. We start at the upper left corner and work our way down to the bottom right. As we move down, we pay close attention to the upper, left, and upper-left diagonal cells relative to our current position. Starting with the loops, our function should look like


```r
find_lcs <- function(x,y){
  
  x <- unlist(strsplit(x,''))
  y <- unlist(strsplit(y,''))
  
  m <- length(x)
  n <- length(y)
  
  s <- matrix(0,length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
  b <- matrix('',length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
  
  for (i in seq(2,m+1)){
    for (j in seq(2,n+1)){
  
      ## more code will go here
          
    }
  }
}
```

Pay careful attention to the indexes we're iterating over; this is critical. Position i=2, j=2, means different things with respect to the score and backtracking matrices compared to the actual sequences. This position (2,2), in the matrices, represents letters A and A for the sequences, respectively. But position 2 in either sequence, so x[2] and y[2], would actually index the second letter in either (C and T). Consequently, we can iterate with respect to the matrices' indexes and adjust when we need to index the sequences or vice versa. Here, we'll iterate with respect to the matrices and adjust the indexes for x and y, such that s[i,] corresponds to x[i-1] and s[,j] corresponds to y[j-1]:


```r
i <- 5; j <- 5

rownames(s)[i] == x[i]
```

```
## [1] FALSE
```

```r
rownames(s)[i] == x[i-1]
```

```
## [1] TRUE
```

```r
rownames(s)[j] == x[j]
```

```
## [1] FALSE
```

```r
rownames(s)[j] == x[j-1]
```

```
## [1] TRUE
```

Next, we'll add our scoring criteria, which pertains to the left, upper, and upper-left diagonal scores relative to our current position. Thus, if we are at position i, j, then we care about scores s[i-1,j], s[i,j-1], and s[i-1,j-1], respectively. We also want to add our match bonus, in which we'll add 1 to the score of the upper-left diagonal, assuming there is an actual match in our sequence at this position (which we'll get to next).


```r
find_lcs <- function(x,y){
  
  x <- unlist(strsplit(x,''))
  y <- unlist(strsplit(y,''))
  
  m <- length(x)
  n <- length(y)
  
  s <- matrix(0,length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
  b <- matrix('',length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
  
  for (i in seq(2,m+1)){
    for (j in seq(2,n+1)){
  
      s_up <- s[i-1,j]
      s_left <- s[i,j-1]
      s_diag <- s[i-1,j-1] + 1
      
      ## more code will go here
          
    }
  }
}
```

Now, we'll deal with awarding the correct score for our current position. Recall the piecewise function is as follows,


```r
....5 if the current nucleotides in x and y are the same, such that x[i-1] == y[i-1], 
....  set s[i,j] to the maximum score among S[i,j-1], S[i-1,j-1] + 1, and S[i-1,j], 
....  where "+1" is the bonus for a match; 
....  if x[i-1] != y[i-1], then set s[i,j] to the maximum score among S[i,j-1] and S[i-1,j]
```

Thus, we'll add the following to our function, which will create a vector that stores our 2 or 3 scores, depending on whether there is a match:


```r
find_lcs <- function(x,y){
  
  x <- unlist(strsplit(x,''))
  y <- unlist(strsplit(y,''))
  
  m <- length(x)
  n <- length(y)
  
  s <- matrix(0,length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
  b <- matrix('',length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
  
  for (i in seq(2,m+1)){
    for (j in seq(2,n+1)){
  
      s_up <- s[i-1,j]
      s_left <- s[i,j-1]
      s_diag <- s[i-1,j-1] + 1
      
      if (x[i-1]==y[j-1]) scores <- c(s_up,s_left,s_diag) else scores <- c(s_up,s_left)
      
      score_update <- max(scores)
      
      ## more code will go here
          
    }
  }
}
```

We also need to know which of the 3 cells contributed to our current position's score, so we'll record the index of the max score in our scores vector. We'll call it backtrack update (which will soon be obvious).


```r
find_lcs <- function(x,y){
  
  x <- unlist(strsplit(x,''))
  y <- unlist(strsplit(y,''))
  
  m <- length(x)
  n <- length(y)
  
  s <- matrix(0,length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
  b <- matrix('',length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
  
  for (i in seq(2,m+1)){
    for (j in seq(2,n+1)){
  
      s_up <- s[i-1,j]
      s_left <- s[i,j-1]
      s_diag <- s[i-1,j-1] + 1
      
      if (x[i-1]==y[j-1]) scores <- c(s_up,s_left,s_diag) else scores <- c(s_up,s_left)
      
      backtrack_update <- which.max(scores)
      score_update <- max(scores)
      
      ## more code will go here
          
    }
  }
}
```

Our backtrack_update variable records the position in the scores vector that had the max value, whereas the score_update variable contains the actual max score. We can use the position of the max value as a way of indexing a dictionary that can contain arrows that point to the cell that contributed to the current positions score. We'll name this dictionary backtrack_key and place it outside of the loop:


```r
find_lcs <- function(x,y){
  
  x <- unlist(strsplit(x,''))
  y <- unlist(strsplit(y,''))
  
  m <- length(x)
  n <- length(y)
  
  backtrack_key <- c('|','--','\\')
  
  s <- matrix(0,length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
  b <- matrix('',length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
  
  for (i in seq(2,m+1)){
    for (j in seq(2,n+1)){
  
      s_up <- s[i-1,j]
      s_left <- s[i,j-1]
      s_diag <- s[i-1,j-1] + 1
      
      if (x[i-1]==y[j-1]) scores <- c(s_up,s_left,s_diag) else scores <- c(s_up,s_left)
      
      backtrack_update <- which.max(scores)
      score_update <- max(scores)
      
      ## more code will go here
          
    }
  }
}
```

Finally, we can update our s and b matrices:


```r
find_lcs <- function(x,y){
  
  x <- unlist(strsplit(x,''))
  y <- unlist(strsplit(y,''))
  
  m <- length(x)
  n <- length(y)
  
  backtrack_key <- c('|','--','\\')
  
  s <- matrix(0,length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
  b <- matrix('',length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
  
  for (i in seq(2,m+1)){
    for (j in seq(2,n+1)){
  
      s_up <- s[i-1,j]
      s_left <- s[i,j-1]
      s_diag <- s[i-1,j-1] + 1
      
      if (x[i-1]==y[j-1]) scores <- c(s_up,s_left,s_diag) else scores <- c(s_up,s_left)
      
      backtrack_update <- which.max(scores)
      score_update <- max(scores)
      
      s[i,j] <- score_update
      b[i,j] <- backtrack_key[backtrack_update]
          
    }
  }
}
```

That's the function. It will calculate both the path history and the scores for all nucleotides. Now the question is "what about the backtracking"? We'll deal with that for global alignment, but for now, we'll use the following recursive function to handle it. We also have to update an option to allow for very deep recursions (if that makes no sense, don't worry about it). We'll finish the function up with a list that contains all of the information we'd like returned.


```r
backtrack_lcs <- function(b,x,m,n){
  if (m==0 | n==0) return(NULL)
  if (b[m+1,n+1] == '\\'){
    return(c(x[m],backtrack_lcs(b,x,m-1,n-1)))
  }else if(b[m+1,n+1] == '|'){
    backtrack_lcs(b,x,m-1,n)
  }else{
    backtrack_lcs(b,x,m,n-1)
  }
}

find_lcs <- function(x,y){
  
  options(expressions=10000)
  
  x <- unlist(strsplit(x,''))
  y <- unlist(strsplit(y,''))
  
  m <- length(x)
  n <- length(y)
  
  backtrack_key <- c('|','--','\\')
  s <- matrix(0,length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
  b <- matrix('',length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
  
  for (i in seq(2,m+1)){
    for (j in seq(2,n+1)){
      
      s_up <- s[i-1,j]
      s_left <- s[i,j-1]
      s_diag <- s[i-1,j-1] + 1
      
      if (x[i-1]==y[j-1]) scores <- c(s_up,s_left,s_diag) else scores <- c(s_up,s_left)
      
      backtrack_update <- which.max(scores)
      score_update <- max(scores)
      
      s[i,j] <- score_update
      b[i,j] <- backtrack_key[backtrack_update]
      
    }
  }
  
  lcs <- backtrack_lcs(b,x,m,n)
  
  return(list(lcs=paste0(rev(lcs),collapse=''),
              length=s[length(s)],
              score=s,
              backtrack=b))
  
}
```

We can now test it out:


```r
find_lcs('AGCAGACACGTGAT','ATCACCGGTAT')
```

```
## $lcs
## [1] "ACACCGTAT"
## 
## $length
## [1] 9
## 
## $score
##     A T C A C C G G T A T
##   0 0 0 0 0 0 0 0 0 0 0 0
## A 0 1 1 1 1 1 1 1 1 1 1 1
## G 0 1 1 1 1 1 1 2 2 2 2 2
## C 0 1 1 2 2 2 2 2 2 2 2 2
## A 0 1 1 2 3 3 3 3 3 3 3 3
## G 0 1 1 2 3 3 3 4 4 4 4 4
## A 0 1 1 2 3 3 3 4 4 4 5 5
## C 0 1 1 2 3 4 4 4 4 4 5 5
## A 0 1 1 2 3 4 4 4 4 4 5 5
## C 0 1 1 2 3 4 5 5 5 5 5 5
## G 0 1 1 2 3 4 5 6 6 6 6 6
## T 0 1 2 2 3 4 5 6 6 7 7 7
## G 0 1 2 2 3 4 5 6 7 7 7 7
## A 0 1 2 2 3 4 5 6 7 7 8 8
## T 0 1 2 2 3 4 5 6 7 8 8 9
## 
## $backtrack
##      A    T    C    A    C    C    G    G    T    A    T   
##   "" ""   ""   ""   ""   ""   ""   ""   ""   ""   ""   ""  
## A "" "\\" "--" "--" "--" "--" "--" "--" "--" "--" "--" "--"
## G "" "|"  "|"  "|"  "|"  "|"  "|"  "\\" "--" "--" "--" "--"
## C "" "|"  "|"  "\\" "--" "--" "--" "|"  "|"  "|"  "|"  "|" 
## A "" "|"  "|"  "|"  "\\" "--" "--" "--" "--" "--" "--" "--"
## G "" "|"  "|"  "|"  "|"  "|"  "|"  "\\" "--" "--" "--" "--"
## A "" "|"  "|"  "|"  "|"  "|"  "|"  "|"  "|"  "|"  "\\" "--"
## C "" "|"  "|"  "|"  "|"  "\\" "--" "|"  "|"  "|"  "|"  "|" 
## A "" "|"  "|"  "|"  "|"  "|"  "|"  "|"  "|"  "|"  "|"  "|" 
## C "" "|"  "|"  "|"  "|"  "|"  "\\" "--" "--" "--" "|"  "|" 
## G "" "|"  "|"  "|"  "|"  "|"  "|"  "\\" "--" "--" "--" "--"
## T "" "|"  "\\" "|"  "|"  "|"  "|"  "|"  "|"  "\\" "--" "--"
## G "" "|"  "|"  "|"  "|"  "|"  "|"  "|"  "\\" "|"  "|"  "|" 
## A "" "|"  "|"  "|"  "|"  "|"  "|"  "|"  "|"  "|"  "\\" "--"
## T "" "|"  "|"  "|"  "|"  "|"  "|"  "|"  "|"  "\\" "|"  "\\"
```

And now with much longer sequences:


```r
xy <- readr::read_lines('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/6617b2ceae8765ff7b91f4a600fc5460b335279d/lcs_sequences')

find_lcs(xy[1],xy[2])$lcs
```

```
## [1] "CTCTAAGCCAATGGCTCAGGGTGGGTGTAGCGATCCCGCGACAGTAAGCGTCTTGGTAGTTTCATCGCAGCTTCGACCTCGGGTTATCCCGCACCACCCCTTAGACTATTAAGAGCGAATACCCAATAAGTCTTTGCTCACATCCCCAGGGTTAAAACCGCGTGAGAGCGTCTCGACACGTTCTCGTTTATAGGGCGTAATTCATTCTGAGTGTGCGCCAGACTCAATACATCGTAGTCTTCTGCAAACGACTATGAGGGATGGACGCTCGGTACTACTGATTAGGTTGACAGTGATACAAACCAGGATGTATATATCGGTAAAGCCAGTTTACCGACTGCATCGCCGAGTGAAGTTCCACTCATGAAGTAAAATCTTAGTATATGTTATGCCCCCCCTTTTTTATCGAAAGTAGACTAGGTCACTTGTTGTACCGAGCCGCGCTGACATATTAACGCACCACCGACCGTTTTAGTCCTGATTGTTGGGGCATATCAGTGAATTGCGTGAGAAGTGGTCAGGTGGAGCCAGCCGAGAATCGGCTCTAGCCTCAAAACCTACTTCCACGTTCTGCCCCACTGCGAGAGAGAACTCCGGAGTCCTTCGTAGGACCAGGGC"
```

## Global Alignment

Now we can extend the algorithm above for a global alignment problem. Unlike the LCS problem where we only accounted for insertions and deletions, we're now going to penalize for mismatches. The score we give based on a particular nucleotide-nucleotide or amino acid-amino acid pair will be defined by a scoring matrix such as the amino acid scoring matrix BLOSSUM62:


```
##    A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
## A  4  0 -2 -1 -2  0 -2 -1 -1 -1 -1 -2 -1 -1 -1  1  0  0 -3 -2
## C  0  9 -3 -4 -2 -3 -3 -1 -3 -1 -1 -3 -3 -3 -3 -1 -1 -1 -2 -2
## D -2 -3  6  2 -3 -1 -1 -3 -1 -4 -3  1 -1  0 -2  0 -1 -3 -4 -3
## E -1 -4  2  5 -3 -2  0 -3  1 -3 -2  0 -1  2  0  0 -1 -2 -3 -2
## F -2 -2 -3 -3  6 -3 -1  0 -3  0  0 -3 -4 -3 -3 -2 -2 -1  1  3
## G  0 -3 -1 -2 -3  6 -2 -4 -2 -4 -3  0 -2 -2 -2  0 -2 -3 -2 -3
## H -2 -3 -1  0 -1 -2  8 -3 -1 -3 -2  1 -2  0  0 -1 -2 -3 -2  2
## I -1 -1 -3 -3  0 -4 -3  4 -3  2  1 -3 -3 -3 -3 -2 -1  3 -3 -1
## K -1 -3 -1  1 -3 -2 -1 -3  5 -2 -1  0 -1  1  2  0 -1 -2 -3 -2
## L -1 -1 -4 -3  0 -4 -3  2 -2  4  2 -3 -3 -2 -2 -2 -1  1 -2 -1
## M -1 -1 -3 -2  0 -3 -2  1 -1  2  5 -2 -2  0 -1 -1 -1  1 -1 -1
## N -2 -3  1  0 -3  0  1 -3  0 -3 -2  6 -2  0  0  1  0 -3 -4 -2
## P -1 -3 -1 -1 -4 -2 -2 -3 -1 -3 -2 -2  7 -1 -2 -1 -1 -2 -4 -3
## Q -1 -3  0  2 -3 -2  0 -3  1 -2  0  0 -1  5  1  0 -1 -2 -2 -1
## R -1 -3 -2  0 -3 -2  0 -3  2 -2 -1  0 -2  1  5 -1 -1 -3 -3 -2
## S  1 -1  0  0 -2  0 -1 -2  0 -2 -1  1 -1  0 -1  4  1 -2 -3 -2
## T  0 -1 -1 -1 -2 -2 -2 -1 -1 -1 -1  0 -1 -1 -1  1  5  0 -2 -2
## V  0 -1 -3 -2 -1 -3 -3  3 -2  1  1 -3 -2 -2 -3 -2  0  4 -3 -1
## W -3 -2 -4 -3  1 -2 -2 -3 -3 -2 -1 -4 -4 -2 -3 -3 -2 -3 11  2
## Y -2 -2 -3 -2  3 -3  2 -1 -2 -1 -1 -2 -3 -1 -2 -2 -2 -1  2  7
```

Note that this matrix scores for mismatches *and* matches, hence we need not develop separate match and mismatch scoring criteria; we can simply use this matrix as a lookup table. Also, unlike before, we're no longer going to recursively backtrack through our algorithm; instead, we'll use a while loop to backtrack in a fashion similar to doing the problem on paper. This allows us to take full advantage of the fact we have the longest path recorded in our backtrack matrix and also permits us to shift the position of NTs or AAs in our sequences when necessary by adding '-', allowing us to create our correctly oriented alignment strings.

The algorithm is as follows:


```r
DEFINE GLOBAL_ALIGNMENT:

Given sequences x and y, penalty p, and scoring lookup table score

1 Create empty an matrix S for scores and an empty matrix B for the backtracking path
2 Preinitialize edges of S and B based on penalty
3 For i in 2 to length x
..4 For j in 2 to length y
....5 Obtain the score of the upper S[i,j-1], upper-left diagonal S[i-1,j-1], 
....  and left S[i-1,j] cells:
...... set s[i,j] to the maximum score among S[i,j-1] + p, S[i-1,j-1] + score[x[i-1],y[j-1]], 
...... and S[i-1,j] + p
....6 Record the position of the maximum score (upper, upper-left, or left).
..End loop
End loop

7 Perform backtrack
```

The BLOSSUM62 scoring matrix can be found here:


```r
blosum62 <- as.matrix(read.table('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/34c4a34ce49d58b595c3fb2dc77b89a5b6a9b0af/blosum62.dat'))
```

And our backtracking algorithm:


```r
DEFINE GLOBAL_ALIGNMENT_BACKTRACK:

Given sequences x and y, backtracking matrix b, penalty p, and scoring lookup table score

1 Create empty alignment vectors align_x and align_y and 
. preinitialze score s to 0, m to len(x), and n to len(y)
2 While m > 0 OR n > 0
..3 if b[m+1,n+1] == UP
.... update align_x
.... update align_y
.... update score
.... m -= 1
..4 if b[m+1,n+1] == LEFT
.... update align_x
.... update align_y
.... update score
.... n -= 1
..5 else
.... update align_x
.... update align_y
.... update score
.... n -= 1
.... m -= 1
..End loop
End loop
```

A sample sequence pair can be found here:


```r
xy <- readr::read_lines('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/0488087071c0b27cd7a576913d5502bdf95db6e5/global_sequences')
```

## Local Alignment

Often the alignment between two subsequences between our sequences x and y is better than their global alignment. Performing a local alignment is not much different than the algorithms above except for an important change to the recurrence (piecewise function) that dictates the score for a given position. In the global alignment recurrence, we only factor in scores from the upper-left, left, and upper cells. Now, for a local alignment, we will also factor in the sink (cell position 1,1). Thus, the pseudocode will change to:


```r
DEFINE LOCAL_ALIGNMENT:

  ...

....5 Obtain the score of the upper S[i,j-1], upper-left diagonal S[i-1,j-1], 
....  and left S[i-1,j] cells:
...... set s[i,j] to the maximum score among S[i,j-1] + p, S[i-1,j-1] + score[x[i-1],y[j-1]], 
...... S[i-1,j] + p, and 0.

  ...
```

The backtrack algorithm will also change such that if we land on a zero, we set m and n to 0 to end the backtrack. Lastly, be sure to record the position of the max score in the matrix, which should act as a starting point during backtracking.

## Local Alignment: Homework

**Edit the global alignment code with the information above to create a function that calculates the local alignment between strings x and y.** Use a penalty of -5 and a PAM250 scoring matrix:


```r
pam250 <- as.matrix(read.table('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/0b783786e3797d1bb55172e750776e26021224b0/PAM250.dat'))
```

See below for the global alignment code. Fit your algorithm to these sequences:


```r
xy <- readr::read_lines('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/6617b2ceae8765ff7b91f4a600fc5460b335279d/local_sequences')
```

Report the following:

1. The alignment in its own text file (.txt, .dat, etc.)
2. The max score
3. Your code (in a text file -- specifically, not a .doc or .docx)

It's worth commenting your code to help indicate what you were thinking. Note that a correct forward pass will result in the correct max score; a correct forward *and* backward pass is required for the correct alignment. 

A few hints:

1. You'll need to update (1) the piecewise function and (2) the score vector.
2. You no longer want to start your backtrack at the bottom right corner (the sink), so update accordingly.
3. The backtrack will require another if statement.
4. The backtrack may require some debugging to ensure that your indexes are correct.

Also, the code may take as little as 35 seconds to run (and even faster on Python) for the sequences given above. That cmdfun command wrapping the function will speed up your function; it is **not** necessary. If you're approaching 5 minutes, something is likely wrong. The correct score is around 3000.

## Global Alignment Code (R)


```r
backtrack_global <- function(x,y,b,score_matrix,penalty){
  
  m <- length(x)
  n <- length(y)
  
  score <- 0
  align_x <- NULL
  align_y <- NULL
  
  while (m > 0 | n > 0){
    if (b[m+1,n+1] == '|'){
      align_x <- c(align_x,x[m])
      align_y <- c(align_y,'-')
      score <- score + penalty
      m <- m-1
    }else if(b[m+1,n+1] == '--'){
      align_x <- c(align_x,'-')
      align_y <- c(align_y,y[n])
      score <- score + penalty
      n <- n-1
    }else{
      align_x <- c(align_x,x[m]) 
      align_y <- c(align_y,y[n])
      score <- score + score_matrix[x[m],y[n]]
      n <- n-1
      m <- m-1
    }
  }
  
  alignment <- c(paste0(rev(align_x),collapse=''),paste0(rev(align_y),collapse=''))
  
  return(list(score=score,alignment=alignment))
  
}

library(compiler)
global_alignment <- cmpfun(function(x,y,score_matrix,penalty){
  
  x <- unlist(strsplit(x,''))
  y <- unlist(strsplit(y,''))
  
  m <- length(x)
  n <- length(y)
  
  backtrack_key <- c('|','--','\\')
  
  s <- matrix(0,length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
  s[1,] <- cumsum(c(0,rep(penalty,ncol(s)-1)))
  s[,1] <- cumsum(c(0,rep(penalty,nrow(s)-1)))
  
  b <- matrix('',length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
  b[1,] <- '--'
  b[,1] <- '|'
  b[1,1] <- '\\'
  
  for (i in seq(2,m+1)){
    for (j in seq(2,n+1)){
      
      s_up <- s[i-1,j] + penalty
      s_left <- s[i,j-1] + penalty
      s_diag <- s[i-1,j-1] + score_matrix[x[i-1],y[j-1]]
      
      scores <- c(s_up,s_left,s_diag) 
      
      backtrack_update <- which.max(scores)
      score_matrix_update <- max(scores)
      
      s[i,j] <- score_matrix_update
      b[i,j] <- backtrack_key[backtrack_update]
      
    }
  }
  
  return(backtrack_global(x,y,b,score_matrix,penalty))
  
})
```

## Global Alignment Code (Python)


```python
import urllib.request

def scoring_matrix(filename):
    scoring = urllib.request.urlopen(filename).readlines()
    scoring = [i.decode("utf-8").strip('\n')  for i in scoring[1:]]
    keys = [i[0] for i in scoring]
    scoring = [i.split()[1:] for i in scoring]
    scoring_dict = {}
    for ii,i in enumerate(keys):
        scoring_dict[i] = {}
        for ji,j in enumerate(keys):
            scoring_dict[i][j] = int(scoring[ii][ji])
    return scoring_dict
    
def global_alignment(v,w,penalty,matrix_file):
    score_matrix = scoring_matrix(matrix_file)
    s = [[0]*(len(w)+1) for i in range(len(v)+1)]
    s[0] = list(range(0,penalty*len(s[0]),penalty))
    for i in range(1,len(s)):
        s[i][0] = penalty + s[i-1][0]
    path = [['||'] + ['--']*(len(w)) for i in range(len(v)+1)]
    path[0][0] = '\\'
    for i in range(1,len(v)+1):
        for j in range(1,len(w)+1):
            score = score_matrix[v[i-1]][w[j-1]]
            s[i][j] = max(s[i-1][j-1] + score, s[i][j-1] + penalty,s[i-1][j] + penalty)
            if s[i][j] == s[i-1][j] + penalty:
                path[i][j] = '||'
            if s[i][j] == s[i][j-1] + penalty:
                path[i][j] = "--"
            if s[i][j] == s[i-1][j-1] + score:
                path[i][j] = "\\"
    score = 0
    align1 = ''
    align2 = ''
    while i >= 1 or j >= 1:
        if path[i][j] == "||":
            align1 += v[i-1]
            align2 += '-'
            score += penalty
            i -= 1
        elif path[i][j] == "--":
            align1 += '-'
            align2 += w[j-1]
            score += penalty
            j -= 1
        else:
            align1 += v[i-1]
            align2 += w[j-1]
            score += score_matrix[w[j-1]][v[i-1]]
            i -= 1
            j -= 1
    align1 = align1[::-1]
    align2 = align2[::-1]
    print('\n'.join([str(score),align1,align2]))
    return [score, align1,align2]
```
