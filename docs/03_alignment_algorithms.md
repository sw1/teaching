
# Alignment Algorithms {#alignalg}

## Longest Common Subsequence


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

library(compiler)
find_lcs <- cmpfun(function(x,y){
  
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
  
})
```

## Global Alignment (R)


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

```r
blosum62 <- as.matrix(read.table('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/34c4a34ce49d58b595c3fb2dc77b89a5b6a9b0af/blosum62.dat'))
x <- 'ISTHISALL'
y <- 'ALIGNED'
penalty <- -5
global_alignment(x,y,blosum62,-5)
```

```
## $score
## [1] -15
## 
## $alignment
## [1] "ISTHISALL" "-AL-IGNED"
```

## Global Alignment (Python)


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
    
blosum62 = 'https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/34c4a34ce49d58b595c3fb2dc77b89a5b6a9b0af/blosum62.dat'
x = 'ISTHISALL'
y = 'ALIGNED'
penalty = -5
global_alignment(x,y,penalty,blosum62)
```
