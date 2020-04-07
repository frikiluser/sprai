// J Comput Biol. 1997 Fall;4(3):369-83.
// ReAligner: a program for refining DNA sequence multi-alignments.
// Anson EL1, Myers EW.
//
// we modified Eric & Myers ReAligner, allowed by Myers on 10/22/2012
// written by Takamasa Imai

double MIX = 0.5;
int BandSize = 8;
int minimum_ballots = 3;
int maximum_ballots = 11;

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include "col2fqcell.h"

#define FoS 16

int MaxDepth = FoS*32;
//int MaxDepth = 512;
//int MaxDepth = 64;// for large
//int MaxFragSize = 32*4*1024*1024;// for large
int MaxFragSize = 131072;
int MaxFragNum = 1024*1024;// max # of reads
int NumIts = 0;
int NumCons = 0;

int DoneFlag = 0;

int Rows;
#define BUFSIZE 4096
char * buffer;
char * sbuf;
char * qbuf;
char * chrname;
int chr_is=0;
int comment_is=0;

int opt_fastq=0;
int opt_qvs=0;
int opt_consensus=0;
int opt_vertical=0;


base_t * buf4print;
char * buf4printSeq;
char * buf4printDepth;
char * buf4printQual;
char * base_exists;
char * buf4printComment;

double distinguishable=0.70;

typedef struct fragEl{
  base_t el;// typedef struct base_t{char base; char qv;}base_t;
//   u
//  p.n
//   d
  struct fragEl *prev;
  struct fragEl *next;
  struct fragEl *up;
  struct fragEl *down;
}fragEl;

typedef struct ElPool{
  struct ElPool * ext;
  struct fragEl * Elements;
  size_t count;
  size_t size;
}ElPool;

void* 
malloc_or_die(size_t size, const char* message)
{
  void *pval;
  pval = malloc(size); 
  if(pval == NULL){ 
    fputs("Fatal error!\n", stderr);
    fputs("cannot allocate memory: ", stderr);
    fputs(message, stderr);
    fputs("\n", stderr);
    exit(1); 
  } 
  return pval;
}
void* 
realloc_or_die(void*buf, size_t size, const char* message)
{
  void *pval;
  pval = realloc(buf, size); 
  if(pval == NULL){ 
    fputs("Fatal error!\n", stderr);
    fputs("cannot reallocate memory: ", stderr);
    fputs(message, stderr);
    fputs("\n", stderr);
    exit(1); 
  } 
  return pval;
}

ElPool*
init_ElPool(ElPool* elpool, size_t size)
{
  fragEl* FreeElList = NULL;
  elpool->ext = NULL;
  elpool->count = 0;
  elpool->size = size;
  FreeElList = (fragEl *) malloc_or_die(sizeof(fragEl)*size, "FreeElList");
  elpool->Elements = FreeElList;
  return elpool;
}
fragEl *
getEl(ElPool* pool){
  size_t el_index;
  if(pool->count >= pool->size){
    if(pool->ext == NULL){
      pool->ext = malloc_or_die(sizeof(ElPool),"el_pool_expand");
      init_ElPool(pool->ext, pool->size * 1.4);
    }
    return getEl(pool->ext);
  }
  if(pool->Elements == NULL){
    fprintf(stderr, "pool->Elements is NULL\n");
    abort();
  }
  el_index = pool->count;
  pool->Elements[el_index].el.base = -1;
  pool->Elements[el_index].el.qv = ' ';
  pool->Elements[el_index].prev = NULL;
  pool->Elements[el_index].next = NULL;
  pool->Elements[el_index].up = NULL;
  pool->Elements[el_index].down = NULL;
  return &pool->Elements[pool->count++];
}

size_t
freeElPool(ElPool* pool)
{
  size_t ext_size;
  if(pool == NULL)
    return 0;
  free(pool->Elements);
  ext_size = freeElPool(pool->ext);
  ext_size += pool->size ;
  return ext_size;
}
void
flushElPool(ElPool*pool)
{
  size_t ext_size;
  ext_size = freeElPool(pool->ext);
  pool->ext = NULL;
  if(ext_size){
    pool->size += ext_size;
    pool->Elements = realloc_or_die(pool->Elements, sizeof(fragEl)*pool->size, "expanding element pool");
  }
  pool->count = 0;
}

typedef struct colStruc{
  int colInf[6];// Information, not Infimum // colInf[.] <- the number of '.'
  int  colDepth;
  double preCalc[6];
//  p.n
//   f
  struct colStruc *prev;
  struct colStruc *next;
  fragEl frags;
}colStruc;

typedef struct{
  fragEl *ofrag;
  fragEl *lastEl;
  base_t *efrag;//encoded
  base_t *afrag;//ascii
  int len;
  colStruc *scol;
  colStruc *ecol;
  int row;
}frag;

frag    *Frags;// reads
colStruc  *FirstCol;
colStruc  *LastCol;

typedef struct colPool{
  struct colPool *ext;
  colStruc *Cols;
  int count;
  size_t size;
}colPool;

colPool*
init_colPool(colPool* pool, size_t size)
{
  pool->ext = NULL;
  pool->count = 0;
  pool->size = size;
  pool->Cols = (colStruc *) malloc_or_die(sizeof(colStruc)*size, "init_colPool");
  return pool;
}

colStruc *getCol(colPool*pool){
  colStruc *acol = NULL;
  if(pool->count >= pool->size){
    fprintf(stderr,"never come here (getCol() in myrealigner)\n");
    exit(1);
  }
  int i;
  if(pool == NULL){
    fprintf(stderr, "FreeColList is NULL\n");
    abort();
  }
  acol = pool->Cols+pool->count;
  for (i=0; i < 6; ++i)
    acol->colInf[i] = 0;
  acol->colDepth = 0;
  acol->preCalc[5] = 0.0;
  acol->next = NULL;
  acol->prev = NULL;
  acol->frags.up = &(acol->frags);
  acol->frags.down = &(acol->frags);
  acol->frags.prev = (fragEl*) &(acol);
  return &(pool->Cols[pool->count++]);
}
size_t
freeColPool(colPool* pool)
{
  size_t ext_size;
  if(pool == NULL)
    return 0;
  free(pool->Cols);
  ext_size = freeColPool(pool->ext);
  ext_size += pool->size ;
  return ext_size;
}
void
flushColPool(colPool*pool)
{
  size_t ext_size;
  ext_size = freeColPool(pool->ext);
  pool->ext = NULL;
  if(ext_size){
    pool->size += ext_size;
    pool->Cols = realloc_or_die(pool->Cols, sizeof(fragEl)*pool->size, "expanding element pool");
  }
  pool->count = 0;
}


int encode[128];
int * row;
int maxNumFrags=-1;
FILE * input;

int readFrags(ElPool*elpool,colPool*colpool){
  char *s;
  int  i, r, numFrags;
  colStruc  *curCol;
  fragEl  *elPtr, *telPtr;

  for (i=0; i < MaxDepth; ++i)
    row[i] = -1;
  //buffer[BUFSIZE-1] = '\0';
  numFrags = 0;
  elPtr = telPtr = getEl(elpool);
  elPtr->el.base = -1;
  elPtr->el.qv = ' ';
  for(i=1; i < BandSize; ++i){
    // ...epeqer- -> ...epeqer--------  
    //                        bandsize  
    //          t                    t  
    //          e                    e  
    //          l                    l  
    //          p                    p  
    //          t                    t  
    //          r                    r  
    elPtr->next = getEl(elpool);
    elPtr = elPtr->next;
    elPtr->el.base = -1;
    elPtr->el.qv = ' ';
  }

  FirstCol = curCol = getCol(colpool);

  Rows = -1;
  while((s=fgets(sbuf, 2*MaxFragSize-1, input)) != NULL){
    i = strlen(sbuf);
    if(i == 0){
      fprintf(stderr, "strange input format\n");
      //return;
    }
    else if(sbuf[i-1] == '\n'){
      sbuf[--i] = '\0';
    }
    else{
      fprintf(stderr,"Each input line must not be more than %d chars\n",BUFSIZE-1);
      abort();
    }
    if(sbuf[0] == '%'){
//      if(opt_fastq != 1){
//        printf("%s\n",sbuf);
//      }
//      else{
        if(chr_is==0){
          strcpy(chrname, &sbuf[1]);
          chr_is=1;
        }
        else{
          fprintf(stderr, "soutei gai 100\n");
          exit(1);
        }
        //printf("@%s\n",&sbuf[1]);
//      }
      continue;
    }
    else if(sbuf[0] == '#'){
      if(chr_is && comment_is == 0){
        strcpy(buf4printComment,&sbuf[1]);
        comment_is=1;
      }
      else{
        fprintf(stderr, "soutei gai 101\n");
        exit(1);
      }
      continue;
    }
    if(opt_qvs)
    {
      if(sbuf[0] == '\0'){
        buffer[0] = qbuf[0] = '\0';
      }
      else{
        int i;
        for(i=0; sbuf[i] != '\t'; ++i){
          buffer[i] = sbuf[i];
        }
        buffer[i] = '\0';
        int j;
        for(j=0,++i; sbuf[i] != '\0'; ++j,++i){
          qbuf[j] = sbuf[i];
        }
        qbuf[j] = '\0';
        if(strlen(buffer) != strlen(qbuf)){
          fprintf(stderr, "the input format is broken\n");
          fprintf(stderr, "#%s#\n",sbuf);
          fprintf(stderr, "#%s#\n",buffer);
          fprintf(stderr, "#%s#\n",qbuf);
          abort();
        }
        if((int)strlen(buffer)>MaxDepth){
          fprintf(stderr,"too thick: depth %d\n",MaxDepth);
          fprintf(stderr,"%s\n",sbuf);
          fprintf(stderr,"%s\n",buffer);
          fprintf(stderr,"%s\n",qbuf);
          abort();
        }
      }
    }
    else{
      strcpy(buffer,sbuf);
      //strcpy(qbuf,sbuf);
      int i;
      for(i=0; buffer[i] != '\0'; ++i){
        qbuf[i] = '5';// dummy value
      }
      if((int)strlen(buffer)>MaxDepth){
        fprintf(stderr,"too thick: depth %d\n",MaxDepth);
        fprintf(stderr,"%s\n",sbuf);
        fprintf(stderr,"%s\n",buffer);
        fprintf(stderr,"%s\n",qbuf);
        abort();
      }
      //fprintf(stdout, "#%s#\n",buffer);
      //fprintf(stdout, "#%s#\n",qbuf);
    }

    r = 0;
    curCol->next = getCol(colpool);
    curCol->next->prev = curCol;
    curCol = curCol->next;
    int j;
    for (j=0; buffer[j] != '\0'; ++j){
      if(buffer[j] != ' '){
        if(r > Rows)
          Rows = r;
        if(r >= MaxDepth){
          fprintf(stderr, "too deep. Change the MaxDepth value in your realigner.c\n");
          abort();
        }
        if(row[r] == -1){
          row[r] = numFrags;// id in the unit
          Frags[numFrags].scol = curCol;// read's starting column?
          Frags[numFrags].row = r;
          if(numFrags >= maxNumFrags){
            Frags[numFrags].efrag = (base_t *) malloc(MaxFragSize*sizeof(base_t));
            if(Frags[numFrags].efrag == NULL){
              fprintf(stderr, "cannot allocate memory: Frags[numFrags].efrag\n");
              abort();
            }
            Frags[numFrags].afrag = (base_t *) malloc(MaxFragSize*sizeof(base_t));
            if(Frags[numFrags].afrag == NULL){
              fprintf(stderr, "cannot allocate memory: Frags[numFrags].afrag\n");
              abort();
            }
          }
          Frags[numFrags].len = 0;
          elPtr = getEl(elpool);
          for(i=0; i < BandSize; ++i){
            // e1e2e3... -> --------e1e2e3...
            //              bandsize  b a n d s i z e
            //              t                       p
            //              p                       t
            //              t                       r
            //              r                        
            // in the loop 'for (j=0, tptr=ptr=pf->ofrag; j < BandSize; ++j)...' in reAlign
            elPtr->el.base = -1;
            elPtr->el.qv = ' ';
            elPtr->next = getEl(elpool);
            elPtr->next->prev = elPtr;
            elPtr = elPtr->next;
          }
          Frags[numFrags].ofrag = Frags[numFrags].lastEl = elPtr;
          ++numFrags;
          if (numFrags == MaxFragNum){
            fprintf(stderr, "too many frags. Change the MaxFragNum value in your realigner.c\n");
            abort();
          }
        }
        else
          elPtr = Frags[row[r]].lastEl;
        if((i = encode[(int)buffer[j]]) < 0){
          fprintf(stderr,"Illegal char in input line %d\n",r);
          exit(0);
        }
        ++curCol->colInf[i];
        ++curCol->colDepth;
        elPtr->el.base = i;
        elPtr->el.qv = qbuf[j];
        elPtr->up = &(curCol->frags);
        elPtr->down = curCol->frags.down;
        curCol->frags.down = elPtr;
        elPtr->down->up = elPtr;
        elPtr->next = Frags[row[r]].lastEl = getEl(elpool);
        Frags[row[r]].lastEl->prev = elPtr;
        if(i != 0){// not '-'
          Frags[row[r]].afrag[Frags[row[r]].len].base = buffer[j];
          Frags[row[r]].afrag[Frags[row[r]].len].qv = qbuf[j];
          Frags[row[r]].efrag[Frags[row[r]].len].base = i;
          Frags[row[r]].efrag[Frags[row[r]].len].qv = qbuf[j];
          Frags[row[r]].len++;
          if(Frags[row[r]].len >= MaxFragSize){
            fprintf(stderr, "too long frag. Change the MaxFragSize value in your realigner.c\n");
            abort();
          }
        }
      }
      else if(row[r] != -1){
        Frags[row[r]].lastEl->el.base = -1;
        Frags[row[r]].lastEl->el.qv = ' ';
        Frags[row[r]].lastEl->next = telPtr;// terminal pointer (empty)
        Frags[row[r]].lastEl = Frags[row[r]].lastEl->prev;
        Frags[row[r]].ecol = curCol->prev;
        row[r] = -1;
      }
      ++r;
    }
    while (r <= Rows){
      // padding
      if(row[r] != -1){
        Frags[row[r]].lastEl->el.base = -1;
        Frags[row[r]].lastEl->el.qv = ' ';
        Frags[row[r]].lastEl->next = telPtr;
        Frags[row[r]].lastEl = Frags[row[r]].lastEl->prev;
        Frags[row[r]].ecol = curCol->prev;
        row[r] = -1;
      }
      ++r;
    }
    if(curCol->colDepth == 0){
      //fprintf(stderr, "1 unit loaded\n");
      break;
    }
  }
  if(s == NULL)
    DoneFlag = 1;
  curCol->next = LastCol = getCol(colpool);
  LastCol->prev = curCol;
  ++Rows;
  if(Rows>=MaxDepth*2){
    fprintf(stderr, "souteigai 2 (2) chr:%s\n",chrname);
    //abort();
    exit(1);
  }
  return numFrags;// num of the reads in this unit
}

fragEl  **curEl;
base_t    **curChar;

void print_dfq(char *chr, char *seq, char *depth, char *qual, char *base_exists, char *comment){
  int i;
  int limit = strlen(seq);
  //int limit = strlen(base_exists);// 0='\0'
  int from=0;
  int to=limit-1;
  for(i=0; i<limit; ++i){
    if(base_exists[i]==1){
      from = i;
      break;
    }
  }
  for(i=limit-1; i>=0; --i){
    if(base_exists[i]==1){
      to = i;
      break;
    }
  }
  /*
  if(to>=from){
    char tmp;
    printf("@%s\n",chr);

    tmp = seq[to+1];
    seq[to+1] = '\0';
    printf("%s\n",&seq[from]);
    seq[to+1] = tmp;

    tmp = depth[to+1];
    depth[to+1] = '\0';
    printf("+\t%s\n",&depth[from]);
    depth[to+1] = tmp;

    tmp = qual[to+1];
    qual[to+1] = '\0';
    printf("%s\n",&qual[from]);
    qual[to+1] = tmp;
  }
  else{
    fprintf(stderr, "buggy %s %d %d %d\n",chr, from, to, limit);
    abort();
  }
  */
  printf("@%s\n",chr);
  printf("%s\n",seq);
  printf("+\t%s\t%s\n",depth,comment);
  //printf("+\t%s\n",depth);
  printf("%s\n",qual);
  for(i=from+1; i<to; ++i){
    if(base_exists[i]!=1){
      fprintf(stderr, "WARNING: base does not exist: %d %s\n", i+1, chr);
    }
  }
  return;
}

void printAlign(int numFrags){
  int    row, i;
  colStruc  *col, *ecol;
  fragEl  *ptr;

  col = FirstCol;
  while(col != LastCol){
    // chop '-' only columns
    if(col->colDepth == col->colInf[0]){
      if(col != FirstCol)
        col->prev->next = col->next;
      else
        FirstCol = col->next;
      col->next->prev = col->prev;
      ptr=col->frags.down;
      for(i=0; i < col->colDepth; ++i){
        ptr->prev->next = ptr->next;
        ptr->next->prev = ptr->prev;
        ptr = ptr->down;
      }
      ecol = col;
      col = ecol->next;
    }
    else{
      col->frags.el.base = 0;// initialize
      col = col->next;
    }
  }
  for (i=0; i < numFrags; ++i)
    Frags[i].scol->frags.el.base = 1;// TODO
  for(i=0; i < Rows; ++i){
    curEl[i] = NULL;
    curChar[i] = NULL;
  }
  col = FirstCol;
  int be_idx=0;
  int b4p_idx=0;
  if(opt_fastq != 1){
    printf("%%%s\n",chrname);
    printf("#%s\n",buf4printComment);
  }
  while(col != LastCol){
    if(col->frags.el.base){
      for(i=0; i < numFrags; ++i){
        if(col == Frags[i].scol){
          if(curEl[Frags[i].row] == NULL){
            curEl[Frags[i].row] = Frags[i].ofrag;
            curChar[Frags[i].row] = Frags[i].afrag;
          }
          else{
            curEl[Rows] = Frags[i].ofrag;
            curChar[Rows] = Frags[i].afrag;
            ++Rows;
            if(Rows>=MaxDepth*2){
              fprintf(stderr, "souteigai 2 (1) chr %s\n",chrname);
              exit(1);
              //abort();
            }
          }
        }
      }
    }
    int pi=0;
    for(row=0; row < Rows; ++row){
      if(curEl[row] == NULL){
        buf4print[pi].base = ' ';
        buf4print[pi].qv = ' ';
        ++pi;
      }
      else if(curEl[row]->el.base == -1){
        buf4print[pi].base = ' ';
        buf4print[pi].qv = ' ';
        ++pi;
        curEl[row] = NULL;
        curChar[row] = NULL;
      }
      else{
        if(curEl[row]->el.base == 0){
          buf4print[pi].base = '-';
          buf4print[pi].qv = curEl[row]->el.qv;
          ++pi;
        }
        else{
          buf4print[pi].base = (*curChar[row]).base;
          if(buf4print[pi].base == '\0'){
            fprintf(stderr, "base 0\n");
            abort();
          }
          buf4print[pi].qv = (*curChar[row]).qv;
          ++pi;
          ++curChar[row];
        }
        curEl[row] = curEl[row]->next;
      }
    }
    buf4print[pi].base='\0';
    buf4print[pi].qv='\0';
    if(pi == 0){
      fprintf(stderr, "pi 0\n");
      abort();
    }

    if(opt_fastq != 1){
      //printf("%s\n",buf4print);
      if(opt_consensus==1){
        char consensus,tmp2,tmp3;
        col2fqcell(buf4print,&consensus,&tmp2,&tmp3,maximum_ballots,minimum_ballots,distinguishable);
        putchar(consensus);
        putchar('\t');
      }
      int tmp;
      for(tmp = 0; tmp<pi; ++tmp){
        putchar(buf4print[tmp].base);
      }
      if(!opt_vertical){
        putchar('\t');
        for(tmp = 0; tmp<pi; ++tmp){
          putchar(buf4print[tmp].qv);
        }
      }
      putchar('\n');
    }
    else{
      col2fqcell(buf4print, buf4printSeq+b4p_idx, buf4printDepth+b4p_idx, buf4printQual+b4p_idx, maximum_ballots, minimum_ballots,distinguishable);
      if(buf4print[0].base != ' '){
        base_exists[be_idx]=1;
      }
      else{
        base_exists[be_idx]=0;
      }
      if(buf4printSeq[b4p_idx] != ' '){
        ++b4p_idx;
      }
      else{
        // The base of the base read was '-', and col2fqcell did not elect.
        // We ignore this '-'.
      }
      ++be_idx;
    }
    col = col->next;
  }
  buf4printSeq[b4p_idx] = '\0';
  buf4printDepth[b4p_idx] = '\0';
  buf4printQual[b4p_idx] = '\0';
  base_exists[be_idx] = '\0';

  if(opt_fastq != 1){
    printf("\n");
  }
  else{
    print_dfq(chrname,buf4printSeq,buf4printDepth,buf4printQual,base_exists,buf4printComment);
  }
  chr_is=0;
  comment_is=0;
  return;
}

#define DEL 1
#define SUB 0
#define INS 2

double  *mat;
char  *bmat;
int  *bmatPtr;
unsigned long long BmatSize;
int  *shift;

/* reAlign realigns the fnum fragment against the alignment of the rest of the fragments */
// fnum: fragment number = id of the fragment (given in readFrags)
int reAlign (int fnum, ElPool*elpool,colPool*colpool){
  int  i, j, m, n;
  char    *cptr;
  base_t * fel;
  double  min, dval, sval, ival;
  double  *lval, *tval;
  frag    *pf;// pointer of fragment? presenting fragment?
  int    mlen, max, mark;
  colStruc  *col, *minCol, *tcol, *mstart, *mstop;
  fragEl  *ptr, *tptr;
  fragEl  *fptr;

  pf = &Frags[fnum];
  mark = 0;

  /* Strip fragment from structure */
  col = pf->scol;// starting column
  for(ptr=pf->ofrag; ptr->el.base != -1; ptr=ptr->next){// starting element to last element
    ptr->up->down = ptr->down;
    ptr->down->up = ptr->up;
    --(col->colDepth);
    --(col->colInf[(int)ptr->el.base]);
    col = col->next;
  }

  mstart = pf->scol;
  mlen = 1+2*BandSize;
  // set start column
  for(n=0; n < BandSize; ++n){// go back to the n previous colum
    if(mstart == FirstCol){// if you reach FirstColumn, then add a new FirstCol
      FirstCol = mstart->prev = getCol(colpool);
      FirstCol->next = mstart;
    }
    mstart = mstart->prev;
  }
  // done

  // set terminals of a starting element
  // ptr: right, tpter: left
  for(j=0, tptr=ptr=pf->ofrag; j < BandSize; ++j){
    // move to next base, not '-'
    ptr = ptr->next;
    while(ptr->el.base == 0){// '-'
      ++mlen;
      ptr = ptr->next;
    }
    tptr = tptr->prev;// there is no insertion, because pf->ofrag is the terminal of the read
  }
  // done
  // mstart ~ tptr
  // but mstart does not have el, so it needs tptr

  // set score matrix to 0.0
  // mat[0] = MaxFragSize
  // mat[0] is always MaxFragSize
  for(j=1; j <= mlen; ++j)
    mat[j] = 0.0;

  // set stop column. stop column is column next to last column
  for(mstop=pf->ecol->next,j=0; mstop!=LastCol && j<BandSize; ++j,mstop=mstop->next);

  // precalculate deltas of the dynamic programming matrix
  for(col = mstart; col != mstop; col = col->next){
    m = col->colDepth - col->colInf[5];
    if(m != 0){
      max = col->colInf[0];
      for(i=1; i < 5; ++i)
        if(max < col->colInf[i])
          max = col->colInf[i];
      min = m;
      for(i=0; i < 5; ++i){
        col->preCalc[i] = MIX*(1.0-(double)col->colInf[i]/min);
        if(col->colInf[i] != max)
          col->preCalc[i] += (1.0-MIX);
      }
    }
    else{
      for(i=0; i < 5; ++i)
        col->preCalc[i] = 0.0;
    }
  }
  // done

  fel = pf->efrag;// (a pointer to) Fragment's ELement? // encoded character array // this array does not include '-'
  for(i = 1; i <= pf->len; ++i, ++fel){// pf->len is the length of pf->efrag, not including '-' length
    // one '-' lengthens mlen one
    ptr = ptr->next;
    while(ptr->el.base == 0){// '-'
      ++mlen;
      mat[mlen] = MaxFragSize;// out of the band // infinity penalty
      ptr = ptr->next;
    }
    mat[mlen+1] = MaxFragSize;// out of the band // infinity penalty
    shift[i] = 1;// xxxxxxxxxxx      <- i iteration
                 //     xxxxxxxxxxx  <- i+1 iteration
                 // |--|
                 //   `-shift
                 //
                 //  before above, 
                 //  insertions are thought as mismatches like below
                 //
                 //    ss ->  s-s .
                 //    ss ->  s-s .
                 //    ss ->  s-s .
                 // q  \  ->  \   .
                 // q  |  ->   \  .
                 // q   \ ->    \ .

    // one '-' shortens mlen one
    // mstart (col = mstart) is important
    while(tptr->el.base == 0){// '-'
      --mlen;
      ++(shift[i]);// length tptr moved when tptr->el reached next base, not '-'
      tptr = tptr->next;
      mstart = mstart->next;
    }
    tptr = tptr->next;// preparing for the next loop
    col = mstart;
    mstart = mstart->next;// preparing for the next loop
    bmatPtr[i] = mark;// the boundary of bmat
    cptr = &bmat[mark];// backtrack matrix
    mark += mlen;// preparing for the next loop
    if((unsigned long long)mark > BmatSize){
      fprintf(stderr, "too large backtrack matrix. Change the BmatSize value in your realigner\n");
      abort();
    }
    // fill the dynamic programming matrix and get the path (cptr)
    tval = mat;
    lval = &mat[shift[i]];
    // xxtxxxlxxxx      <- i iteration
    //     xxtxxxlxxxx  <- i+1 iteration
    // |--|
    //   `-shift
    for(j=1; j <= mlen && col != LastCol; ++j, col=col->next){
      // (*tval++) = (*(tval++)). pointer moves, not value
      dval = (*tval++) + col->preCalc[0];// Frags[fnum][.] = '-', a deletion for the consensus
      sval = (*lval++) + col->preCalc[(int)(fel->base)];// a match or substitution
      if (fel->base != 5)
        ival = *lval+1.0;// *lval+(1-MIX)*1+MIX*1
      else
        ival = *lval;// 'N' means it has no penalty
      if(sval <= dval && sval <= ival){
        mat[j] = sval;
        *cptr = SUB;
      }
      else if(dval <= ival){
        mat[j] = dval;
        *cptr = DEL;
      }
      else{
        mat[j] = ival;
        *cptr = INS;
      }
      ++cptr;
    }
  }

  cptr = &bmat[bmatPtr[pf->len]];// cptr is 0-origin
  for(n=1,col=mstart->prev;col!=pf->ecol && n<=mlen; ++n,col=col->next);// n>bandwidth, because tptr+bandwidth <= pf->ecol

  min = mat[n];
  minCol = col;
  cptr = &bmat[bmatPtr[pf->len]];
  j = n+1;
  col = minCol->next;
  tcol = minCol->prev;// temporary column? terminal column?
  for(i=n-1; i > 0; --i){
    if(j <= mlen && col != LastCol){
      if(mat[j] < min || (mat[j] == min && cptr[n-1]==DEL)){// score is the strict minimum or (minumum & not a terminal deletion
        n = j;
        min = mat[j];
        minCol = col;
      }
      ++j;
      col = col->next;
    }
    // else{ empty mat[j];}
    if (mat[i] < min || (mat[i] == min && cptr[n-1]==DEL)){
      n = i;
      min = mat[i];
      minCol = tcol;
    }
    tcol = tcol->prev;
  }
  // now, mat[n] is the strict minimum or not strict minimum but not a terminal deletion

  // let's traceback!
  ptr = pf->lastEl;
  mlen = j = n-1;// j is an offset from bmatPtr[pf->len]
  i = pf->len;
  fel = &(pf->efrag[pf->len-1]);
  col = minCol;
  while(i > 0){
    if(bmat[bmatPtr[i]+j] == SUB){
      ptr->el.base = m = fel->base;
      ptr->el.qv = fel->qv;
      --fel;
      ++(col->colDepth);
      ++(col->colInf[m]);
      ptr->up = &(col->frags);
      ptr->down = col->frags.down;
      ptr->up->down = ptr;
      ptr->down->up = ptr;
      col = col->prev;
      j = j+shift[i]-1;
      --i;
    }
    else if(bmat[bmatPtr[i]+j] == DEL){
      ptr->el.base = 0;
      if(ptr->prev->el.qv >= '!'){
        if(ptr->next->el.qv >= '!'){
          ptr->el.qv = (ptr->prev->el.qv+ptr->next->el.qv)/2;//((a-k)+(b-k))/2+k = (a+b)/2
        }
        else{
          ptr->el.qv = ptr->prev->el.qv;
        }
      }
      else if(ptr->next->el.qv >= '!'){
        ptr->el.qv = ptr->next->el.qv;
      }
      else{
        fprintf(stderr, "sth buggy 200\n");
        abort();
      }
      ++(col->colDepth);
      ++(col->colInf[0]);
      ptr->up = &(col->frags);
      ptr->down = col->frags.down;
      ptr->up->down = ptr;
      ptr->down->up = ptr;
      col = col->prev;
      --j;
    }
    else{
      tcol = getCol(colpool);
      tcol->prev = col;
      tcol->next = col->next;
      col->next->prev = tcol;
      col->next = tcol;
      ++(tcol->colDepth);
      ptr->el.base = m = fel->base;
      ptr->el.qv = fel->qv;
      --fel;
      ++(tcol->colInf[m]);
      ptr->down = ptr->up = &(tcol->frags);// frags is the fragEl head of the column
      tcol->frags.down = tcol->frags.up = ptr;
      tcol->frags.prev = (fragEl *) tcol;
      fptr = col->frags.down;
      for(n=0; n < col->colDepth; ++n){
        if(fptr->next->el.base != -1){
          ++(tcol->colDepth);
          ++(tcol->colInf[0]);
          tptr = getEl(elpool);
          tptr->prev = fptr;
          tptr->next = fptr->next;
          tptr->next->prev = tptr;
          tptr->prev->next = tptr;
          tptr->el.base = 0;
          tptr->el.qv = (fptr->el.qv+fptr->next->el.qv)/2;
          tptr->up = &(tcol->frags);
          tptr->down = tcol->frags.down;
          tptr->up->down = tptr;
          tptr->down->up = tptr;
        }
        fptr = fptr->down;
      }
      j = j+shift[i];
      --i;
    }
    if(ptr == pf->ofrag && i > 0){
      pf->ofrag = getEl(elpool);
      pf->ofrag->prev = ptr->prev;
      ptr->prev->next = pf->ofrag;
      ptr->prev = pf->ofrag;
      pf->ofrag->next = ptr;
    }
    ptr = ptr->prev;
  }
  pf->ofrag = ptr->next;
  while(ptr->el.base != -1){
    ptr->el.base = -1;
    ptr->el.qv = ' ';
    ptr = ptr->prev;
  }
  if(col != NULL)
    pf->scol = col->next;
  else
    pf->scol = FirstCol;
  if (bmat[bmatPtr[pf->len]+mlen] == INS)// when last column is an insertion, the inserted column is the last column(pf->ecol)
    pf->ecol = (colStruc *) pf->lastEl->down->prev;
  else
    pf->ecol = minCol;

  return 0;
}

int opt_not_realign=0;

void useErr(char *name)
{
  fprintf(stderr,"usage: %s <in.vertical>\n",name);
  fprintf(stderr,"\t[-l (max read length): (default: %d)]\n",MaxFragSize);
  fprintf(stderr,"\t[-w (bandwidth): (default: 8)]\n");
  fprintf(stderr,"\t[-m (mix): (1.0-mix)*delta_a+mix*delta_c (default: 0.5)]\n");
  fprintf(stderr,"\t[-q: use quality values if exist in the input]\n");
  fprintf(stderr,"\t[-f: outputs in fastq]\n");
  fprintf(stderr,"\t[-b (the minimum of ballots): decides whether to elect or not (default: 3)]\n");
  fprintf(stderr,"\t[-B (the maximum of ballots): scales quality values (default: 11)]\n");
  fprintf(stderr,"\t[-c: outputs consensus bases in vertical like format]\n");
  fprintf(stderr,"\t[-v: outputs only realigned bases in vertical format]\n");
  fprintf(stderr,"\t[-n: does NOT realign but take consensus]\n");
  fprintf(stderr,"\t[-d (0.5-1.0 default:0.7): when voting, if(1st/(1st+2nd) >= d) then change a base else does not]\n");
  exit(1);
}

int main(int argc, char **argv){
  ElPool elpool;
  colPool colpool;
  int hitnum=0;
  {
    int r;
    while((r=getopt(argc,argv,"w:m:fqb:B:cnd:vl:")) != -1){
      switch(r){
        case 'l':
          MaxFragSize = atoi(optarg);
          if(MaxFragSize < 1){
            fprintf(stderr,"Illegal read length: %d\n",MaxFragSize);
            useErr(argv[0]);
          }
          hitnum+=2;
          break;
        case 'w':
          BandSize = atoi(optarg);
          if(BandSize < 1){
            fprintf(stderr,"Illegal band size\n");
            useErr(argv[0]);
          }
          hitnum+=2;
          break;
        case 'm':
          MIX=atof(optarg);
          if(MIX < 0.0 || MIX > 1.0){
            fprintf(stderr, "0.0<=mix<=1.0\n");
            useErr(argv[0]);
          }
          hitnum+=2;
          break;
        case 'f':
          opt_fastq=1;
          hitnum+=1;
          break;
        case 'q':
          opt_qvs=1;
          hitnum+=1;
          break;
        case 'b':
          minimum_ballots=atoi(optarg);
          if(minimum_ballots < 1 || minimum_ballots > 1024){
            fprintf(stderr, "1<=-v(minimum_ballots)<=1024\n");
            useErr(argv[0]);
          }
          hitnum+=2;
          break;
        case 'B':
          maximum_ballots=atoi(optarg);
          if(maximum_ballots < 1 || maximum_ballots > 1024){
            fprintf(stderr, "1<=-v(maximum_ballots)<=1024\n");
            useErr(argv[0]);
          }
          hitnum+=2;
          break;
        case 'c':
          opt_consensus=1;
          hitnum+=1;
          break;
        case 'n':
          opt_not_realign=1;
          hitnum+=1;
          break;
        case 'd':
          distinguishable = atof(optarg);
          if(distinguishable < 0.5 || distinguishable > 1.0){
            fprintf(stderr, "0.5<=distinguishable<=1.0\n");
            useErr(argv[0]);
          }
          hitnum+=2;
          break;
        case 'v':
          opt_vertical=1;
          opt_fastq=0;
          hitnum+=1;
          break;
        default:
          useErr(argv[0]);
          break;
      }
    }
  }
  if(argc != 2+hitnum){
    useErr(argv[0]);
  }
  if(strcmp(argv[1+hitnum],"-") == 0){
    input = stdin;
  }
  else{
    input = fopen(argv[1+hitnum],"r");
    if(input == NULL){
      fprintf(stderr, "cannot open %s\n",argv[1+hitnum]);
      abort();
    }
  }
//  if(opt_not_realign == 1 && opt_fastq != 1){
//    fprintf(stderr, "-n option must be given with -f option\n");
//    return 1;
//  }
//  if(opt_qvs && maximum_ballots<1){
//    fprintf(stderr, "give a value to the -B option\n");
//    useErr(argv[0]);
//  }
  int  i;
  int    numFrags, flag;
  int    score, max, oldScore, n;
  colStruc  *col;
  fragEl  *ep;
  buffer = (char*)malloc(sizeof(char)*BUFSIZE);
  if(buffer == NULL){
    fprintf(stderr, "cannot allocate memory: buffer\n");
    abort();
  }
  qbuf = (char*)malloc(sizeof(char)*BUFSIZE);
  if(qbuf == NULL){
    fprintf(stderr, "cannot allocate memory: qbuf\n");
    abort();
  }
  sbuf = (char*)malloc(sizeof(char)*2*MaxFragSize);
  if(sbuf == NULL){
    fprintf(stderr, "cannot allocate memory: sbuf\n");
    abort();
  }
  chrname = (char*)malloc(sizeof(char)*BUFSIZE);
  if(chrname == NULL){
    fprintf(stderr, "cannot allocate memory: buffer\n");
    abort();
  }
  buf4print = (base_t*)malloc(sizeof(base_t)*MaxDepth);
  if(buf4print == NULL){
    fprintf(stderr, "cannot allocate memory: buf4print\n");
    abort();
  }
  buf4printSeq = (char*)malloc(sizeof(char)*MaxFragSize);
  if(buf4printSeq == NULL){
    fprintf(stderr, "cannot allocate memory: buf4printSeq\n");
    abort();
  }
  buf4printDepth = (char*)malloc(sizeof(char)*MaxFragSize);
  if(buf4printDepth == NULL){
    fprintf(stderr, "cannot allocate memory: buf4printDepth\n");
    abort();
  }
  buf4printQual = (char*)malloc(sizeof(char)*MaxFragSize);
  if(buf4printQual == NULL){
    fprintf(stderr, "cannot allocate memory: buf4printQual\n");
    abort();
  }
  buf4printComment = (char*)malloc(sizeof(char)*MaxFragSize);
  if(buf4printComment == NULL){
    fprintf(stderr, "cannot allocate memory: buf4printComment\n");
    abort();
  }
  base_exists = (char*)malloc(sizeof(char)*MaxFragSize);
  if(base_exists == NULL){
    fprintf(stderr, "cannot allocate memory: base_exists\n");
    abort();
  }
  row = (int *) malloc(MaxDepth*sizeof(int));
  if(row == NULL){
    fprintf(stderr, "cannot allocate memory: row\n");
    abort();
  }
  /* initialize matrices needed for calculations */
  BmatSize = (unsigned long long)MaxFragSize;
  BmatSize *=(unsigned long long)(4*BandSize+2);
  //printf("%llu\n",BmatSize);
  bmat = (char *) malloc(BmatSize*sizeof(char));
  if(bmat == NULL){
    fprintf(stderr, "cannot allocate memory: bmat\n");
    abort();
  }
  mat = (double *) malloc(MaxFragSize*sizeof(double));
  if(mat == NULL){
    fprintf(stderr, "cannot allocate memory: mat\n");
    abort();
  }
  mat[0] = MaxFragSize;
  bmatPtr = (int *) malloc(MaxFragSize*sizeof(int));
  if(bmatPtr == NULL){
    fprintf(stderr, "cannot allocate memory: bmatPtr\n");
    abort();
  }
  shift = (int *) malloc(MaxFragSize*sizeof(int));
  if(shift == NULL){
    fprintf(stderr, "cannot allocate memory: shift\n");
    abort();
  }
  Frags = (frag *) malloc(MaxFragNum*sizeof(frag));
  if(Frags == NULL){
    fprintf(stderr, "cannot allocate memory: Frags\n");
    abort();
  }
  curEl = (fragEl **) malloc(sizeof(fragEl *)*MaxDepth*2);
  if(curEl == NULL){
    fprintf(stderr, "cannot allocate memory: curEl\n");
    abort();
  }
  curChar = (base_t **) malloc(sizeof(base_t *)*MaxDepth*2);
  if(curChar == NULL){
    fprintf(stderr, "cannot allocate memory: curChar\n");
    abort();
  }
  {
    unsigned long long MaxIns = 31;
    init_ElPool(&elpool, MaxFragSize); /* should be insufficient and trigger espansion at some point*/
    size_t colsize = (unsigned long long)(MaxFragSize*(1+MaxIns));
    //colsize = 1000000; /* should be insufficient and trigger espansion at some point*/
    init_colPool(&colpool,colsize);
  }

  for (i = 0; i < 128; i++){
    encode[i] = 5;
  }
  encode['-'] = 0;
  encode['a'] = encode['A'] = 1;
  encode['c'] = encode['C'] = 2;
  encode['g'] = encode['G'] = 3;
  encode['t'] = encode['T'] = 4;
  encode['n'] = encode['N'] = 5;

  while(!DoneFlag){
    buf4printComment[0] = '\0';
    numFrags = readFrags(&elpool, &colpool);
    if(numFrags == 0){
      if(!DoneFlag)
        fprintf(stderr, "no frag\n");
      chr_is=0;
      comment_is=0;
      continue;
    }
    maxNumFrags = (numFrags > maxNumFrags) ? numFrags : maxNumFrags;

    /* tack 2 blank columns on end of alignment to allow movement */
    LastCol->next = getCol(&colpool);
    LastCol->next->prev = LastCol;
    LastCol = LastCol->next;
    LastCol->next = getCol(&colpool);
    LastCol->next->prev = LastCol;
    LastCol = LastCol->next;

    score = 0;
    // get the sum of scores
    for(col=FirstCol; col != LastCol; col=col->next){
      if(col->colDepth > 0){
        max = col->colInf[0];// num of deletions
        for(i = 1; i < 5; ++i)
          if(col->colInf[i] > max)
            max = col->colInf[i];
        score += (col->colDepth-max);
      }
    }

    oldScore = score+1;
    flag = 0;
    if(opt_not_realign == 0){
      while(oldScore > score){
        oldScore = score;
        ++flag;
        for(i=0; i < numFrags; i++){
          reAlign(i, &elpool, &colpool);// reAlign each read
        }

        // calculate the score again
        score = 0;
        n = 0;
        for(col=FirstCol; n < BandSize; ++n,col=col->next){
          if(col->colDepth > 0){
            max = col->colInf[0];
            for(i = 1; i < 5; ++i)
              if(col->colInf[i] > max)
                max = col->colInf[i];
              score += (col->colDepth-max);
          }
        }
        for( ; col != LastCol; col=col->next){
          max = col->colInf[0];
          if (col->colDepth == max)  /* if column of blanks, remove */
          {
            col->prev->next = col->next;
            col->next->prev = col->prev;
            ep=col->frags.down;
            for(i=0; i < col->colDepth; ++i){
              ep->prev->next = ep->next;
              ep->next->prev = ep->prev;
              ep = ep->down;
            }
          }
          else{
            for(i = 1; i < 5; ++i)
              if (col->colInf[i] > max)
                max = col->colInf[i];
            score += (col->colDepth-max);
          }
        }
      }
    }

    printAlign(numFrags);
    //fprintf(stderr,"After %d iterations\n",flag);
    NumIts += flag;
    ++NumCons;

    flushElPool(&elpool);
    flushColPool(&colpool);
  }

  free(bmat);
  free(mat);
  free(bmatPtr);
  free(shift);
  for(i=0; i < maxNumFrags; ++i){
    free(Frags[i].afrag);
    free(Frags[i].efrag);
  }
  free(Frags);

  //fprintf(stderr,"%d %d\t", NumIts,NumCons);
  //fprintf(stderr,"Total %d its in %d contigs for ave of %5.2f its/con\n", NumIts,NumCons,(double)NumIts/(double)NumCons);
  free(buffer);
  free(chrname);
  free(buf4print);
  free(buf4printSeq);
  free(buf4printDepth);
  free(buf4printQual);
  free(buf4printComment);
  free(base_exists);
  free(row);
  free(curEl);
  free(curChar);
  free(sbuf);
  free(qbuf);
  if(strcmp(argv[1+hitnum],"-") == 0){
  }
  else{
    fclose(input);
  }
  
  return 0;
}

