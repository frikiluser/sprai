#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define LAUX 32

int LSEQ = 4096;
int LBUF;/* = (3*LSEQ); */
int NUM_ALIGNMENTS = 512;

int maxeval = 1.0;

typedef struct cigar_t{
  int num;
  char sym;
}cigar_t;

typedef struct cmaux_t{
  char sym;
  int num;
}cmaux_t;

typedef struct sam_t{
  char * qname;
  int flag;
  char * rname;
  int pos;
  int mapq;
  struct cigar_t * cigar;
  int cigar_length;
  int cigar_capacity;
  char * rnext;
  int pnext;
  int tlen;
  char * seq;
  char * qual;
  double as;/* bit score */
  double ev;/* e-value */
  double pident;/* percentage of identical matches */
}sam_t;

typedef struct region{
  int left;
  int right;
}region_t;

typedef struct blast{
  int sstart,send,qstart,qend;
  double bitscore,evalue,pident;
}blast_t;

int opt_invsam=0;
int opt_unique_ref=0;
int opt_unique_que=0;
double valid_evalue=-1.0;
int opt_chimeric_filter=0;
int opt_dummy_qv=0;
int opt_unique_bfmt7=0;
int opt_m5=0;
int opt_graph=0;
int opt_top_hit_only=0;
int opt_pctidy=0;
int opt_repeat=0;

void print_sam(sam_t * s){
  printf("%s",s->qname);
  printf("\t");
  printf("%d",s->flag);
  printf("\t");
  printf("%s",s->rname);
  printf("\t");
  printf("%d",s->pos);
  printf("\t");
  printf("%d",s->mapq);
  printf("\t");
  if(opt_graph){
    printf("*");
  }
  else{
    int i;
    for(i=0; i<s->cigar_length; ++i){
      printf("%d%c",s->cigar[i].num,s->cigar[i].sym);
    }
  }
  printf("\t");
  printf("%s",s->rnext);
  printf("\t");
  printf("%d",s->pnext);
  printf("\t");
  printf("%d",s->tlen);
  printf("\t");
  if(opt_graph){
    printf("*");
  }
  else{
    printf("%s",s->seq);
  }
  printf("\t");
  if(opt_graph){
    printf("*");
  }
  else{
    printf("%s",s->qual);
  }
  printf("\t");
  printf("AS:i:%d",(int)(s->as+0.499));/* already rounded? */
  printf("\t");
  if(s->ev > 0.0){
    printf("EV:Z:%1.0e",s->ev);
  }
  else{
    printf("EV:Z:%1.1f",s->ev);
  }
  printf("\t");
  printf("PI:Z:%2.2f",s->pident);
  printf("\n");
  return;
}

void init_sam(sam_t * s, size_t LBUF, int LSEQ){
 /* Because the string element are read with sscanf from a string
 * of at most LBUF long, each string should have a capcity of LBUF.
 * This structure is not allocated for many instance.
 * LSEQ is used to determine the size of cigar, which is not
 * just a string.
 */
  s->qname = (char*)malloc(LBUF);
  if(s->qname == NULL){
    fprintf(stderr,"cannot allocate memory: qname\n");
    abort();
  }
  s->rname = (char*)malloc(LBUF);
  if(s->rname == NULL){
    fprintf(stderr,"cannot allocate memory: rname\n");
    abort();
  }
  s->cigar = (cigar_t*)malloc(sizeof(cigar_t)*LSEQ);
  if(s->cigar == NULL){
    fprintf(stderr,"cannot allocate memory: cigar\n");
    abort();
  }
  s->cigar_capacity = LSEQ;
  s->rnext = (char*)malloc(LBUF);
  if(s->rnext == NULL){
    fprintf(stderr,"cannot allocate memory: rnext\n");
    abort();
  }
  s->seq = (char*)malloc(LBUF);
  if(s->seq == NULL){
    fprintf(stderr,"cannot allocate memory: seq\n");
    abort();
  }
  s->qual = (char*)malloc(LBUF);
  if(s->qual == NULL){
    fprintf(stderr,"cannot allocate memory: qual\n");
    abort();
  }
  return;
}

sam_t * realloc_sam(sam_t * s, size_t LBUF, int LSEQ){
 /* Because the string element are read with sscanf from a string
 * of at most LBUF long, each string should have a capcity of LBUF.
 * This structure is not allocated for many instance.
 * LSEQ is used to determine the size of cigar, which is not
 * just a string.
 */
  void *tmp_p;
  tmp_p = realloc(s->qname, LBUF);
  if(tmp_p == NULL){
    fprintf(stderr,"cannot reallocate memory: qname\n");
    fprintf(stderr,"LBUF: %lu, LSEQ: %d\n", (unsigned long) LBUF, LSEQ);
    exit(EXIT_FAILURE);
  }
  s->qname = (char*)tmp_p;

  tmp_p = realloc(s->rname, LBUF);
  if(tmp_p == NULL){
    fprintf(stderr,"cannot reallocate memory: rname\n");
    fprintf(stderr,"LBUF: %lu, LSEQ: %d\n", (unsigned long) LBUF, LSEQ);
    exit(EXIT_FAILURE);
  }
  s->rname = (char*)tmp_p;

  tmp_p = realloc(s->cigar, sizeof(cigar_t)*LSEQ);
  if(tmp_p == NULL){
    fprintf(stderr,"cannot reallocate memory: cigar\n");
    fprintf(stderr,"LBUF: %lu, LSEQ: %d\n", (unsigned long) LBUF, LSEQ);
    exit(EXIT_FAILURE);
  }
  s->cigar = (cigar_t*)tmp_p;
  s->cigar_capacity = LSEQ;

  tmp_p = realloc(s->rnext, LBUF);
  if(tmp_p == NULL){
    fprintf(stderr,"cannot reallocate memory: rnext\n");
    fprintf(stderr,"LBUF: %lu, LSEQ: %d\n", (unsigned long) LBUF, LSEQ);
    exit(EXIT_FAILURE);
  }
  s->rnext = (char*) tmp_p;
   
  tmp_p = realloc(s->seq, LBUF);
  if(tmp_p == NULL){
    fprintf(stderr,"cannot reallocate memory: seq\n");
    fprintf(stderr,"LBUF: %lu, LSEQ: %d\n", (unsigned long) LBUF, LSEQ);
    exit(EXIT_FAILURE);
  }
  s->seq = (char*) tmp_p;

  tmp_p = realloc(s->qual, LBUF);
  if(tmp_p == NULL){
    fprintf(stderr,"cannot reallocate memory: qual\n");
    fprintf(stderr,"LBUF: %lu, LSEQ: %d\n", (unsigned long) LBUF, LSEQ);
    exit(EXIT_FAILURE);
  }
  s->qual = (char*) tmp_p;

  return s;
}

void reset_cigar(sam_t * s){
  s->cigar_length=0;
}

void reset_sam(sam_t * s){
  s->qname[0] = '\0';
  s->rname[0] = '\0';
  s->rnext[0] = '\0';
  s->cigar_length=0;
  s->seq[0] = '\0';
  s->qual[0] = '\0';
  return;
}

void free_sam(sam_t * s){
  free(s->qname);
  free(s->rname);
  free(s->cigar);
  free(s->rnext);
  free(s->seq);
  free(s->qual);
  return;
}

void blast_print_sam(sam_t * sam, cmaux_t * cmaux, blast_t * blast, char** quals);
void aln2cm(sam_t * sam, char * q, char * s, cmaux_t * cmaux);

char fqname[1024];

int header=0;

int main(int argc, char ** argv)
{
  int hitnum=0;
  char * in_blastn;
  char * prev_ref_name;
  char * prev_que_name;
  sam_t sam;
  char * query_name;
  FILE * fp;
  char * buf;
  char * bufq;
  region_t * regions;
  int region_index=0;
  int scanfret;
  blast_t blast;
  char strand[255];
  char ** quals = NULL;
  int n_units=0;
  LBUF = (3*LSEQ);

  {
    int result;
    while((result=getopt(argc,argv,"iuUe:c:dbxgtprl:")) != -1){
      switch(result){
        case 'l':
          LSEQ=atoi(optarg);
          if(LSEQ < 1){
            fprintf(stderr,"LSEQ>=1: %d is given.",LSEQ);
            exit(1);
          }
          hitnum+=2;
          break;
        case 'r':
          opt_repeat=1;
          header=1;
          ++hitnum;
          break;
        case 'p':
          opt_pctidy=1;
          opt_top_hit_only=1;
          header=1;
          ++hitnum;
          break;
        case 't':
          opt_top_hit_only=1;
          ++hitnum;
          break;
        case 'g':
          opt_graph=1;
          ++hitnum;
          break;
        case 'i':
          opt_invsam=1;
          ++hitnum;
          break;
        case 'u':
          opt_unique_ref=1;
          ++hitnum;
          break;
        case 'U':
          opt_unique_que=1;
          ++hitnum;
          break;
        case 'e':
          valid_evalue=atof(optarg);
          if(valid_evalue < 0.0){
            fprintf(stderr,"e>=0.0: %f is given.",valid_evalue);
            abort();
          }
          hitnum+=2;
          break;
        case 'c':
          opt_chimeric_filter=atoi(optarg);
          if(opt_chimeric_filter < 0){
            fprintf(stderr,"c: %d is given.",opt_chimeric_filter);
            fprintf(stderr,"\tmust be >= 0");
            abort();
          }
          hitnum+=2;
          break;
        case 'd':
          opt_dummy_qv=1;
          ++hitnum;
          break;
        case 'b':
          opt_unique_bfmt7=1;
          ++hitnum;
          break;
        case 'x':
          opt_m5=1;
          opt_unique_bfmt7=1;
          ++hitnum;
          break;
        case '?':
          printf("humei\n");
          break;
        default:
          break;
      }
    }
  }
  LBUF = (3*LSEQ);

  if(argc != 2+hitnum){
    char msg[] = "off";
    fprintf(stderr, "USAGE: <this> <in.blastn | - >\n");
    if(opt_invsam){
      strcpy(msg,"on");
    }
    fprintf(stderr, "\t-i: regards blasted queries as references in the output sam (default: %s)\n",msg);
    fprintf(stderr, "\t-U: avoids double voting\n");
    fprintf(stderr, "\t-e <valid_evalue_threshold> : discards input records with more than the threshold\n");
    fprintf(stderr, "\t-d: give dummy qvs to save memory\n");
    fprintf(stderr, "\t-c <trim length> : trim both <trim length> bases of alignment against chimeric reads (default: 0)\n");
    return 1;
  }
  if(opt_unique_ref + opt_unique_que > 1){
    fprintf(stderr, "-u and -U are incompatible\n");
    return 1;
  }
  if(opt_chimeric_filter && opt_graph){
    opt_chimeric_filter = 0;
    fprintf(stderr, "WARNING: because opt_graph was on, opt_chimeric_fileter was forced to be 0.\n");
  }

  in_blastn = argv[1+hitnum];
  prev_ref_name = (char*)malloc(LBUF);
  if(prev_ref_name == NULL){
    fprintf(stderr,"cannot allocate memory: prev_ref_name\n");
    abort();
  }
  prev_ref_name[0] = '\0';
  prev_que_name = (char*)malloc(LBUF);
  if(prev_que_name == NULL){
    fprintf(stderr,"cannot allocate memory: prev_que_name\n");
    abort();
  }
  prev_que_name[0] = '\0';
  init_sam(&sam, LBUF, LSEQ);
  reset_sam(&sam);
  sam.qname[0]='\0';
  sam.mapq = 255;
  sam.rnext[0]='*';
  sam.rnext[1]='\0';
  sam.pnext=0;
  sam.tlen=0;
  query_name = (char*)malloc(LBUF);
  if(query_name == NULL){
    fprintf(stderr,"cannot allocate memory: query_name\n");
    abort();
  }
  if(in_blastn[0] == '-'){
    fp = stdin;
  }
  else{
    fp = fopen(in_blastn,"r");
  }
  if(fp == NULL){
    fprintf(stderr,"cannot open the file %s\n", in_blastn);
    abort();
  }
  buf = (char*)malloc(LBUF);
  if(buf == NULL){
    fprintf(stderr, "cannot allocate memory: buf\n");
    abort();
  }
  bufq = (char*)malloc(LBUF); /* to avoid overflow, the buffer should be as long as the line */
  if(bufq == NULL){
    fprintf(stderr, "cannot allocate memory: bufq\n");
    abort();
  }
  regions = (region_t*)malloc(sizeof(region_t)*NUM_ALIGNMENTS);
  if(regions == NULL){
    fprintf(stderr, "cannot allocate memory: regions\n");
    abort();
  }
  while(fgets(buf,LBUF,fp) != NULL){
    cmaux_t cmaux;
    int line_length;
    line_length = strlen(buf);
    if(line_length >= LBUF - 1){
      /* The line is incompletely read if the buffer 
       * is fully occuppied and not ending with newline */
      while(line_length >= LBUF - 1 && buf[line_length - 1] != '\n'){
        char*newbuf;
        newbuf = (char*)realloc(buf, LBUF*2);
        if(!newbuf){
          fputs("realloc for buf failed!\n", stderr);
          exit(EXIT_FAILURE);
        }
        buf = newbuf;
        fgets(buf+line_length, LBUF + 1, fp);
        line_length = strlen(buf);
        LBUF *= 2;
      }
      {
        char *newbufq = (char*)realloc(bufq, LBUF);
        if(newbufq == NULL){
          fprintf(stderr, "cannot reallocate memory: bufq\n");
          exit(EXIT_FAILURE);
        }
        bufq = newbufq;
        LSEQ = LBUF / 2; 
/* now that the length of the buffer was measured, 
 * the line should contain qseq and sseq of the same length.
 * Thus the sequence should be less than half of the line length */ 
        realloc_sam(&sam, LBUF, LSEQ);
      }
    }
    if(buf[0] == '#'){
      prev_ref_name[0] = '\0';
      prev_que_name[0] = '\0';
      region_index=0;
      n_units=0;
      continue;
    }
    if(opt_top_hit_only){
      if(n_units){
        continue;
      }
    }
    if(!opt_m5){
      /* -outfmt '7 qseqid sstart send sacc qstart qend bitscore evalue pident qseq sseq' */
      scanfret = sscanf(buf,"%s %d %d %s %d %d %lf %lf %lf %s %s", sam.qname, &blast.sstart, &blast.send, sam.rname, &blast.qstart, &blast.qend, &blast.bitscore, &blast.evalue, &blast.pident, sam.seq, sam.qual);/* sam.qual <- sseq */
      if(opt_graph){
        if(scanfret == 9){
          /* set a dummy base */
          sam.seq[0] = 'A';
          sam.seq[1] = '\0';
          sam.qual[0] = 'A';
          sam.qual[1] = '\0';
        }
        else if(scanfret != 11){
          fprintf(stderr, "sth strange: scanfret %d\n",scanfret);
          fprintf(stderr, "buf: %s\n", buf);
          exit(1);
        }
      }
      else if(scanfret != 11){
        fprintf(stderr, "sth strange: scanfret %d\n",scanfret);
        fprintf(stderr, "buf: %s\n", buf);
        exit(2);
      }
    }
    else{
      /* -outfmt '7 qseqid sacc bitscore pident qstart qend sstrand sstart send' */
      scanfret = sscanf(buf,"%s %s %lf %lf %d %d %s %d %d", sam.qname, sam.rname, &blast.bitscore, &blast.pident, &blast.qstart, &blast.qend, strand, &blast.sstart, &blast.send);
      if(scanfret != 9){
        fprintf(stderr, "sth strange: scanfret %d\n",scanfret);
        fprintf(stderr, "buf: %s\n", buf);
        abort();
      }
    }
    if(strcmp(prev_que_name,sam.qname) == 0){
      ++n_units;
      if(n_units>NUM_ALIGNMENTS){
        /* discard this record */
        continue;
      }
    }
    else{
      n_units = 1;
    }
    if(blast.qend-blast.qstart < 0){
      fprintf(stderr, "unexpected blast qstt qend %d %d\n",blast.qstart,blast.qend);
      exit(1);
    }
    if(opt_repeat){
      int sstt = blast.sstart;
      int send = blast.send;
      int qstt = blast.qstart;
      int qend = blast.qend;
      int len = strlen(sam.seq);
      /* qstt <= qend */
      if(strcmp(sam.qname,sam.rname) == 0){
        if((qstt <= sstt && sstt <= qend) || (qstt <= send && send <= qend)){
          /* diagonal */
          continue;
        }
      }
      if(header){
        header = 0;
        printf("%s\n","len");
      }
      printf("%d\n",len);
      continue;
    }
    if(strcmp(prev_que_name,sam.qname) == 0){
      ++n_units;
      if(n_units>NUM_ALIGNMENTS){
        /* discard this record */
        continue;
      }
    }
    else{
      n_units = 1;
    }
    if(valid_evalue >= 0.0){
      if(blast.evalue > valid_evalue){
        continue;
      }
    }
    if(opt_pctidy){
      int len = strlen(sam.seq);
      int i;
      int match=0;
      int mm=0;
      int indel=0;
      double pctidy;
      for(i=0; i<len; ++i){
        if(sam.seq[i] == sam.qual[i]){
          ++match;
        }
        else if(sam.seq[i] == '-' || sam.qual[i] == '-'){
          ++indel;
        }
        else{
          ++mm;
        }
      }
      pctidy = (double)match/(double)len;
      if(header){
        header = 0;
        printf("%s\t%s\t%s\t%s\t%s\n","pctidy","match","len","indel","mm");
      }
      printf("%.9f\t%d\t%d\t%d\t%d\n",pctidy,match,len,indel,mm);
      continue;
    }
    if(opt_unique_ref){
      if(strcmp(prev_ref_name, sam.rname)==0 && strcmp(prev_que_name,sam.qname)==0){
        /* check overlap */
        int i=0;
        int cf=0;
        for(i=0; i<region_index; ++i){
          if(regions[i].left <= blast.sstart && blast.sstart <= regions[i].right){
            cf=1;
            break;
          }
          if(regions[i].left <= blast.send && blast.send <= regions[i].right){
            cf=1;
            break;
          }
        }
        if(cf==1){
          continue;
        }
      }
      else{
        /* another ref */
        region_index=0;
        strcpy(prev_ref_name, sam.rname);
        strcpy(prev_que_name, sam.qname);
      }
    }
    else if(opt_unique_que){
      if(strcmp(prev_que_name, sam.qname)==0 && strcmp(prev_ref_name,sam.rname)==0){
        /* check overlap */
        int i=0;
        int cf=0;
        for(i=0; i<region_index; ++i){
          if(regions[i].left <= blast.qstart && blast.qstart <= regions[i].right){
            cf=1;
            break;
          }
          if(regions[i].left <= blast.qend && blast.qend <= regions[i].right){
            cf=1;
            break;
          }
        }
        if(cf==1){
          continue;
        }
      }
      else{
        /* another que */
        region_index=0;
        strcpy(prev_que_name, sam.qname);
        strcpy(prev_ref_name, sam.rname);
      }
    }
    else{
      region_index=0;
    }

    if(opt_unique_ref){/* blast.qend - blast.qstart >= 0 is always true */
      if(blast.send - blast.sstart < 0){
        /* reverse */
        sam.flag = 0x10;
        sam.pos = blast.send;
        regions[region_index].left=blast.send;
        regions[region_index].right=blast.sstart;
        ++region_index;
      }
      else{
        /* forward */
        sam.flag = 0x0;
        sam.pos = blast.sstart;
        regions[region_index].left=blast.sstart;
        regions[region_index].right=blast.send;
        ++region_index;
      }
    }
    else if(opt_unique_que){/* blast.qend - blast.qstart >= 0 is always true */
      regions[region_index].left=blast.qstart;
      regions[region_index].right=blast.qend;
      ++region_index;
      if(blast.send - blast.sstart < 0){
        /* reverse */
        sam.flag = 0x10;
        sam.pos = blast.send;
      }
      else{
        /* forward */
        sam.flag = 0x0;
        sam.pos = blast.sstart;
      }
    }
    else{/* blast.qend - blast.qstart >= 0 is always true */
      if(blast.send - blast.sstart < 0){
        /* reverse */
        sam.flag = 0x10;
        sam.pos = blast.send;
      }
      else{
        /* forward */
        sam.flag = 0x0;
        sam.pos = blast.sstart;
      }
    }
    if(region_index >= NUM_ALIGNMENTS){
      region_t * tmp;
      NUM_ALIGNMENTS = NUM_ALIGNMENTS * 2;
      tmp = (region_t*)realloc(regions,sizeof(region_t)*NUM_ALIGNMENTS);
      if(tmp == NULL){
        fprintf(stderr, "cannot reallocate memory: regions. NUM_ALIGNMENTS %d\n",NUM_ALIGNMENTS);
        abort();
      }
      else{
        regions = tmp;
      }
    }

    if(opt_unique_bfmt7){
      if(!opt_m5){
        /* -outfmt '7 qseqid sstart send sacc qstart qend bitscore evalue pident qseq sseq' */
        printf("%s\t%d\t%d\t%s\t%d\t%d\t%e\t%e\t%f\t%s\t%s\n", sam.qname, blast.sstart, blast.send, sam.rname, blast.qstart, blast.qend, blast.bitscore, blast.evalue, blast.pident, sam.seq, sam.qual);
      }
      else{
        /* -outfmt '7 qseqid sacc bitscore pident qstart qend sstrand sstart send' */
        printf("%s\t%s\t%e\t%f\t%d\t%d\t%s\t%d\t%d\n", sam.qname, sam.rname, blast.bitscore, blast.pident, blast.qstart, blast.qend, strand, blast.sstart, blast.send);
      }
      continue;
    }

    if(abs(blast.bitscore) > 255.0){
      sam.mapq = 255;
    }
    else{
      sam.mapq = (int)(abs(blast.bitscore)+0.499);
    }
    cmaux.sym='\0'; cmaux.num=0;
    sam.cigar_length=0;

    /* trim terminals agains chimera */
    if(opt_chimeric_filter > 0){
      if(strcmp(sam.qname,sam.rname) != 0){
        int trim = opt_chimeric_filter;
        int len_sseq = strlen(sam.qual);
        int h_trim = trim;
        int t_trim = trim;
        if(len_sseq != strlen(sam.seq)){
          fprintf(stderr, "lengths of que and sbj differ\n");
          abort();
        }
        if(len_sseq <= trim*2){/* sam.qual <- sseq */
          continue;
        }

        while(sam.seq[h_trim] == '-' || sam.qual[h_trim] == '-'){
          ++h_trim;
        }
        while(sam.seq[len_sseq-1-t_trim] == '-' || sam.qual[len_sseq-1-t_trim] == '-'){
          ++t_trim;
        }
        if(h_trim > len_sseq-t_trim){
          continue;
        }

        /* adapt blast.qstart */
        {
          int nhyphen=0;
          int i;
          for(i=0; i<h_trim; ++i){
            if(sam.seq[i] == '-')
              ++nhyphen;
          }
          blast.qstart += (h_trim-nhyphen);
        }

        /* adapt blast.qend */
        {
          int nhyphen=0;
          int i;
          int endpos=strlen(sam.seq)-1;
          for(i=0; i<t_trim; ++i){
            if(sam.seq[endpos-i] == '-')
              ++nhyphen;
          }
          blast.qend -= (t_trim-nhyphen);
        }
        
        /* adapt blast.sstart */
        {
          int nhyphen=0;
          int i;
          for(i=0; i<h_trim; ++i){
            if(sam.qual[i] == '-')
              ++nhyphen;
          }
          if(sam.flag == 0x0){
            blast.sstart += (h_trim-nhyphen);
          }
          else if(sam.flag == 0x10){
            blast.sstart -= (h_trim-nhyphen);
          }
        }

        /* adapt blast.send */
        {
          int nhyphen=0;
          int i;
          int endpos=strlen(sam.qual)-1;
          for(i=0; i<t_trim; ++i){
            if(sam.qual[endpos-i] == '-')
              ++nhyphen;
          }
          if(sam.flag == 0x0){
            blast.send -= (t_trim-nhyphen);
          }
          else if(sam.flag == 0x10){
            blast.send += (t_trim-nhyphen);
          }
        }

        /* adapt sam.pos */
        if(sam.flag == 0x0){
          int nhyphen=0;
          int i;
          for(i=0; i<h_trim; ++i){
            if(sam.qual[i] == '-')
              ++nhyphen;
          }
          sam.pos += (h_trim-nhyphen);
        }
        else if(sam.flag == 0x10){
          int nhyphen=0;
          int i;
          int endpos=strlen(sam.qual)-1;
          for(i=0; i<t_trim; ++i){
            if(sam.qual[endpos-i] == '-')
              ++nhyphen;
          }
          sam.pos += (t_trim-nhyphen);
        }
        else{
          fprintf(stderr, "souteigai flag\n");
          abort();
        }

        strcpy(bufq, &sam.qual[h_trim]);
        bufq[strlen(bufq)-t_trim]='\0';
        strcpy(sam.qual, bufq);
        strcpy(bufq, &sam.seq[h_trim]);
        bufq[strlen(bufq)-t_trim]='\0';
        strcpy(sam.seq, bufq);
      }
      else{
      }
    }

    aln2cm(&sam, sam.seq, sam.qual, &cmaux);

    sam.rnext[0] = '*';
    sam.rnext[1] = '\0';
    sam.pnext=0;
    sam.tlen=0;
    sam.as=abs(blast.bitscore);
    sam.ev=blast.evalue;
    sam.pident=blast.pident;
    blast_print_sam(&sam, &cmaux, &blast, quals);
  }

  if(in_blastn[0] != '-'){
    fclose(fp);
  }
  free_sam(&sam);
  free(query_name);
  free(buf);
  free(bufq);
  free(prev_ref_name);
  free(prev_que_name);
  free(regions);
  return 0;
}

void blast_print_sam(sam_t * sam, cmaux_t * cmaux, blast_t * blast, char** quals){
  if(cmaux->num > 0){
    if(sam->cigar_length >= sam->cigar_capacity){
      fputs("cigar capacity exceeded in blast_print_sam()!\n", stderr);
      exit(EXIT_FAILURE);
    }
    sam->cigar[sam->cigar_length].sym = cmaux->sym;
    sam->cigar[sam->cigar_length].num = cmaux->num;
    ++sam->cigar_length;
  }

  if(opt_invsam == 1){
    char * tmp = sam->rname;
    sam->rname = sam->qname;
    sam->qname = tmp;
    sam->pos = blast->qstart;
    {
      int i;
      int loop;
      for(i=0,loop=sam->cigar_length; i<loop; ++i){
        if(sam->cigar[i].sym == 'D'){
          sam->cigar[i].sym = 'I';
        }
        else if(sam->cigar[i].sym == 'I'){
          sam->cigar[i].sym = 'D';
        }
        else if(sam->cigar[i].sym == 'M'){
          /* do nothing */
        }
        else{
          fprintf(stderr, "souteigai cigar\n");
          abort();
        }
      }
    }
    /* swap query for subject */
    tmp = sam->seq;
    sam->seq = sam->qual;
    sam->qual = tmp;
  }

  {
    int loop=strlen(sam->seq);
    int i,j;
    char * tmp;
    for(i=0,j=0; i<loop; ++i){
      if(sam->seq[i] != '-'){
        sam->qual[j++] = sam->seq[i];
      }
    }
    sam->qual[j]='\0';
    tmp = sam->seq;
    sam->seq = sam->qual;
    sam->qual = tmp;
  }

  /* append quality values */
  sam->qual[0] = '*';
  sam->qual[1] = '\0';
  if(opt_invsam != 1 && sam->flag & 0x10){
    {
      int i;
      int loop;
      int len;
      for(i=0,len=sam->cigar_length,loop=len/2; i<loop; ++i){
        char tmp = sam->cigar[i].sym;
        int tmp2 = sam->cigar[i].num;
        sam->cigar[i].sym = sam->cigar[len-1-i].sym;
        sam->cigar[i].num = sam->cigar[len-1-i].num;
        sam->cigar[len-1-i].sym = tmp;
        sam->cigar[len-1-i].num = tmp2;
      }
      for(i=0,len=strlen(sam->seq),loop=len/2; i<loop; ++i){
        char tmp = sam->seq[i];
        sam->seq[i] = sam->seq[len-1-i];
        sam->seq[len-1-i] = tmp;
      }
      for(i=0,loop=strlen(sam->seq); i<loop; ++i){
        /* complementary nucleotide */
        char t = sam->seq[i];
        switch(t){
          case 'a':
          case 'A':
            t = 'T';
            break;
          case 'c':
          case 'C':
            t = 'G';
            break;
          case 'g':
          case 'G':
            t = 'C';
            break;
          case 't':
          case 'T':
            t = 'A';
            break;
          case 'N':
          case 'n':
            t = 'N';
            break;
          default:
            fprintf(stderr, "souteigai : %c\n",t);
            abort();
            break;
        }
        sam->seq[i] = t;
      }
    }
  }

  /* RNEXT is used for qstt,qend,sstt,send in this program. */
  /* these are trimmed values. */
  /* 1-origin */
  if(opt_invsam){
    sprintf(sam->rnext,"%d,%d,%d,%d",blast->sstart,blast->send,blast->qstart,blast->qend);
  }
  else{
    sprintf(sam->rnext,"%d,%d,%d,%d",blast->qstart,blast->qend,blast->sstart,blast->send);
  }

  print_sam(sam);
  return;
}

void aln2cm(sam_t * sam, char * q, char * s, cmaux_t * cmaux){
  int len = strlen(q);
  int i;
  for(i=0; i<len; ++i){
    int op;
    if(q[i] == '-'){
      op = 'D';
    }
    else if(s[i] == '-'){
      op = 'I';
    }
    else{
      op = 'M';
    }
    if(cmaux->sym == op){
      ++cmaux->num;
    }
    else{
      if(cmaux->num > 0){
        if(sam->cigar_length >= sam->cigar_capacity){
          fputs("cigar capacity exceeded in aln2cm()!\n", stderr);
          exit(EXIT_FAILURE);
        }
        sam->cigar[sam->cigar_length].sym = cmaux->sym;
        sam->cigar[sam->cigar_length].num = cmaux->num;
        ++sam->cigar_length;
      }
      cmaux->sym = op;
      cmaux->num = 1;
    }
  }
  return;
}
