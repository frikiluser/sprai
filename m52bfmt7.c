#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int LBUF = 1024*1024*1024;
int LNAME = 1024;

typedef struct blast{
  int sstart,send,qstart,qend;
  double bitscore,evalue,pident;
}blast_t;


void print_bfmt7(blast_t * b){
  return;
}

void mem(char * var_name){
  fprintf(stderr,"cannot allocate memory: %s\n",var_name);
  exit(1);
}

int opt_swap_query_and_subject=0;

int main(int argc, char ** argv)
{

  int hitnum=0;
  {
    int result;
    while((result=getopt(argc,argv,"s")) != -1){
      switch(result){
        case 's':
          opt_swap_query_and_subject=1;
          ++hitnum;
          break;
//        case 'c':
//          opt_chimeric_filter=atoi(optarg);
//          if(opt_chimeric_filter < 0){
//            fprintf(stderr,"c: %d is given.",opt_chimeric_filter);
//            fprintf(stderr,"\tmust be >= 0");
//            abort();
//          }
//          hitnum+=2;
//          break;
        case '?':
          printf("humei\n");
          break;
        default:
          break;
      }
    }
  }

  if(argc != 2+hitnum){
    fprintf(stderr, "USAGE: <this> <in.m5 | - >\n");
    fprintf(stderr, "\t-s: swap query and subject names\n");
    return 1;
  }
  char * in_m5 = argv[1+hitnum];

  FILE * fp;
  if(in_m5[0] == '-'){
    fp = stdin;
  }
  else{
    fp = fopen(in_m5,"r");
  }
  if(fp == NULL){
    fprintf(stderr,"cannot open the file %s\n", in_m5);
    abort();
  }

  char * buf = (char*)malloc(sizeof(char)*LBUF);
  if(buf == NULL){
    mem("buf");
    fprintf(stderr, "cannot allocate memory: buf\n");
    exit(1);
  }
  char * qseq = (char*)malloc(sizeof(char)*LBUF);
  if(qseq == NULL){
    mem("qseq");
    fprintf(stderr, "cannot allocate memory: qseq\n");
    exit(1);
  }
  char * sseq = (char*)malloc(sizeof(char)*LBUF);
  if(sseq == NULL){
    mem("sseq");
    fprintf(stderr, "cannot allocate memory: sseq\n");
    exit(1);
  }
  char * ali = (char*)malloc(sizeof(char)*LBUF);
  if(ali == NULL){
    mem("ali");
    fprintf(stderr, "cannot allocate memory: ali\n");
    exit(1);
  }
  char * qname = (char*)malloc(sizeof(char)*LNAME);
  if(qname==NULL){
    mem("qname");
    fprintf(stderr,"cannot allocate memory: qname\n");
    exit(1);
  }
  char * sname = (char*)malloc(sizeof(char)*LNAME);
  if(sname==NULL){
    mem("sname");
    fprintf(stderr,"cannot allocate memory: sname\n");
    exit(1);
  }
  int qlen,qstt,qend,slen,sstt,send,score,nMatch,nMismatch,nIns,nDel,mapqv;
  char qstrand,sstrand;

  int scanfret;

  while(fgets(buf,LBUF,fp) != NULL){
    // m5:
    // qname qlen qstt(0-origin) qend(1-origin) qstrand sname slen sstt(0-origin) send(1-origin) sstrand score nMatch nMismatch nIns nDel mapqv qseq ali sseq
    // [stt,end)
    scanfret = sscanf(buf,"%s %d %d %d %s %s %d %d %d %s %d %d %d %d %d %d %s %s %s", qname,&qlen,&qstt,&qend,&qstrand,sname,&slen,&sstt,&send,&sstrand,&score,&nMatch,&nMismatch,&nIns,&nDel,&mapqv,qseq,ali,sseq);
    if(scanfret != 19){
      fprintf(stderr, "sth strange: scanfret %d\n",scanfret);
      exit(1);
    }
    // -outfmt '7 qseqid sstart send sacc qstart qend bitscore evalue pident qseq sseq'
    double dummy_evalue=0.0;
    int qnamelen = strlen(qname);
    {
      int i;
      for(i=qnamelen-1; i>=0; --i){
        if(qname[i] == '/'){
          qname[i] = '\0';
          break;
        }
      }
    }
    if(sstrand == '-'){// qstrand is always '+'
      int tmp = send;
      send = sstt+1;// 1-originize
      sstt = tmp-1;// 0-originize
    }
    printf("%s %d %d %s %d %d %d %f %f %s %s\n",qname,sstt+1,send,sname,qstt+1,qend,score,dummy_evalue,(double)(nMatch)/(double)(nMatch+nMismatch+nIns+nDel),qseq,sseq);
  }
  free(buf);
  free(qseq);
  free(sseq);
  free(ali);
  free(qname);
  free(sname);
  return 0;
}
