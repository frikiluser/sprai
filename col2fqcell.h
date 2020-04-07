
typedef struct base_t{
  char base;
  char qv;
}base_t;

void set_vals(int col_index, int coded_base, int * ballot, char * max_qvs, base_t * col){
  ++ballot[coded_base];
  max_qvs[coded_base] = (max_qvs[coded_base] < (col[col_index].qv-'!')) ? (col[col_index].qv-'!') : max_qvs[coded_base];
}

// for each col in vertical scrolls
void col2fqcell(base_t * col, char * seq, char * depth, char * qual, int maximum_ballots, int minimum_ballots, double distinguishable){
  if(maximum_ballots<2){
    fprintf(stderr, "maximum_ballots must be >= 2. %d was given.\n",maximum_ballots);
    abort();
  }
  if(minimum_ballots<0){
    fprintf(stderr, "minimum_ballots must be >= 0. %d was given.\n",minimum_ballots);
    abort();
  }

  int i;
  int ballot[6];
  char max_qvs[6];
  for(i=0; i<6; ++i){
    ballot[i]=0;
    max_qvs[i]=0;
  }

  int total_ballots = 0;

  for(i=0; col[i].base != '\0'; ++i){
    if(total_ballots <= maximum_ballots){
      switch(col[i].base){
        case ' ':
          break;
        case 'a':
          set_vals(i,1,ballot,max_qvs,col);
          ++total_ballots;
          break;
        case 'A':
          set_vals(i,1,ballot,max_qvs,col);
          ++total_ballots;
          break;
        case 'c':
          set_vals(i,2,ballot,max_qvs,col);
          ++total_ballots;
          break;
        case 'C':
          set_vals(i,2,ballot,max_qvs,col);
          ++total_ballots;
          break;
        case 'g':
          set_vals(i,3,ballot,max_qvs,col);
          ++total_ballots;
          break;
        case 'G':
          set_vals(i,3,ballot,max_qvs,col);
          ++total_ballots;
          break;
        case 't':
          set_vals(i,4,ballot,max_qvs,col);
          ++total_ballots;
          break;
        case 'T':
          set_vals(i,4,ballot,max_qvs,col);
          ++total_ballots;
          break;
        case 'n':
          set_vals(i,5,ballot,max_qvs,col);
          ++total_ballots;
          break;
        case 'N':
          set_vals(i,5,ballot,max_qvs,col);
          ++total_ballots;
          break;
        case '-':
          set_vals(i,0,ballot,max_qvs,col);
          ++total_ballots;
          break;
        default:
          fprintf(stderr, "arienai: %c\n",col[i].base);
          abort();
      //++ballot[(int)encode[(int)buf4print[i]]];
      }
    }
  }

  int number_of_ballots = 0;
  {
    int i;
    for(i=0; i<6; ++i){
      number_of_ballots += ballot[i];
    }
  }
  if(number_of_ballots < 1){
    fprintf(stderr, "sth buggy: consensus\n");
    abort();
  }
  *depth = (number_of_ballots < 93) ? (char)(number_of_ballots+33) : '~';

  // to elect or not
  if(number_of_ballots < minimum_ballots){// do not elect
    // keep
    *seq = col[0].base;
    *qual = col[0].qv;
    if(*seq == '-'){
      *seq = ' ';// this column will not be printed in output fastq
    }
    if(*qual < '!' && *seq != ' '){
      // error
      int i;
      for(i=0; col[i].base != '\0'; ++i){
        fprintf(stderr,"%c",col[i].base);
      }
      fprintf(stderr,"\n");
      for(i=0; col[i].qv != '\0'; ++i){
        fprintf(stderr,"%c",col[i].qv);
      }
      fprintf(stderr,"\n");
      abort();
    }
    return;
  }
  else{ // do elect
    int max_ballot=0;
    for(i=0; i<6; ++i){
      max_ballot = (max_ballot < ballot[i]) ? (ballot[i]) : max_ballot;
    }
    int second_ballot=0;
    {
      for(i=0; i<6; ++i){
        if(ballot[i] < max_ballot && second_ballot < ballot[i]){
          second_ballot = ballot[i];
        }
      }
      int num_top=0;
      for(i=0; i<6; ++i){
        if(ballot[i] == max_ballot){
          ++num_top;
        }
      }
      if(num_top > 1){
        second_ballot = max_ballot;
      }
    }
//    double distinguishable = 0.70;
    double rate = (double)max_ballot/(double)(max_ballot+second_ballot);
    if(rate>=distinguishable){
      for(i=0; i<6; ++i){
        if(ballot[i] == max_ballot){
          *seq = "-ACGTN"[i];
          *qual = max_qvs[i]+'!';
        }
      }
      return;
    }
    else{
      // keep
      *seq = col[0].base;
      *qual = '!';// set to minimum
      //*qual = col[0].qv;
      if(*seq == '-'){
        *seq = ' ';// this column will not be printed in output fastq
      }
      if(*qual < '!' && *seq != ' '){
        fprintf(stderr,"sth strange base\n");
        int i;
        for(i=0; col[i].base != '\0'; ++i){
          fprintf(stderr,"%c",col[i].base);
        }
        fprintf(stderr,"\n");
        for(i=0; col[i].qv != '\0'; ++i){
          fprintf(stderr,"%c",col[i].qv);
        }
        fprintf(stderr,"\n");
        abort();
      }
      return;
    }
  }

  fprintf(stderr,"never come here\n");
  abort();
  return;
}

/*
void col2fqcell_bak(base_t * col, char * seq, char * depth, char * qual, int maximum_ballots, int minimum_ballots){
  if(maximum_ballots<2){
    fprintf(stderr, "maximum_ballots must be >= 2. %d was given.\n",maximum_ballots);
    abort();
  }
  if(minimum_ballots<0){
    fprintf(stderr, "minimum_ballots must be >= 0. %d was given.\n",minimum_ballots);
    abort();
  }

  int i;
  int ballot[6];
  int sum_qvs[6];
  char max_qvs[6];
  //int pi = strlen(col);
  for(i=0; i<6; ++i){
    ballot[i]=0;
    sum_qvs[i]=0;
    max_qvs[i]=0;
  }
  void set_vals(int col_index, int coded_base){
    ++ballot[coded_base];
    sum_qvs[coded_base] += (int)(col[col_index].qv-'!');
    max_qvs[coded_base] = (max_qvs[coded_base] < (col[col_index].qv-'!')) ? (col[col_index].qv-'!') : max_qvs[coded_base];
  }

  for(i=0; col[i].base != '\0'; ++i){
    switch(col[i].base){
      case 'a':
        set_vals(i,1);
        break;
      case 'A':
        set_vals(i,1);
        break;
      case 'c':
        set_vals(i,2);
        break;
      case 'C':
        set_vals(i,2);
        break;
      case 'g':
        set_vals(i,3);
        break;
      case 'G':
        set_vals(i,3);
        break;
      case 't':
        set_vals(i,4);
        break;
      case 'T':
        set_vals(i,4);
        break;
      case 'n':
        set_vals(i,5);
        break;
      case 'N':
        set_vals(i,5);
        break;
      case '-':
        set_vals(i,0);
        break;
      case ' ':
        break;
      default:
        fprintf(stderr, "arienai: %c\n",col[i].base);
        abort();
    //++ballot[(int)encode[(int)buf4print[i]]];
    }
  }

  int number_of_ballots = 0;
  {
    int i;
    for(i=0; i<6; ++i){
      number_of_ballots += ballot[i];
    }
  }
  if(number_of_ballots < 1){
    fprintf(stderr, "sth buggy: consensus\n");
    abort();
  }
  *depth = (number_of_ballots < 93) ? (char)(number_of_ballots+33) : '~';

  // to elect or not
  if(number_of_ballots < minimum_ballots){
    // do not change the base and qv of the base read.
    *seq = col[0].base;
    // *depth = (char)(1+33);
    *qual = col[0].qv;
    if(*seq == '-'){
      *seq = ' ';// this column will not be printed in output fastq
    }
    if(*qual < '!' && *seq != ' '){
      fprintf(stderr,"kita #%c#, #%c#\n",*seq,*qual);
      int i;
      for(i=0; col[i].base != '\0'; ++i){
        fprintf(stderr,"%c",col[i].base);
      }
      fprintf(stderr,"\n");
      for(i=0; col[i].qv != '\0'; ++i){
        fprintf(stderr,"%c",col[i].qv);
      }
      fprintf(stderr,"\n");
      //abort();
    }
    return;
  }
  else{
    // suppose 1-p ~ 1
    // pi(p(~x)) * pi(1-p(x)) ~ pi(p(~x))
    // {x| min{sum(qv of ~x)}}
    //    = {x| min{sum(qvall) - sum(qv of x)}}
    //    = {x| max{sum(qv of x)}}
    int maxsumqv=0;
    for(i=0; i<6; ++i){
      maxsumqv = (maxsumqv < sum_qvs[i]) ? (sum_qvs[i]) : maxsumqv;
    }
    int second_sum_qv=0;
    for(i=0; i<6; ++i){
      if(sum_qvs[i] < maxsumqv){
        second_sum_qv = (second_sum_qv < sum_qvs[i]) ? sum_qvs[i] : second_sum_qv;
      }
    }

    int num_top=0;
    for(i=0; i<6; ++i){
      if(sum_qvs[i] == maxsumqv){
        ++num_top;
      }
    }

    if(num_top > 1){
      second_sum_qv = maxsumqv;
    }

    for(i=0; i<6; ++i){
      if(sum_qvs[i] == maxsumqv){
        *seq = "-ACGTN"[i];
        // *seq = Decode[i];
        int q_cand = (maxsumqv-second_sum_qv);
        //int q_cand = (maxsumqv-second_sum_qv)/(maximum_ballots-1);// 1/(maximum_ballots-1) is a scaling factor
        if(q_cand > 1){
          if(q_cand < 93){
            *qual = (char)(q_cand+33);
          }
          else{
            *qual = '~';
          }
        }
        else{
          *qual = (char)(1+33);
        }
      }
    }
    if(num_top==1){
      return;
    }
    else{
      // do not change the base (but change qv)
      *seq = col[0].base;
      *qual = (char)(1+33);// (0-0)/(maximum_ballots-1) = 0; -> 1
      return;
    }
  }

  fprintf(stderr,"never come here\n");
  abort();
  return;
}
*/
