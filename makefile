APPNAME = 'sprai'
VERSION = '0.9.9.23'

PREFIX=$(PWD)
COMPILED= \
bfmt72s \
nss2v_v3 \
myrealigner \
m52bfmt7 \


SCRIPTS= \
ca_ikki_v5.pl \
ezez4makefile_v4.pl \
ezez4qsub_vx1.pl \
ezez_vx1.pl \
dumbbell_filter.pl \
fa2fq.pl \
fq2fa.pl \
fq2idfq.pl \
fqfilt.pl \
get_top_20x_fa.pl \
partition_fa.pl \
get_target_fasta_records.pl \
dfq2fq_v2.pl \
extract_fq.pl \
bfmtx2m4.pl \
check_redundancy.pl \
check_circularity.pl \
#mira_ikki.pl \
#sprai_dagcon_v3.py \
#sprai_dagcon.cfg \
#re2cons.pl \
#ezez4qsub_v9.pl \
#ezez_v8.pl \


all: $(COMPILED)

bfmt72s: bfmt72s.c
	$(CC) -Wall -O3 -g -o $@ $<

nss2v_v3: nss2v_v3.c
	$(CC) -Wall -O3 -g -o $@ $<

myrealigner: myrealigner.c
	$(CC) -Wall -O3 -g -o $@ $^

m52bfmt7: m52bfmt7.c
	$(CC) -Wall -O3 -g -o $@ $<


install: $(COMPILED) $(SCRIPTS)
	chmod 766 $^
	cp -p $^ $(PREFIX)/bin/

