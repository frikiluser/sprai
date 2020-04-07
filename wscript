APPNAME = 'sprai'
VERSION = '0.9.9.23'

srcdir = '.'
blddir = 'build'

def options(opt):
	opt.load('compiler_c perl')

def configure(conf):
	conf.load('compiler_c perl')

#	if conf.check_cc(lib = 'z'):
#		pass
#	else:
#		conf.fatal("zlib does not exist.")

	conf.check_perl_version()

#	if conf.check_perl_module('Compress::Zlib') is None:
#		conf.fatal("Perl module 'Compress::Zlib' does not exist.")
#	else:
#		pass

	conf.check_perl_module('Statistics::Descriptive')

def build(bld):
	bld.install_files(
        '${PREFIX}/bin', [
        'ca_ikki_v5.pl',
        'dumbbell_filter.pl',
#        'ezez4qsub_v9.pl',
#        'ezez4qsub_v8_iter.pl',
        'ezez4qsub_vx1.pl',
#        'ezez_v8.pl',
#        'ezez_v7_iter.pl',
        'ezez_vx1.pl',
        'fa2fq.pl',
        'fq2fa.pl',
        'fq2idfq.pl',
        'fqfilt.pl',
        'get_top_20x_fa.pl',
#        'mira_ikki.pl',
        'partition_fa.pl',
        'ezez4makefile_v4.pl',
        'get_target_fasta_records.pl',
        'dfq2fq_v2.pl',
        'extract_fq.pl',
        'check_redundancy.pl',
        'check_circularity.pl',
        'bfmtx2m4.pl',
#        'sprai_dagcon_v3.py',
#        'sprai_dagcon.cfg',
#        're2cons.pl',
#        'count_chars.pl'
        ], chmod=0o755)

	bld.program(
		source='bfmt72s.c',
		target='bfmt72s',

		includes='.',

#		lib=['z'],

		install_path = '${PREFIX}/bin',
		cflags       = ['-O3', '-Wall'],
		dflags       = ['-g'],
	)

	bld.program(
		source='nss2v_v3.c',
		target='nss2v_v3',

		includes='.',

		install_path = '${PREFIX}/bin',
		cflags       = ['-O3', '-Wall'],
		dflags       = ['-g'],
	)

	bld.program(
		source='myrealigner.c',
		target='myrealigner',

		includes='.',

		install_path = '${PREFIX}/bin',
		cflags       = ['-O3', '-Wall'],
		dflags       = ['-g'],
	)

	bld.program(
		source='m52bfmt7.c',
		target='m52bfmt7',

		includes='.',

		install_path = '${PREFIX}/bin',
		cflags       = ['-O3', '-Wall'],
		dflags       = ['-g'],
	)

def dist(ctx):
	ctx.algo = 'tar.gz'
	ctx.files = ctx.path.ant_glob([
        'ca_ikki_v5.pl',
        'dfq2fq_v2.pl',
        'dumbbell_filter.pl',
#        'ezez4qsub_v9.pl',
#        'ezez4qsub_v8_iter.pl',
        'ezez4qsub_vx1.pl',
#        'ezez_v8.pl',
#        'ezez_v7_iter.pl',
        'ezez_vx1.pl',
        'fa2fq.pl',
        'fq2fa.pl',
        'fq2idfq.pl',
        'fqfilt.pl',
        'get_top_20x_fa.pl',
#        'mira_ikki.pl',
        'partition_fa.pl',
        'ezez4makefile_v4.pl',
        'get_target_fasta_records.pl',
        'doc/_build/html/**',
        'configure',
        'pbasm.spec',
#        'myasm.spec',
        'ec.spec',
#        'ec_iter.spec',
        'bfmt72s.c',
        'col2fqcell.h',
        'LICENSE.txt',
        'ChangeLog.txt',
        'myrealigner.c',
        'nss2v_v3.c',
        'waf',
        'wscript',
        'm52bfmt7.c',
        'bfmtx2m4.pl',
        'extract_fq.pl',
        'check_redundancy.pl',
        'check_circularity.pl',
#        'count_chars.pl',
#        'sprai_dagcon_v3.py',
#        'sprai_dagcon.cfg',
#        're2cons.pl',
        'makefile'
        ])
