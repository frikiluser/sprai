#!/usr/bin/python3
# encoding: utf-8
# WARNING! Do not edit! http://waf.googlecode.com/git/docs/wafbook/single.html#_obtaining_the_waf_file

import os
from waflib import Utils
from waflib.Tools import ccroot,ar
from waflib.Configure import conf
@conf
def find_sxx(conf):
	v=conf.env
	cc=None
	if v['CXX']:cc=v['CXX']
	elif'CXX'in conf.environ:cc=conf.environ['CXX']
	if not cc:cc=conf.find_program('CC',var='CXX')
	if not cc:cc=conf.find_program('c++',var='CXX')
	if not cc:conf.fatal('Could not find a Sun C++ compiler')
	cc=conf.cmd_to_list(cc)
	try:
		conf.cmd_and_log(cc+['-flags'])
	except Exception:
		conf.fatal('%r is not a Sun compiler'%cc)
	v['CXX']=cc
	v['CXX_NAME']='sun'
@conf
def sxx_common_flags(conf):
	v=conf.env
	v['CXX_SRC_F']=[]
	v['CXX_TGT_F']=['-c','-o']
	if not v['LINK_CXX']:v['LINK_CXX']=v['CXX']
	v['CXXLNK_SRC_F']=[]
	v['CXXLNK_TGT_F']=['-o']
	v['CPPPATH_ST']='-I%s'
	v['DEFINES_ST']='-D%s'
	v['LIB_ST']='-l%s'
	v['LIBPATH_ST']='-L%s'
	v['STLIB_ST']='-l%s'
	v['STLIBPATH_ST']='-L%s'
	v['SONAME_ST']='-Wl,-h,%s'
	v['SHLIB_MARKER']='-Bdynamic'
	v['STLIB_MARKER']='-Bstatic'
	v['cxxprogram_PATTERN']='%s'
	v['CXXFLAGS_cxxshlib']=['-Kpic','-DPIC']
	v['LINKFLAGS_cxxshlib']=['-G']
	v['cxxshlib_PATTERN']='lib%s.so'
	v['LINKFLAGS_cxxstlib']=['-Bstatic']
	v['cxxstlib_PATTERN']='lib%s.a'
def configure(conf):
	conf.find_sxx()
	conf.find_ar()
	conf.sxx_common_flags()
	conf.cxx_load_tools()
	conf.cxx_add_flags()
	conf.link_add_flags()